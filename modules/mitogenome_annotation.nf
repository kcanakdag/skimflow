// Annotation of GetOrganelle mitogenome FASTA files.
//
// GetOrganelle emits an assembled mitochondrial FASTA. MitoZ and MITOS2
// consume that FASTA for gene annotation, so this module runs downstream of
// MITO.out.fasta.

process MITOS2_DB {
    tag "${refseqver}"
    container params.mitos2_container
    publishDir "${params.outdir}/mitogenome_annotation", mode: 'copy'

    cpus 1
    memory '2 GB'
    time '1h'

    input:
    val refseqver

    output:
    path 'mitos2_refdata', emit: refdir

    script:
    """
    set -euo pipefail

    mkdir -p mitos2_refdata
    python3 -c "import urllib.request; urllib.request.urlretrieve('https://zenodo.org/records/4284483/files/${refseqver}.tar.bz2?download=1', '${refseqver}.tar.bz2')"
    tar -xjf ${refseqver}.tar.bz2 -C mitos2_refdata
    """
}

process MITOZ_ANNOTATE {
    tag "${meta.id}"
    container params.mitoz_container
    publishDir "${params.outdir}/mitogenome_annotation/mitoz", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 4
    memory '8 GB'
    time '2h'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.mitoz"),            emit: dir
    tuple val(meta), path("${meta.id}.mitoz.log.txt"),    emit: log
    tuple val(meta), path("${meta.id}.mitoz.status.txt"), emit: status
    tuple val(meta), path("${meta.id}.mitoz.gbk"),        emit: genbank, optional: true
    tuple val(meta), path("${meta.id}.mitoz_mqc.tsv"),    emit: summary

    script:
    """
    set -euo pipefail

    topo="${params.mitogenome_topology}"
    case "\$topo" in
        auto|linear|circular) ;;
        *) echo "mitogenome_topology must be one of: auto, linear, circular" >&2; exit 2 ;;
    esac

    awk -v global_topo="\$topo" '
        /^>/ {
            n++
            topo = global_topo
            if (topo == "auto") {
                line = tolower(\$0)
                topo = (line ~ /circular/) ? "circular" : "linear"
            }
            printf(">mito%d topology=%s\\n", n, topo)
            next
        }
        { print }
    ' ${fasta} > ${meta.id}.mitoz.input.fasta

    mkdir -p ${meta.id}.mitoz

    set +e
    mitoz annotate \\
        --genetic_code ${params.mitoz_genetic_code} \\
        --clade ${params.mitoz_clade} \\
        --outprefix ${meta.id} \\
        --thread_number ${task.cpus} \\
        --fastafiles ${meta.id}.mitoz.input.fasta \\
        > ${meta.id}.mitoz.log.txt 2>&1
    status=\$?
    set -e

    printf "%s\\n" "\$status" > ${meta.id}.mitoz.status.txt

    [ -f mitoz.log ] && cp mitoz.log ${meta.id}.mitoz/mitoz.log
    [ -d mt_annotation ] && cp -a mt_annotation ${meta.id}.mitoz/
    find . -maxdepth 1 -type d -name "${meta.id}*.result" -exec cp -a {} ${meta.id}.mitoz/ \\;

    gb=\$(find ${meta.id}.mitoz -type f \\( -name "*.gbf" -o -name "*.gb" -o -name "*.gbk" -o -name "*.genbank" \\) | head -1 || true)
    if [ -n "\$gb" ]; then
        cp "\$gb" ${meta.id}.mitoz.gbk
        gbk="yes"
        cds=\$(grep -Ec '^     CDS[[:space:]]' "\$gb" || true)
        trna=\$(grep -Ec '^     tRNA[[:space:]]' "\$gb" || true)
        rrna=\$(grep -Ec '^     rRNA[[:space:]]' "\$gb" || true)
    else
        gbk="no"
        cds=0
        trna=0
        rrna=0
    fi

    features=\$((cds + trna + rrna))
    if [ "\$status" -ne 0 ]; then
        status_label="exit_\$status"
    elif [ "\$gbk" = "yes" ]; then
        status_label="ok"
    else
        status_label="no_genbank"
    fi

cat > ${meta.id}.mitoz_mqc.tsv <<EOF
# id: "mitoz_annotation"
# section_name: "MitoZ annotation"
# description: "MitoZ annotation summary for GetOrganelle mitogenome FASTA files."
# plot_type: "table"
Sample	Status	CDS	tRNA	rRNA	Features	GenBank
${meta.id}	\$status_label	\$cds	\$trna	\$rrna	\$features	\$gbk
EOF

    exit 0
    """
}

process MITOS2_ANNOTATE {
    tag "${meta.id}"
    container params.mitos2_container
    publishDir "${params.outdir}/mitogenome_annotation/mitos2", mode: 'copy', saveAs: { fn -> "${meta.id}/${fn}" }

    cpus 4
    memory '8 GB'
    time '2h'

    input:
    tuple val(meta), path(fasta)
    path refdir

    output:
    tuple val(meta), path("${meta.id}.mitos2"),            emit: dir
    tuple val(meta), path("${meta.id}.mitos2.log.txt"),    emit: log
    tuple val(meta), path("${meta.id}.mitos2.status.txt"), emit: status
    tuple val(meta), path("${meta.id}.mitos2.gff"),        emit: gff, optional: true
    tuple val(meta), path("${meta.id}.mitos2.bed"),        emit: bed, optional: true
    tuple val(meta), path("${meta.id}.mitos2_mqc.tsv"),    emit: summary

    script:
    """
    set -euo pipefail

    topo="${params.mitogenome_topology}"
    case "\$topo" in
        auto)
            if grep -iq 'circular' ${fasta}; then
                linear_arg=""
            else
                linear_arg="--linear"
            fi
            ;;
        linear)
            linear_arg="--linear"
            ;;
        circular)
            linear_arg=""
            ;;
        *)
            echo "mitogenome_topology must be one of: auto, linear, circular" >&2
            exit 2
            ;;
    esac

    mkdir -p ${meta.id}.mitos2

    set +e
    runmitos.py \\
        --input ${fasta} \\
        --code ${params.mitos2_genetic_code} \\
        --outdir ${meta.id}.mitos2 \\
        --refdir ${refdir} \\
        --refseqver ${params.mitos2_refseqver} \\
        \$linear_arg \\
        ${params.mitos2_extra_args} \\
        > ${meta.id}.mitos2.log.txt 2>&1
    status=\$?
    set -e

    printf "%s\\n" "\$status" > ${meta.id}.mitos2.status.txt

    gff=\$(find ${meta.id}.mitos2 -type f \\( -name "*.gff" -o -name "*.gff3" \\) | head -1 || true)
    bed=\$(find ${meta.id}.mitos2 -type f -name "*.bed" | head -1 || true)
    if [ -n "\$gff" ]; then
        cp "\$gff" ${meta.id}.mitos2.gff
        gff_out="yes"
        cds=\$(awk 'BEGIN{c=0} !/^#/ && tolower(\$3)=="cds"{c++} END{print c+0}' "\$gff")
        trna=\$(awk 'BEGIN{c=0} !/^#/ && tolower(\$3) ~ /trna/{c++} END{print c+0}' "\$gff")
        rrna=\$(awk 'BEGIN{c=0} !/^#/ && tolower(\$3) ~ /rrna/{c++} END{print c+0}' "\$gff")
        features=\$(awk 'BEGIN{c=0} !/^#/ && NF>=3{c++} END{print c+0}' "\$gff")
    else
        gff_out="no"
        cds=0
        trna=0
        rrna=0
        features=0
    fi
    if [ -n "\$bed" ]; then
        cp "\$bed" ${meta.id}.mitos2.bed
    fi

    if [ "\$status" -ne 0 ]; then
        status_label="exit_\$status"
    elif [ "\$gff_out" = "yes" ]; then
        status_label="ok"
    else
        status_label="no_gff"
    fi

cat > ${meta.id}.mitos2_mqc.tsv <<EOF
# id: "mitos2_annotation"
# section_name: "MITOS2 annotation"
# description: "MITOS2 annotation summary for GetOrganelle mitogenome FASTA files."
# plot_type: "table"
Sample	Status	CDS	tRNA	rRNA	Features	GFF
${meta.id}	\$status_label	\$cds	\$trna	\$rrna	\$features	\$gff_out
EOF

    exit 0
    """
}

workflow MITOGENOME_ANNOTATION {
    take:
    mito_fasta_ch

    main:
    MITOZ_ANNOTATE(mito_fasta_ch)

    if (params.mitos2_refdir) {
        mitos2_refdir_ch = Channel.value(file(params.mitos2_refdir, checkIfExists: true))
    } else {
        db_req_ch = mito_fasta_ch.map { meta, fasta -> params.mitos2_refseqver }.first()
        MITOS2_DB(db_req_ch)
        mitos2_refdir_ch = MITOS2_DB.out.refdir.collect(flat: false).map { it[0] }
    }
    MITOS2_ANNOTATE(mito_fasta_ch, mitos2_refdir_ch)

    emit:
    summary = MITOZ_ANNOTATE.out.summary.mix(MITOS2_ANNOTATE.out.summary)
}
