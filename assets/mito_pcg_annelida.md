# Reference proteins for self-seeding

`mito_pcg_annelida.faa` holds the 13 mitochondrial protein-coding genes of
eight annelids (104 proteins). `FIND_SEED` (modules/mito_seed.nf) searches each
sample's MEGAHIT contigs with them (tblastn, mitochondrial genetic code) to
find the sample's own mitochondrial contigs, which then seed GetOrganelle.
The proteins only locate contigs; no reference sequence enters the seed.

Headers are `<gene>|<species>|<protein accession>`. Proteins are the CDS
translations of the NCBI RefSeq mitogenomes below (efetch
`rettype=fasta_cds_aa`, fetched 2026-11-03), one genome per genus, chosen to
span the annelid tree. The genera of the test set (Pilargis, Ancistrosyllis,
Synelmis, Antonbruunia, Trypanobia) are left out so the validation stays
independent.

| Genus | RefSeq | Species |
| --- | --- | --- |
| Owenia | NC_028712.1 | Owenia fusiformis |
| Chaetopterus | NC_028710.1 | Chaetopterus variopedatus |
| Platynereis | NC_052845.1 | Platynereis cf. australis |
| Phyllodoce | NC_087881.1 | Phyllodoce medipapillata |
| Eusyllis | NC_031402.1 | Eusyllis blomstrandi |
| Ramisyllis | NC_065765.1 | Ramisyllis kingghidorahi |
| Lumbricus | NC_001673.1 | Lumbricus terrestris |
| Sabella | NC_056279.1 | Sabella spallanzanii |

For other phyla, pass a matching protein set with `--mito_seed_ref` and the
right `--mito_seed_gencode` (mitochondrial proteins are conserved enough that
the annelid set also finds fish contigs under code 2, as in the smoke test).
