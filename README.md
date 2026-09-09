# skimflow

Nextflow pipeline that turns Illumina genome skims into assemblies, mitogenomes, and BUSCO marker scores. Built for animal biodiversity work in the Bleidorn group, University of Göttingen.

## What it does

From paired-end (or single-end) Illumina skim reads, skimflow produces:

- **Trimmed reads** with fastp
- **Genome size + coverage estimate** with RESPECT (the "is this data good enough?" gate)
- **De novo assembly** with MEGAHIT (short reads) or Flye (long reads)
- **Mitogenome** with GetOrganelle, annotated with MitoZ and MITOS2
- **Per-gene mitochondrial multi-FASTAs** harvested from each sample's MITOS2 annotation and pooled across all samples into one file per gene (e.g. `COX1.fasta` with every species' COX1), plus a gene-occupancy matrix for downstream phylogenetics
- **Phylogenetic trees** with MAFFT + trimAl + IQ-TREE: each mito gene is aligned and trimmed, then concatenated into a partitioned supermatrix for a maximum-likelihood species tree (both nucleotide and amino-acid, plus optional per-gene trees)
- **BUSCO marker scores** against the chosen lineage
- **Per-sample metrics table** (`summary/skimflow_metrics.csv`) joining read QC, genome size, assembly stats, mitogenome length, gene occupancy, and BUSCO into one row per sample
- **One-page MultiQC HTML report** that summarises everything above

## Quick start

### Local

```bash
# 1. one-time setup
podman build -t local/respect:0.2 -f containers/respect/Containerfile containers/respect
./scripts/fetch_test_data.sh

# 2. get a free academic Gurobi WLS licence from https://www.gurobi.com/academia/
#    and save it to ~/gurobi.lic

# 3. run
nextflow run . -profile test,podman
```

Other engines: `-profile test,docker` or `-profile test,apptainer`.

### GWDG HPC (Göttingen)

Run these commands in an SSH session on GWDG. The launcher uses the `gwdg`
profile (Slurm + Apptainer). FASTQ and CSV paths must point to files that are
visible on the cluster; absolute paths are safest for the one-liner.

```bash
ssh u<youruser>@glogin-p3.hpc.gwdg.de
# NHR standard96: ssh u<youruser>@glogin.hpc.gwdg.de
# NHR standard96s: ssh u<youruser>@glogin-p2.hpc.gwdg.de
```

Before the first run, save your Gurobi WLS licence (free academic licence from
https://www.gurobi.com/academia/) as `~/gurobi.lic` and lock down its
permissions:

```bash
chmod 600 ~/gurobi.lic
```

The `gwdg` profile runs every step through Apptainer and binds this file into
the RESPECT container at `/opt/gurobi/gurobi.lic` automatically; the licence is
never copied into the image. If your licence lives elsewhere, pass
`--gurobi-lic /path/to/gurobi.lic` or set `GUROBI_LIC` in `.env`. Runs without
RESPECT (`--skip-respect`, or any long-read-only run) need no licence.

The one-liner clones or updates the repo on GWDG, then starts the run.

**Single sample (no CSV needed):**

```bash
curl -fsSL https://raw.githubusercontent.com/kcanakdag/skimflow/main/bootstrap_gwdg.sh \
    | bash -s -- --r1 /path/to/sample_R1.fastq.gz --r2 /path/to/sample_R2.fastq.gz
```

`sample_id` is derived from the R1 filename (e.g. `Pmisa_R1.fastq.gz` becomes `Pmisa`); override with `--sample-id`. Drop `--r2` for single-end.

**Multi-sample with a CSV:**

```bash
curl -fsSL https://raw.githubusercontent.com/kcanakdag/skimflow/main/bootstrap_gwdg.sh \
    | bash -s -- --input my_samplesheet.csv
```

Optional kraken2 decontamination. First put a kraken2 database on the cluster
(see "Preparing a kraken2 database on GWDG" under Optional steps below), then
point the run at it:

```bash
curl -fsSL https://raw.githubusercontent.com/kcanakdag/skimflow/main/bootstrap_gwdg.sh \
    | bash -s -- --input my_samplesheet.csv --kraken2-db ~/.project/<projectid>/kraken2_db
```

If the repo is already cloned, run from inside the clone instead:

```bash
cd ~/projects/skimflow
./scripts/run_gwdg.sh \
    --r1 /path/to/sample_R1.fastq.gz \
    --r2 /path/to/sample_R2.fastq.gz \
    --sample-id sample_name
```

NHR users can add their partition to either command:

```bash
curl -fsSL https://raw.githubusercontent.com/kcanakdag/skimflow/main/bootstrap_gwdg.sh \
    | bash -s -- --partition standard96 --input my_samplesheet.csv
```

Optional: copy the example environment file for repeated settings such as
`GUROBI_LIC`, `KRAKEN2_DB`, `MITOS2_REFDIR`, or proxy settings:

```bash
cp .env.example .env
$EDITOR .env
set -a
source .env
set +a
```

Use `--gurobi-lic /path/to/gurobi.lic` if your licence is not at
`~/gurobi.lic`. Use `--skip-respect` to run without RESPECT/Gurobi. Use
`./scripts/run_gwdg.sh --help` for all options.

## Input options

For a single sample, point at the FASTQs directly with `--r1` / `--r2` (the GWDG quick-start above shows this). No samplesheet required.

For multiple samples, use a CSV with one row per sample:

```
sample_id,fastq_1,fastq_2,species_id,expected_genome_size_bp
demo,test_data/demo_R1.fastq.gz,test_data/demo_R2.fastq.gz,Bostrychus_sinensis,1000000000
```

- Empty `fastq_2` → single-end mode.
- `expected_genome_size_bp` is optional; if set, it appears alongside RESPECT's data-derived estimate in the MultiQC report as a sanity check.
- Relative paths are resolved against the repo root.
- For a long-read sample, give the row a `long_reads` column (and optional `lr_type`) instead of `fastq_1`/`fastq_2`; see Optional steps below.

Working example: `assets/samplesheet_test.csv`.

## Common flags

The GWDG launcher uses shell-style flag names (`--sample-id`). Raw Nextflow
commands use the parameter names from `nextflow.config` (`--sample_id`).
Everything after `--` in `./scripts/run_gwdg.sh ... -- ...` is passed directly
to Nextflow, so use the Nextflow spelling there.

| GWDG launcher | Nextflow | Meaning |
| --- | --- | --- |
| `--r1 PATH` | `--r1 PATH` | R1 FASTQ for one short-read sample. |
| `--r2 PATH` | `--r2 PATH` | R2 FASTQ; omit for single-end reads. |
| `--input CSV` | `--input CSV` | Samplesheet for one or more samples. |
| `--sample-id NAME` | `--sample_id NAME` | Sample name for direct `--r1` or `--long-reads` runs. |
| `--species NAME` | `--species NAME` | Species label for direct runs. |
| `--expected-size BP` | `--expected_size BP` | Expected genome size shown beside RESPECT's estimate. |
| `--long-reads PATH` | `--long_reads PATH` | One long-read FASTQ; skips RESPECT and GetOrganelle. |
| `--lr-type nanopore|pacbio` | `--lr_type nanopore|pacbio` | Long-read technology for Flye mode selection. |
| `--kraken2-db DIR` | `--kraken2_db DIR` | Enable Kraken2 decontamination for short reads. |
| `--mitos2-refdir DIR` | `--mitos2_refdir DIR` | Use an existing MITOS2 reference-data directory. |
| `--skip-respect` | `--skip_respect true` | Skip RESPECT and Gurobi. |
| `--gurobi-lic PATH` | `GUROBI_LIC=/path/to/gurobi.lic` | Use a Gurobi licence outside `~/gurobi.lic`. |

GWDG-only launcher flags:

| Flag | Meaning |
| --- | --- |
| `--partition NAME` | Slurm partition, e.g. `scc-cpu`, `standard96`, or `standard96s`. |
| `--account ID` | Slurm account override; usually unnecessary. |
| `--filesystem NAME` | Workspace filesystem, default `ceph-ssd`. |
| `--workspace-name NAME` | Workspace name, default `genome-skim`. |
| `--workspace-days N` | Initial workspace duration, default `30`. |

Common Nextflow-only overrides:

| Flag | Meaning |
| --- | --- |
| `--outdir DIR` | Results directory, default `results/`. |
| `--save_trimmed_reads` | Publish fastp's trimmed FASTQs under `read_qc/`; off by default (large intermediates). |
| `--busco_lineage NAME` | BUSCO lineage, default `metazoa_odb10`. |
| `--busco_db DIR` | Existing BUSCO download directory. |
| `--organelle_type NAME` | GetOrganelle target, default `animal_mt`. |
| `--organelle_db DIR` | Existing GetOrganelle database directory. |
| `--mitoz_clade NAME` | MitoZ clade, default `Arthropoda`. |
| `--mitoz_genetic_code N` | MitoZ mitochondrial genetic code, default `5`. |
| `--mitos2_genetic_code N` | MITOS2 mitochondrial genetic code, default `5`. |
| `--mitos2_refseqver NAME` | MITOS2 reference database version, default `refseq89m`. |
| `--mitos2_extra_args STR` | Extra flags passed to `runmitos`, default `--noplots`. |
| `--mitogenome_topology auto|linear|circular` | Topology hint for MitoZ/MITOS2, default `auto`. |
| `--flye_mode FLAG` | Override Flye mode, e.g. `--nano-raw` or `--pacbio-hifi`. |
| `--filtlong_min_length N` | Filtlong minimum read length, default `1000`. |
| `--filtlong_keep_percent N` | Filtlong retained-read percentage, default `90`. |
| `--megahit_preset NAME` | MEGAHIT preset, default `meta-sensitive`. |
| `--skip_phylo` | Skip the MAFFT/trimAl/IQ-TREE phylogeny stage. |
| `--phylo_gene_trees false` | Build only the supermatrix tree, not the per-gene trees; default `true`. |
| `--phylo_min_taxa_per_gene N` | Drop a gene aligned in fewer than N samples, default `4`. |
| `--mafft_args STR` | MAFFT flags, default `--auto`. |
| `--trimal_mode NAME` | trimAl method (`automated1`, `gappyout`, `strict`, ...), default `automated1`. |
| `--iqtree_bootstrap N` | IQ-TREE UFBoot replicates, default `1000`. |
| `--iqtree_alrt N` | IQ-TREE SH-aLRT replicates, default `1000`. |

## Optional steps

### Long-read assembly (Flye)

A sample is either short-read or long-read, never both. Long-read samples run
Filtlong then Flye and are scored by BUSCO; they do NOT run RESPECT genome-size
estimation or GetOrganelle mitogenome extraction (those are Illumina k-mer
tools).

Single long-read sample:

```bash
nextflow run . -profile podman --long_reads reads.fastq.gz --lr_type nanopore
```

`--lr_type` is `nanopore` (default) or `pacbio`. By default nanopore uses Flye's
`--nano-hq` (recommended for Guppy5+ / R10 chemistry) and pacbio uses
`--pacbio-raw`. Override with `--flye_mode '--nano-raw'` (legacy chemistry) or
`--flye_mode '--pacbio-hifi'` (HiFi reads).

In a samplesheet, give a row a `long_reads` column (and optional `lr_type`)
instead of `fastq_1`/`fastq_2`. Short and long rows can coexist in one sheet;
each is routed to the right assembler and both feed BUSCO.

### Decontamination (kraken2)

Short-read only, opt-in. Classify reads against a kraken2 DB and keep the
unclassified reads (the target) as the clean set feeding MEGAHIT, RESPECT, and
GetOrganelle:

```bash
nextflow run . -profile podman --input my.csv --kraken2_db /path/to/k2_db
```

When decontam is on, RESPECT runs on the decontaminated read set; keeping
unclassified reads can slightly shift the k-mer spectrum (small for animal
targets against a microbial-heavy DB such as PlusPFP). kraken2 loads the whole
DB into RAM (~17 to 20 GB for the 16 GB index), so the assigned node needs that
much free memory; `scc-cpu` and `standard96` both have plenty.

#### Preparing a kraken2 database on GWDG

A kraken2 database is a directory of large index files. Workspaces are
temporary and meant for active job data, so keep the database in permanent
project storage instead (it resolves via `~/.project/<projectid>/`). Download a
prebuilt index once from a login node, which has internet access:

```bash
mkdir -p ~/.project/<projectid>/kraken2_db
cd ~/.project/<projectid>/kraken2_db
# copy a current download link from https://benlangmead.github.io/aws-indexes/k2
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_YYYYMMDD.tar.gz
tar -xzf k2_standard_16gb_*.tar.gz && rm k2_standard_16gb_*.tar.gz
```

The capped Standard 16 GB index is a good default; PlusPF / PlusPFP add fungi,
protozoa, and plasmids at the cost of size and RAM. Once the directory holds the
`hash.k2d`, `opts.k2d`, and `taxo.k2d` files, point runs at it with
`--kraken2-db ~/.project/<projectid>/kraken2_db`, or set `KRAKEN2_DB` in `.env`
so every run picks it up without repeating the flag.

### Mitogenome annotation (MitoZ / MITOS2)

When GetOrganelle produces a mitogenome FASTA, skimflow annotates it with both
MitoZ and MITOS2. The defaults are for invertebrate mitochondrial genomes:
MitoZ clade `Arthropoda`, genetic code `5`, and MITOS2 `refseq89m`.

For vertebrates, override the clade and genetic code:

```bash
nextflow run . -profile podman --input my.csv \
    --mitoz_clade Chordata \
    --mitoz_genetic_code 2 \
    --mitos2_genetic_code 2
```

MITOS2 reference data is downloaded automatically on the first run. To use an
existing copy, point `--mitos2_refdir` at a directory containing entries such
as `refseq89m/`.

On GWDG, pass overrides as extra Nextflow arguments after `--`:

```bash
./scripts/run_gwdg.sh --input my.csv -- \
    --mitoz_clade Chordata \
    --mitoz_genetic_code 2 \
    --mitos2_genetic_code 2
```

### Phylogenetics (MAFFT / trimAl / IQ-TREE)

The harvested per-gene FASTAs feed a phylogenetics stage that runs by default
(disable with `--skip_phylo`). Each gene is aligned with MAFFT and trimmed with
trimAl on its own (whole-mitogenome alignment is avoided because gene order
varies across annelids and the control region is unalignable), then the trimmed
alignments are concatenated into a partitioned supermatrix and a maximum-
likelihood species tree is built with IQ-TREE (ModelFinder per partition,
`-B` ultrafast bootstrap, `-alrt` SH-aLRT). Both sequence types run: nucleotide
(13 PCG + 2 rRNA) and amino acid (13 PCG). Per-gene trees are also built unless
`--phylo_gene_trees false`.

The stage is fail-soft: any tree with fewer than four taxa is skipped rather
than aborting the run, so small or low-coverage batches still complete. Open the
resulting `.treefile` in FigTree or iTOL. Common overrides:

```bash
nextflow run . -profile podman --input my.csv \
    --trimal_mode gappyout \
    --iqtree_bootstrap 1000 \
    --phylo_min_taxa_per_gene 4
```

### Coverage-titration benchmark

`benchmarks/run_titration.sh` down-samples one high-coverage sample to a ladder
of target coverages (rasusa, several seeds each), runs the full pipeline at
every level, and collects each level's `summary/skimflow_metrics.csv` row into
one `titration.csv`. Use it to find the minimum skim depth for reliable
mitogenome / gene recovery. It is a standalone helper (not a Nextflow entry) and
runs rasusa inside the chosen container engine:

```bash
benchmarks/run_titration.sh \
    --r1 test_data/pmisa_small_R1.fastq.gz \
    --r2 test_data/pmisa_small_R2.fastq.gz \
    --genome-size 1200000000 \
    --coverages 0.05,0.1,0.25,0.5,1,2,5 --seeds 3 --engine podman \
    --organelle-db results/mitogenome/getorganelle_db \
    --busco-db results/markers/busco_downloads \
    -- -resume
```

Pass `--organelle-db` / `--busco-db` (and `-- -resume`) so the reference
databases are not re-downloaded for every level. `benchmarks/run_titration.sh --help`
lists all options.

## Output layout

Everything lands under `results/`:

```
results/
├── read_qc/                 fastp JSON + HTML (trimmed reads only with --save_trimmed_reads)
├── decontam/<sample>/       kraken2 cleaned reads + report (when --kraken2_db is set)
├── genome_size/<sample>/    RESPECT estimates
├── long_read_qc/<sample>/   Filtlong-filtered reads (long-read samples)
├── assembly/<sample>/       MEGAHIT (short) or Flye (long) contigs + assembly_stats.tsv
├── mitogenome/<sample>/     GetOrganelle FASTA + log (short-read samples)
├── mitogenome_annotation/   MitoZ / MITOS2 annotation outputs
├── genes/                   per-gene mito multi-FASTAs pooled across samples + occupancy matrix
├── phylogeny/{nt,aa}/       per-gene alignments/trims, gene trees, supermatrix + species tree
├── markers/<sample>/        BUSCO short_summary + full output
├── summary/                 skimflow_metrics.csv (one row per sample, all stages)
├── report/                  multiqc_report.html + multiqc_report_data/
├── pipeline_report.html     Nextflow execution report (resources per task)
├── pipeline_timeline.html   Nextflow timeline (shows the parallel fan-out)
└── pipeline_trace.txt       Nextflow trace (per-task CPU/RAM/wall-clock; feeds benchmarking)
```

GetOrganelle mitogenome output is:

- `results/mitogenome/<sample>/<sample>.mito.fasta` when GetOrganelle produced
  a mitochondrial `path_sequence` FASTA. This is the mitogenome sequence that
  MitoZ and MITOS2 annotate downstream.
- `results/mitogenome/<sample>/<sample>.mito.log.txt` for every attempted
  short-read sample. If no `<sample>.mito.fasta` is present, check this log; it
  usually means there were too few mitochondrial reads or GetOrganelle could not
  resolve a confident assembly.
- `results/mitogenome/getorganelle_db/` is the downloaded GetOrganelle
  seed/label database, shared by samples. It is not a sample result.

Long-read samples do not run GetOrganelle in this pipeline.

Annotation output is:

- `results/mitogenome_annotation/mitoz/<sample>/` for MitoZ outputs, including
  `<sample>.mitoz.log.txt` and `<sample>.mitoz.status.txt`.
- `results/mitogenome_annotation/mitos2/<sample>/` for MITOS2 outputs,
  including `<sample>.mitos2.gff`, `<sample>.mitos2.bed`,
  `<sample>.mitos2.log.txt`, and `<sample>.mitos2.status.txt` when produced.
- `results/report/multiqc_report.html` includes mitogenome assembly and
  annotation summary tables.

### Per-gene mitochondrial harvest

After MITOS2 annotates each sample, skimflow harvests the individual
mitochondrial genes (13 protein-coding genes plus the `16S` and `12S` rRNAs)
and pools them across every sample into one multi-FASTA per gene. This is the
input shape phylogenetics wants, and it stays useful even when a full
mitogenome assembly is too fragmentary to use directly.

- `results/genes/nt/<GENE>.fasta` is the nucleotide alignment input for one
  gene, with one record per sample that has it (e.g. `COX1.fasta` holds every
  species' COX1). Records are labelled by species; when two samples share a
  species label the sample id is appended to keep headers unique.
- `results/genes/aa/<GENE>.faa` is the matching amino-acid file for each
  protein-coding gene (no file for the rRNAs).
- `results/genes/occupancy.tsv` is the gene-by-sample matrix: the harvested
  nucleotide length per gene per sample (`0` where a gene is absent), with
  `species` and `flags` columns alongside. The MultiQC report shows the same
  matrix as its "Mito gene occupancy" section, there rendering absent genes as
  blank cells.
- `results/genes/per_sample/<sample>/` keeps each sample's own harvested gene
  FASTAs and occupancy row before pooling.

Harvest runs per sample with `errorStrategy 'ignore'`, so a sample MITOS2 could
not annotate simply contributes nothing rather than aborting the batch. With
near-zero-coverage test data no genes are harvested and the aggregator still
emits a header-only `occupancy.tsv`; real data is where this step produces
sequence.

### Phylogeny output

Under `results/phylogeny/`, one subtree per sequence type (`nt/` and `aa/`):

- `<type>/alignments/<GENE>.aln.fasta` and `<type>/trimmed/<GENE>.trimmed.fasta`
  are the MAFFT alignment and trimAl-trimmed alignment for each gene.
- `<type>/gene_trees/<GENE>.treefile` is the per-gene ML tree (with `.iqtree`
  report and `.contree`); a gene with fewer than four taxa gets a
  `<GENE>.skipped.txt` instead.
- `<type>/supermatrix.fasta` and `<type>/partitions.nex` are the concatenated
  alignment and its per-gene partition file; `<type>/concat_stats.tsv` records
  taxa/genes/length.
- `<type>/<type>_species.treefile` is the partitioned supermatrix species tree
  (with `.contree`, `.iqtree`, and run log). Open it in FigTree or iTOL.

### Per-sample metrics

`results/summary/skimflow_metrics.csv` has one row per sample joining read
counts and %Q30 (fastp), RESPECT genome size + coverage, assembly contig stats
(n contigs, total bp, N50), mitogenome length, MitoZ / MITOS2 feature counts,
number of target genes recovered, and BUSCO C/S/D/F/M. A stage that was skipped
or produced no result leaves its columns blank. A curated subset of the same
table appears in the MultiQC report as the "skimflow per-sample metrics"
section.

## Tools

| Step | Tool |
| --- | --- |
| Read QC | fastp |
| Decontamination (optional) | kraken2 |
| Genome size + coverage | RESPECT (Sarmashghi et al. 2021, PLOS Comput Biol); needs a Gurobi licence |
| Short-read assembly | MEGAHIT |
| Long-read QC + assembly | Filtlong + Flye |
| Mitogenome | GetOrganelle |
| Mitogenome annotation | MitoZ / MITOS2 |
| Per-gene mito harvest | MITOS2-annotation parser (bundled scripts) |
| Alignment + trimming | MAFFT + trimAl |
| Phylogenetics | IQ-TREE (supermatrix + per-gene, nt + aa) |
| Markers | BUSCO |
| Report | MultiQC |

All steps run in containers (Podman / Docker / Apptainer). No conda required at runtime.

## License

MIT. See `LICENSE` if added.

## Citation

If skimflow is useful for your work, please cite the underlying tools and link this repository.
