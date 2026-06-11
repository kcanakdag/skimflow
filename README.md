# skimflow

Nextflow pipeline that turns Illumina genome skims into assemblies, mitogenomes, and BUSCO marker scores. Built for animal biodiversity work in the Bleidorn group, University of Göttingen.

## What it does

From paired-end (or single-end) Illumina skim reads, skimflow produces:

- **Trimmed reads** with fastp
- **Genome size + coverage estimate** with RESPECT (the "is this data good enough?" gate)
- **De novo assembly** with MEGAHIT (short reads) or Flye (long reads)
- **Mitogenome** with GetOrganelle
- **BUSCO marker scores** against the chosen lineage
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
# NHR standard96s: ssh u<youruser>@glogin-p3.hpc.gwdg.de
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
demo,reads/demo_R1.fastq.gz,reads/demo_R2.fastq.gz,Bostrychus_sinensis,1000000000
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
| `--busco_lineage NAME` | BUSCO lineage, default `metazoa_odb10`. |
| `--busco_db DIR` | Existing BUSCO download directory. |
| `--organelle_type NAME` | GetOrganelle target, default `animal_mt`. |
| `--organelle_db DIR` | Existing GetOrganelle database directory. |
| `--mitoz_clade NAME` | MitoZ clade, default `Arthropoda`. |
| `--mitoz_genetic_code N` | MitoZ mitochondrial genetic code, default `5`. |
| `--mitos2_genetic_code N` | MITOS2 mitochondrial genetic code, default `5`. |
| `--mitogenome_topology auto|linear|circular` | Topology hint for MitoZ/MITOS2, default `auto`. |
| `--flye_mode FLAG` | Override Flye mode, e.g. `--nano-raw` or `--pacbio-hifi`. |
| `--filtlong_min_length N` | Filtlong minimum read length, default `1000`. |
| `--filtlong_keep_percent N` | Filtlong retained-read percentage, default `90`. |
| `--megahit_preset NAME` | MEGAHIT preset, default `meta-sensitive`. |

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

## Output layout

Everything lands under `results/`:

```
results/
├── qc/<sample>/             fastp JSON + trimmed reads
├── decontam/<sample>/       kraken2 cleaned reads + report (when --kraken2_db is set)
├── genome_size/<sample>/    RESPECT estimates
├── long_read_qc/<sample>/   Filtlong-filtered reads (long-read samples)
├── assembly/<sample>/       MEGAHIT (short) or Flye (long) contigs
├── mitogenome/<sample>/     GetOrganelle FASTA + log (short-read samples)
├── mitogenome_annotation/   MitoZ / MITOS2 annotation outputs
├── markers/<sample>/        BUSCO short_summary + full output
└── report/                  multiqc_report.html
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

## Tools

| Step | Tool |
| --- | --- |
| Read QC | fastp |
| Decontamination (optional) | kraken2 |
| Genome size + coverage | RESPECT (Sayyari et al. 2022); needs a Gurobi licence |
| Short-read assembly | MEGAHIT |
| Long-read QC + assembly | Filtlong + Flye |
| Mitogenome | GetOrganelle |
| Mitogenome annotation | MitoZ / MITOS2 |
| Markers | BUSCO |
| Report | MultiQC |

All steps run in containers (Podman / Docker / Apptainer). No conda required at runtime.

## License

MIT. See `LICENSE` if added.

## Citation

If skimflow is useful for your work, please cite the underlying tools and link this repository.
