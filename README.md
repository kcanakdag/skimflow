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

SSH in first. SCC users (default partition `scc-cpu`) go to Emmy Phase 3:

```bash
ssh u<youruser>@glogin-p3.hpc.gwdg.de
# NHR users (standard96 / standard96s): ssh u<youruser>@glogin.hpc.gwdg.de
```

One-liner that clones the repo, builds the container inside a compute job, allocates a workspace, and submits to SLURM. Two input styles:

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

Run `./scripts/run_gwdg.sh --help` from inside a clone to see all flags (`--partition`, `--filesystem`, `--gurobi-lic`, `--expected-size`, `--species`, ...).

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
DB into RAM (~17 to 20 GB for the 16 GB index).

## Output layout

Everything lands under `results/`:

```
results/
├── qc/<sample>/             fastp JSON + trimmed reads
├── decontam/<sample>/       kraken2 cleaned reads + report (when --kraken2_db is set)
├── genome_size/<sample>/    RESPECT estimates
├── long_read_qc/<sample>/   Filtlong-filtered reads (long-read samples)
├── assembly/<sample>/       MEGAHIT (short) or Flye (long) contigs
├── mitogenome/<sample>/     GetOrganelle output
├── markers/<sample>/        BUSCO short_summary + full output
└── report/                  multiqc_report.html
```

## Tools

| Step | Tool |
| --- | --- |
| Read QC | fastp |
| Decontamination (optional) | kraken2 |
| Genome size + coverage | RESPECT (Sayyari et al. 2022); needs a Gurobi licence |
| Short-read assembly | MEGAHIT |
| Long-read QC + assembly | Filtlong + Flye |
| Mitogenome | GetOrganelle |
| Markers | BUSCO |
| Report | MultiQC |

All steps run in containers (Podman / Docker / Apptainer). No conda required at runtime.

## License

MIT. See `LICENSE` if added.

## Citation

If skimflow is useful for your work, please cite the underlying tools and link this repository.
