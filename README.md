# skimflow

Nextflow pipeline that turns Illumina genome skims into assemblies, mitogenomes, and BUSCO marker scores. Built for animal biodiversity work in the Bleidorn group, University of Göttingen.

## What it does

From paired-end (or single-end) Illumina skim reads, skimflow produces:

- **Trimmed reads** with fastp
- **Genome size + coverage estimate** with RESPECT (the "is this data good enough?" gate)
- **De novo assembly** with SPAdes
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

Working example: `assets/samplesheet_test.csv`.

## Output layout

Everything lands under `results/`:

```
results/
├── qc/<sample>/             fastp JSON + trimmed reads
├── genome_size/<sample>/    RESPECT estimates
├── assembly/<sample>/       SPAdes contigs, scaffolds, log
├── mitogenome/<sample>/     GetOrganelle output
├── markers/<sample>/        BUSCO short_summary + full output
└── report/                  multiqc_report.html
```

## Tools

| Step | Tool |
| --- | --- |
| Read QC | fastp |
| Genome size + coverage | RESPECT (Sayyari et al. 2022); needs a Gurobi licence |
| Assembly | SPAdes |
| Mitogenome | GetOrganelle |
| Markers | BUSCO |
| Report | MultiQC |

All steps run in containers (Podman / Docker / Apptainer). No conda required at runtime.

## License

MIT. See `LICENSE` if added.

## Citation

If skimflow is useful for your work, please cite the underlying tools and link this repository.
