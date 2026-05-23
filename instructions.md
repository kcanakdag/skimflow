Hi Kerem,

I hope you are doing great. As we discussed the other day, you can start your project now. You already have an idea of the overall goal, but here is some more detailed information.

The objective is to build a Nextflow pipeline to assemble genome skims, and extract both the mitochondrial genome and phylogenetic markers. It would also be great to include a genome size estimation step (using RESPECT) and generate a final summary report.

Some previous work was created in the department that you can use for inspiration:

https://github.com/ThiloSchulze/eukaryotic-genome-assembly 

https://github.com/ThiloSchulze/mitogenome-extraction 

https://github.com/fethalen/Patchwork 

The general structure of the pipeline and some programs you can use are in the attached md file. I will come back to you with some .fastq files that you can use for testing the pipeline.  You can come and work in the department, but you can also work from home, it is totally fine . Don't hesitate to ask any questions, we can always meet and discuss things in more detail :) 

Have a nice rest of the week, Mateo




## Pipeline Flowchart

1. **Read QC & Preprocessing:** Quality filtering and adapter trimming.
2. **Genome Size Estimation:** K-mer based calculation of the nuclear genome size.
3. **De Novo Assembly:** Assembling reads into contigs.
4. **Mitogenome Extraction:** Pulling out mitochondrial sequences and annotating them.
5. **Marker Extraction:** Identifying conserved phylogenetic markers (e.g., BUSCOs, rRNA).
6. **Reporting:** Collating results into a final user-friendly report.

## Included Tools & Software
This pipeline relies on the following standard bioinformatics tools:
* **QC & Trimming:** `fastp` / `FastQC`
* **Genome Size Estimation:** `RESPECT`, `KMC`
* **Assembly:** `SPAdes`
* **Mitogenome Extraction:** `GetOrganelle` / `MitoZ`
* **Marker Extraction:** `BUSCO`, `BLAST+`
* **Reporting:** `MultiQC`
