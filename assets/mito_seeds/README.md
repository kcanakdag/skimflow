# GetOrganelle seeds for the test set

COX1 sequences used as `mito_seed` for the three Pilargidae in
`assets/samplesheet_test.csv`. GetOrganelle's default `animal_mt` seeds
recruit no true mitochondrial reads from these taxa, so without a seed the
mitogenome step returns nothing (see `modules/mitogenome.nf`).

Each file is the COX1 CDS region cut from the species' published mitogenome
(NCBI nuccore, fetched with efetch `seq_start`/`seq_stop`). A single gene
stands in for the COX1 barcode a lab usually has for its specimens; the rest
of the mitogenome is assembled from the reads.

| File | Source | Region |
|---|---|---|
| `Ancistrosyllis_sp_COX1.fasta` | OR123445.1 | 4848..6398 |
| `Pilargis_wolfi_COX1.fasta` | OR123447.1 | 2666..4216 |
| `Synelmis_amoureuxi_COX1.fasta` | OR123442.1 | 2491..3213 and 5673..6497 (COX1 is split in this record, kept as two records) |

Antonbruunia and Trypanobia keep the default seeds; they assemble without one.
