# 03 Trim and align

This stage is orchestrated by the Snakemake `trim_align` rule. It runs the Celatlas cutadapt and
STAR wrappers with the reference and library paths resolved from the sample sheet.
