# genomic-phenotypic-data-qc
Quality control and filtering of simulated genotype and phenotype data. The R script simulates the SNP data, PLINK and bcftools filter it, and ASReml-R checks the phenotypes.

## Contents

- `scripts/simulation.R` – simulates chickpea SNP genotypes (VCF) and phenotypes in three environments (AlphaSimR)
- `scripts/qc.txt` – genotype QC and filtering with bcftools and PLINK (run from `data/`)
- `scripts/phenotypes.R` – phenotype QC, BLUEs and BLUPs with ASReml-R (run from the repository root)
- `data/snp.vcf` – simulated SNP data (200 samples, 60,000 SNPs)
- `data/pheno.csv` – simulated phenotypes
