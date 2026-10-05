# Simulated chickpea VCF for a PLINK QC demo
rm(list = ls())
library(AlphaSimR)
library(data.table)
set.seed(20260717)

# ---- 1. Parameters ----
nInd <- 200; nChr <- 8
nQtl <- 25                 # per chr (200 total)
nSnp <- 7500               # per chr (60,000 total)
genLen <- 1                # Morgans per chr
PhyLen <- 1e8              # ~100 Mb per chr

currentNe <- 1000
histNe  <- c(2000, 5000, 10000)
histGen <- c(100, 1000, 10000)

MeanP <- 2500; VarGen <- 150000; VarGxE <- 150000; VarE <- 450000
envP     <- c(0.15, 0.50, 0.85)
envNames <- c("Env1_low", "Env2_mid", "Env3_high")
cat("Expected h2:", round(VarGen / (VarGen + VarGxE + VarE), 3), "\n")

# QC problems to inject
baseMiss <- 0.02
nBadSamples <- 15;   badSampMin <- 0.15; badSampMax <- 0.60
nBadMarkers <- 3000; badMarkMin <- 0.10; badMarkMax <- 0.45

outdir <- "sim_qc_demo"
dir.create(outdir, showWarnings = FALSE)
vcfFile <- file.path(outdir, "sim_chickpea.vcf")

# ---- 2. Founder haplotypes (double mutation rate on failure) ----
segSites <- nQtl + nSnp
founderPop <- NULL
mutRate <- 1e-7
for (attempt in 1:5) {
  cat("Attempt", attempt, "mutRate", format(mutRate, scientific = TRUE), "\n")
  founderPop <- tryCatch(
    runMacs2(nInd = nInd, nChr = nChr, segSites = segSites, Ne = currentNe,
             histNe = histNe, histGen = histGen, genLen = genLen, bp = PhyLen,
             mutRate = mutRate),
    error = function(e) { message("MaCS failed: ", conditionMessage(e)); NULL })
  if (!is.null(founderPop)) break
  mutRate <- mutRate * 2
}
if (is.null(founderPop)) stop("MaCS failed to produce ", segSites, " sites per chr after 5 attempts")

# ---- 3. Trait, SNP chip, panel ----
SP <- SimParam$new(founderPop)
SP$addTraitAG(nQtlPerChr = nQtl, mean = MeanP, var = VarGen, varGxE = VarGxE, name = "Simulated_trait")
SP$addSnpChip(nSnpPerChr = nSnp, name = "Chickpea_60K")
Panel <- newPop(founderPop, simParam = SP)

# ---- 4. Phenotypes in 3 environments ----
phenoLong <- rbindlist(lapply(seq_along(envP), function(i) {
  pop <- setPheno(Panel, varE = VarE, p = envP[i], simParam = SP)
  data.table(Genotype = pop@id, Environment = envNames[i], Env_p = envP[i],
             gv = gv(pop)[, 1], phenotype = pheno(pop)[, 1])
}))

realisedH2 <- phenoLong[, .(genetic_variance = var(gv), phenotypic_variance = var(phenotype),
                            realised_h2 = var(gv) / var(phenotype)), by = Environment]
print(realisedH2)

rG <- cor(as.matrix(dcast(phenoLong, Genotype ~ Environment, value.var = "gv")[, -1]))
print(round(rG, 3))

fwrite(phenoLong,  file.path(outdir, "phenotypes_3env.csv"))
fwrite(realisedH2, file.path(outdir, "realised_heritability.csv"))
fwrite(setNames(as.data.table(as.table(rG)), c("Environment1", "Environment2", "Genetic_correlation")),
       file.path(outdir, "genetic_correlations.csv"))

# ---- 5. SNP genotypes and map ----
geno <- pullSnpGeno(Panel, snpChip = 1, simParam = SP)
map  <- getSnpMap(snpChip = 1, simParam = SP)
stopifnot(ncol(geno) == nrow(map))

mapDT <- data.table(originalColumn = seq_len(nrow(map)),
                    chr = as.integer(as.character(map$chr)),
                    site = as.integer(map$site),
                    geneticPosition = as.numeric(map$pos))
chrCounts <- mapDT[, .(Number_of_SNPs = .N), by = chr][order(chr)]
print(chrCounts)
stopifnot(!anyNA(mapDT$chr), identical(chrCounts$chr, seq_len(nChr)),
          all(chrCounts$Number_of_SNPs == nSnp))
fwrite(chrCounts, file.path(outdir, "snp_counts_by_chromosome.csv"))

mafTrue <- colMeans(geno) / 2
mafTrue <- pmin(mafTrue, 1 - mafTrue)
cat("True MAF: min", round(min(mafTrue), 4), "median", round(median(mafTrue), 4),
    "max", round(max(mafTrue), 4), "| <0.05:", round(mean(mafTrue < 0.05), 4), "\n")

# ---- 6. Inject missingness ----
genoQC <- geno
genoQC[sample.int(length(genoQC), round(baseMiss * length(genoQC)))] <- NA_integer_

badSamples <- sample(seq_len(nInd), nBadSamples)
badSampleMiss <- numeric(nBadSamples)
for (k in seq_along(badSamples)) {
  badSampleMiss[k] <- runif(1, badSampMin, badSampMax)
  genoQC[badSamples[k], sample(seq_len(ncol(genoQC)), round(badSampleMiss[k] * ncol(genoQC)))] <- NA_integer_
}

badMarkers <- sample(seq_len(ncol(genoQC)), nBadMarkers)
badMarkerMiss <- numeric(nBadMarkers)
for (k in seq_along(badMarkers)) {
  badMarkerMiss[k] <- runif(1, badMarkMin, badMarkMax)
  genoQC[sample(seq_len(nInd), round(badMarkerMiss[k] * nInd)), badMarkers[k]] <- NA_integer_
}

sampleCallRate <- 1 - rowMeans(is.na(genoQC))
markerCallRate <- 1 - colMeans(is.na(genoQC))
cat("Sample call rate:", round(range(sampleCallRate), 3), "| median", round(median(sampleCallRate), 3), "\n")
cat("Marker call rate:", round(range(markerCallRate), 3), "| median", round(median(markerCallRate), 3), "\n")

sampleTruth <- data.table(Genotype = Panel@id, call_rate = round(sampleCallRate, 5),
                          truly_bad = seq_len(nInd) %in% badSamples, target_missing_rate = NA_real_)
sampleTruth[badSamples, target_missing_rate := badSampleMiss]
fwrite(sampleTruth, file.path(outdir, "truth_sample_quality.csv"))

# ---- 7. VCF positions (Morgans -> bp, strictly increasing per chr) ----
mapDT[, bp := pmax(1L, as.integer(round(geneticPosition / genLen * PhyLen)))]
setorder(mapDT, chr, bp, site, originalColumn)
mapDT[, bp := cummax(bp - seq_len(.N)) + seq_len(.N), by = chr]
if (max(mapDT$bp) > PhyLen) warning("Some positions exceed nominal chromosome length")
mapDT[, vcfID := paste0("SNP_", chr, "_", bp)]
stopifnot(!anyDuplicated(mapDT$vcfID))
genoQC <- genoQC[, mapDT$originalColumn, drop = FALSE]

markerMiss <- rep(NA_real_, ncol(geno))
markerMiss[badMarkers] <- badMarkerMiss
oc <- mapDT$originalColumn
fwrite(data.table(CHROM = mapDT$chr, POS = mapDT$bp, ID = mapDT$vcfID,
                  call_rate = round(markerCallRate[oc], 5), truly_bad = oc %in% badMarkers,
                  MAF_true = round(mafTrue[oc], 5), target_missing_rate = markerMiss[oc]),
       file.path(outdir, "truth_marker_quality.csv"))

# ---- 8. Write VCF ----
stopifnot(all(genoQC %in% c(0:2, NA)))
gt <- matrix(c("0/0", "0/1", "1/1")[t(genoQC) + 1L], ncol = nInd, dimnames = list(NULL, Panel@id))
gt[is.na(gt)] <- "./."

vcfBody <- cbind(data.table(`#CHROM` = mapDT$chr, POS = mapDT$bp, ID = mapDT$vcfID,
                            REF = "A", ALT = "T", QUAL = ".", FILTER = "PASS",
                            INFO = ".", FORMAT = "GT"),
                 as.data.table(gt))

contigs <- mapDT[, .(len = max(bp) + 1000L), by = chr][order(chr)]
writeLines(c("##fileformat=VCFv4.2",
             paste0("##fileDate=", format(Sys.Date(), "%Y%m%d")),
             "##source=AlphaSimR_simulation_teaching_dataset",
             "##reference=simulated_chickpea_genome",
             paste0("##contig=<ID=", contigs$chr, ",length=", contigs$len, ">"),
             '##FILTER=<ID=PASS,Description="All filters passed">',
             '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'),
           vcfFile)
fwrite(vcfBody, vcfFile, sep = "\t", append = TRUE, col.names = TRUE, quote = FALSE)

print(vcfBody[, .N, by = "#CHROM"])
stopifnot(nrow(vcfBody) == nChr * nSnp, identical(sort(unique(vcfBody[["#CHROM"]])), seq_len(nChr)))

# ---- 9. Shell commands for the demo ----
writeLines(c(
  "cd sim_qc_demo",
  "",
  "# Inspect",
  "bcftools view -h sim_chickpea.vcf | less",
  "bcftools view -H sim_chickpea.vcf | head",
  "bcftools view -H sim_chickpea.vcf | wc -l",
  "bcftools query -l sim_chickpea.vcf | wc -l",
  "bcftools query -f '%CHROM\\n' sim_chickpea.vcf | sort -n | uniq -c   # expect 7500 per chr 1-8",
  "",
  "# Compress, index, stats",
  "bgzip -f sim_chickpea.vcf",
  "tabix -f -p vcf sim_chickpea.vcf.gz",
  "bcftools stats sim_chickpea.vcf.gz > vcf_before_qc.stats",
  "",
  "# PLINK QC",
  "plink --vcf sim_chickpea.vcf.gz --double-id --allow-extra-chr --make-bed --out raw",
  "plink --bfile raw --missing --freq --out raw",
  "sort -k6,6gr raw.imiss | head -20",
  "plink --bfile raw --mind 0.20 --make-bed --out qc1_samples",
  "plink --bfile qc1_samples --geno 0.10 --maf 0.05 --make-bed --out qc2_final",
  "",
  "# Comparisons: wrong filter order; zero missingness allowed",
  "plink --bfile raw --geno 0.10 --maf 0.05 --mind 0.20 --make-bed --out qc_wrong_order",
  "plink --bfile raw --geno 0 --maf 0.05 --make-bed --out qc_geno0",
  "",
  "# PCA before/after QC and counts",
  "plink --bfile raw --pca 10 --out raw_pca",
  "plink --bfile qc2_final --pca 10 --out qc2_pca",
  "wc -l raw.bim raw.fam qc1_samples.bim qc1_samples.fam qc2_final.bim qc2_final.fam"
), file.path(outdir, "plink_qc_commands.sh"))

cat("Done:", nrow(vcfBody), "SNPs on", nChr, "chromosomes. Files in", normalizePath(outdir), "\n")
