# Phenotype QC, BLUEs and BLUPs (ASReml-R)
# Run from the project root (folder containing data/ and results/)
library(asreml)
library(ggplot2)

out <- "results"
dir.create(out, showWarnings = FALSE)
save_plot <- function(p, file, w = 8, h = 6) {
  print(p)
  ggsave(file.path(out, file), p, width = w, height = h, dpi = 300)
}

# ---- 1. Load data ----
dat <- read.csv("data/pheno.csv", stringsAsFactors = FALSE)
miss <- setdiff(c("Genotype", "Environment", "phenotype"), names(dat))
if (length(miss)) stop("Missing columns: ", paste(miss, collapse = ", "))

dat$Genotype    <- factor(dat$Genotype)
dat$Environment <- factor(dat$Environment)
dat$phenotype   <- as.numeric(dat$phenotype)

# ---- 2. Basic QC ----
cat("Obs:", nrow(dat), "| Genotypes:", nlevels(dat$Genotype),
    "| Environments:", nlevels(dat$Environment), "\n")
print(colSums(is.na(dat)))

env_summary <- do.call(rbind, lapply(split(dat$phenotype, dat$Environment), function(x) data.frame(
  N = sum(!is.na(x)), Missing = sum(is.na(x)),
  Mean = mean(x, na.rm = TRUE), SD = sd(x, na.rm = TRUE),
  Minimum = min(x, na.rm = TRUE), Maximum = max(x, na.rm = TRUE))))
env_summary <- cbind(Environment = rownames(env_summary), env_summary, row.names = NULL)
print(env_summary)
write.csv(env_summary, file.path(out, "phenotype_summary_by_environment.csv"), row.names = FALSE)

# ---- 3. Flag outliers within environment (1.5 x IQR; flagged, not removed) ----
q <- function(p) ave(dat$phenotype, dat$Environment, FUN = function(x) quantile(x, p, na.rm = TRUE))
iqr <- q(0.75) - q(0.25)
flag <- dat$phenotype < q(0.25) - 1.5 * iqr | dat$phenotype > q(0.75) + 1.5 * iqr
flag[is.na(flag)] <- FALSE

outliers <- dat[flag, c("Genotype", "Environment", "phenotype")]
print(table(outliers$Environment))
print(outliers)
write.csv(outliers, file.path(out, "potential_phenotype_outliers.csv"), row.names = FALSE)

# ---- 4. QC plots ----
theme_set(theme_bw(base_size = 14))
tilt <- theme(axis.text.x = element_text(angle = 30, hjust = 1))

save_plot(ggplot(dat, aes(Environment, phenotype, fill = Environment)) +
            geom_boxplot(width = 0.65, outlier.colour = "red", outlier.size = 2) +
            labs(title = "Phenotype by environment", y = "Phenotype") +
            tilt + theme(legend.position = "none"),
          "phenotype_boxplot.png")

save_plot(ggplot(dat, aes(phenotype)) +
            geom_histogram(bins = 25, colour = "black", fill = "steelblue") +
            facet_wrap(~ Environment, scales = "free") +
            labs(title = "Phenotype distributions", x = "Phenotype", y = "Frequency"),
          "phenotype_histograms.png", w = 10)

save_plot(ggplot(dat, aes(sample = phenotype)) + stat_qq() + stat_qq_line(colour = "red") +
            facet_wrap(~ Environment, scales = "free") +
            labs(title = "Phenotype Q-Q plots", x = "Theoretical", y = "Observed"),
          "phenotype_qq_plots.png", w = 10)

dat$Outlier_status <- ifelse(flag, "Potential outlier", "Not flagged")
save_plot(ggplot(dat, aes(Environment, phenotype, colour = Outlier_status)) +
            geom_jitter(width = 0.18, height = 0, alpha = 0.75) +
            scale_colour_manual(values = c("Not flagged" = "grey50", "Potential outlier" = "red")) +
            labs(title = "Potential outliers", y = "Phenotype", colour = NULL) +
            tilt + theme(legend.position = "bottom"),
          "potential_outliers_plot.png")

# ---- 5. Models (drop missing only) ----
adat <- droplevels(dat[complete.cases(dat[, c("phenotype", "Genotype", "Environment")]), ])
cat("Obs used:", nrow(adat), "\n")
na_opt <- na.method(x = "include", y = "include")

# BLUE: Genotype fixed
blue_model <- asreml(fixed = phenotype ~ Environment + Genotype, data = adat, na.action = na_opt)
print(summary(blue_model))
print(wald(blue_model))
BLUEs <- predict(blue_model, classify = "Genotype")$pvals
BLUEs <- BLUEs[order(BLUEs$predicted.value, decreasing = TRUE), ]
print(head(BLUEs, 10))
write.csv(BLUEs, file.path(out, "genotype_BLUEs.csv"), row.names = FALSE)

# BLUP: Genotype random
blup_model <- asreml(fixed = phenotype ~ Environment, random = ~ Genotype,
                     data = adat, na.action = na_opt)
print(summary(blup_model)$varcomp)
BLUPs <- predict(blup_model, classify = "Genotype")$pvals
BLUPs <- BLUPs[order(BLUPs$predicted.value, decreasing = TRUE), ]
print(head(BLUPs, 10))
write.csv(BLUPs, file.path(out, "genotype_BLUPs.csv"), row.names = FALSE)

# ---- 6. Combine and compare ----
cols <- c("Genotype", "predicted.value", "std.error")
res <- merge(setNames(BLUEs[, cols], c("Genotype", "BLUE", "BLUE_SE")),
             setNames(BLUPs[, cols], c("Genotype", "BLUP", "BLUP_SE")),
             by = "Genotype", all = TRUE)
res <- res[order(res$BLUP, decreasing = TRUE), ]
print(head(res, 10))
write.csv(res, file.path(out, "genotype_BLUEs_and_BLUPs.csv"), row.names = FALSE)

save_plot(ggplot(res, aes(BLUE, BLUP)) + geom_point(colour = "navy", alpha = 0.7) +
            geom_abline(linetype = "dashed", colour = "red") +
            labs(title = "BLUEs vs BLUPs", subtitle = "BLUPs shrink toward the mean"),
          "BLUE_BLUP_comparison.png", w = 7)

# ---- 7. Residual diagnostics ----
diag_df <- data.frame(Fitted = fitted(blup_model), Residual = residuals(blup_model))
save_plot(ggplot(diag_df, aes(Fitted, Residual)) + geom_point(alpha = 0.6, colour = "steelblue") +
            geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
            labs(title = "Residuals vs fitted"),
          "BLUP_residuals_vs_fitted.png", w = 7)
save_plot(ggplot(diag_df, aes(sample = Residual)) + stat_qq() + stat_qq_line(colour = "red") +
            labs(title = "Residual Q-Q plot", x = "Theoretical", y = "Residual"),
          "BLUP_residual_qq_plot.png", w = 7)

cat("Done. Results in:", out, "\n")
