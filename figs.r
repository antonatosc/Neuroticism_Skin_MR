load("D:/data/skin_neuroticism/clusters/scripts/res.RData")
get_labels_and_metrics = function(df){
  df$OR = exp(df$Beta)
  df$loci = df$Beta - 1.96*df$SE
  df$upci = df$Beta + 1.96*df$SE
  df$LO_CI = exp(df$loci)
  df$UP_CI = exp(df$upci)
  df$label = paste0(round(df$OR, 3), " (", round(df$LO_CI, 3), "-", round(df$UP_CI, 3), ")")
  return(df)
}
###psoriasis
psoriasis_exposure_res = list()
for(i in names(diseases_as_exposure[["psoriasis"]])){
  psoriasis_exposure_res[["ivw"]] = rbind(psoriasis_exposure_res[["ivw"]], data.frame(
    Exposure = diseases_as_exposure[["psoriasis"]][[i]]$results$id.exposure[3],
    Outcome = diseases_as_exposure[["psoriasis"]][[i]]$results$id.outcome[3],
    Method = "IVW",
    NSNPs = diseases_as_exposure[["psoriasis"]][[i]]$results$nsnp[3],
    Beta = diseases_as_exposure[["psoriasis"]][[i]]$results$b[3],
    SE = diseases_as_exposure[["psoriasis"]][[i]]$results$se[3],
    Pval = diseases_as_exposure[["psoriasis"]][[i]]$results$pval[3]
  ))
  psoriasis_exposure_res[["mregger"]] = rbind(psoriasis_exposure_res[["mregger"]], data.frame(
    Exposure = diseases_as_exposure[["psoriasis"]][[i]]$results$id.exposure[1],
    Outcome = diseases_as_exposure[["psoriasis"]][[i]]$results$id.outcome[1],
    Method = "MR-Egger",
    NSNPs = diseases_as_exposure[["psoriasis"]][[i]]$results$nsnp[1],
    Beta = diseases_as_exposure[["psoriasis"]][[i]]$results$b[1],
    SE = diseases_as_exposure[["psoriasis"]][[i]]$results$se[1],
    Pval = diseases_as_exposure[["psoriasis"]][[i]]$results$pval[1]
  ))
  psoriasis_exposure_res[["median"]] = rbind(psoriasis_exposure_res[["median"]], data.frame(
    Exposure = diseases_as_exposure[["psoriasis"]][[i]]$results$id.exposure[2],
    Outcome = diseases_as_exposure[["psoriasis"]][[i]]$results$id.outcome[2],
    Method = "Weighted median",
    NSNPs = diseases_as_exposure[["psoriasis"]][[i]]$results$nsnp[2],
    Beta = diseases_as_exposure[["psoriasis"]][[i]]$results$b[2],
    SE = diseases_as_exposure[["psoriasis"]][[i]]$results$se[2],
    Pval = diseases_as_exposure[["psoriasis"]][[i]]$results$pval[2]
  ))
  psoriasis_exposure_res[["presso"]] = rbind(psoriasis_exposure_res[["presso"]] , data.frame(
    Exposure = diseases_as_exposure[["psoriasis"]][[i]]$results$id.exposure[2],
    Outcome = diseases_as_exposure[["psoriasis"]][[i]]$results$id.outcome[2],
    Method = "MR-PRESSO",
    NSNPs = diseases_as_exposure[["psoriasis"]][[i]]$results$nsnp[3] - length(diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
    Beta = ifelse(is.na(diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[2]),
                  diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[1],
                  diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[2]),
    SE = ifelse(is.na(diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$Sd[2]),
                diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$Sd[1],
                diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$Sd[2]),
    Pval = ifelse(is.na(diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[2]),
                  diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[1],
                  diseases_as_exposure[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[2])   
  ))
}
psoriasis_exposure_res = do.call(rbind, psoriasis_exposure_res)
psoriasis_exposure_res = get_labels_and_metrics(psoriasis_exposure_res)
psoriasis_exposure_res[psoriasis_exposure_res$Exposure == "psoriasis", ]$Exposure = "Psoriasis"
psoriasis_exposure_res[psoriasis_exposure_res$Outcome == "neuro", ]$Outcome = "Neuroticism"
psoriasis_exposure_res[psoriasis_exposure_res$Outcome == "worry_cluster", ]$Outcome = "Worry"
psoriasis_exposure_res[psoriasis_exposure_res$Outcome == "depressed_cluster", ]$Outcome = "Depressed affect"
psoriasis_exposure_res[psoriasis_exposure_res$Outcome == "SESA", ]$Outcome = "SESA"
psoriasis_exposure_res$Outcome = factor(psoriasis_exposure_res$Outcome, levels = c("Neuroticism", "Depressed affect", "SESA", "Worry"))
library(ggplot2)
library(dplyr)
psoriasis_exp_plot = psoriasis_exposure_res$Method <- factor(psoriasis_exposure_res$Method, 
                                                             levels = c("MR-PRESSO", "Weighted median", "MR-Egger", "IVW"))
psoriasis_exposure_plot = ggplot(psoriasis_exposure_res,
                                 aes(x = OR,
                                     y = reorder(Outcome, desc(Outcome)),
                                     color = Method)) +
  geom_point(size = 3,
             position = position_dodge(0.5)) +
  geom_errorbar(aes(xmin = LO_CI,
                    xmax = UP_CI),
                position = position_dodge(0.5),
                width = .2) +
  geom_vline(xintercept = 1,
             linetype = "dotted",
             color = "black") +
  scale_color_manual(values = c("IVW" = "#990438",
                                "MR-Egger" = "#FFC107",
                                "Weighted median" = "#004D40",
                                "MR-PRESSO" = "#1E88E5"),
                     breaks = c("IVW", "Weighted median", "MR-Egger", "MR-PRESSO")) +
  theme_classic(base_size = 16) +
  labs(x = "Odds Ratio (95% Confidence Interval)",
       y = "",
       title = "Psoriasis as exposure") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(colour = "white", fill = "white"),
    legend.position = "none"
  )
psoriasis_outcome_res = list()
for(i in names(diseases_as_outcome[["psoriasis"]])){
  psoriasis_outcome_res[["ivw"]] = rbind(psoriasis_outcome_res[["ivw"]], data.frame(
    Exposure = diseases_as_outcome[["psoriasis"]][[i]]$results$id.exposure[3],
    Outcome = diseases_as_outcome[["psoriasis"]][[i]]$results$id.outcome[3],
    Method = "IVW",
    NSNPs = diseases_as_outcome[["psoriasis"]][[i]]$results$nsnp[3],
    Beta = diseases_as_outcome[["psoriasis"]][[i]]$results$b[3],
    SE = diseases_as_outcome[["psoriasis"]][[i]]$results$se[3],
    Pval = diseases_as_outcome[["psoriasis"]][[i]]$results$pval[3]
  ))
  psoriasis_outcome_res[["simex"]] = rbind(psoriasis_outcome_res[["simex"]], data.frame(
    Exposure = diseases_as_outcome[["psoriasis"]][[i]]$results$id.exposure[1],
    Outcome = diseases_as_outcome[["psoriasis"]][[i]]$results$id.outcome[1],
    Method = "MR-Egger/SIMEX",
    NSNPs = diseases_as_outcome[["psoriasis"]][[i]]$results$nsnp[1],
    Beta = diseases_as_outcome[["psoriasis"]][[i]]$SIMEX$beta[1],
    SE = diseases_as_outcome[["psoriasis"]][[i]]$SIMEX$se[1],
    Pval = diseases_as_outcome[["psoriasis"]][[i]]$SIMEX$pval[1]
  ))
  psoriasis_outcome_res[["median"]] = rbind(psoriasis_outcome_res[["median"]], data.frame(
    Exposure = diseases_as_outcome[["psoriasis"]][[i]]$results$id.exposure[2],
    Outcome = diseases_as_outcome[["psoriasis"]][[i]]$results$id.outcome[2],
    Method = "Weighted median",
    NSNPs = diseases_as_outcome[["psoriasis"]][[i]]$results$nsnp[2],
    Beta = diseases_as_outcome[["psoriasis"]][[i]]$results$b[2],
    SE = diseases_as_outcome[["psoriasis"]][[i]]$results$se[2],
    Pval = diseases_as_outcome[["psoriasis"]][[i]]$results$pval[2]
  ))
  psoriasis_outcome_res[["presso"]] = rbind(psoriasis_outcome_res[["presso"]] , data.frame(
    Exposure = diseases_as_outcome[["psoriasis"]][[i]]$results$id.exposure[2],
    Outcome = diseases_as_outcome[["psoriasis"]][[i]]$results$id.outcome[2],
    Method = "MR-PRESSO",
    NSNPs = diseases_as_outcome[["psoriasis"]][[i]]$results$nsnp[3] - length(diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
    Beta = ifelse(is.na(diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[2]),
                  diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[1],
                  diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[2]),
    SE = ifelse(is.na(diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$Sd[2]),
                diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$Sd[1],
                diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$Sd[2]),
    Pval = ifelse(is.na(diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[2]),
                  diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[1],
                  diseases_as_outcome[["psoriasis"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[2])   
  ))
}
psoriasis_outcome_res = do.call(rbind, psoriasis_outcome_res)
psoriasis_outcome_res = get_labels_and_metrics(psoriasis_outcome_res)
psoriasis_outcome_res[psoriasis_outcome_res$Outcome == "psoriasis", ]$Outcome = "Psoriasis"
psoriasis_outcome_res[psoriasis_outcome_res$Exposure == "neuro", ]$Exposure = "Neuroticism"
psoriasis_outcome_res[psoriasis_outcome_res$Exposure == "worry_cluster", ]$Exposure = "Worry"
psoriasis_outcome_res[psoriasis_outcome_res$Exposure == "depressed_cluster", ]$Exposure = "Depressed affect"
psoriasis_outcome_res[psoriasis_outcome_res$Exposure == "SESA", ]$Exposure = "SESA"
psoriasis_outcome_res$Exposure = factor(psoriasis_outcome_res$Exposure, levels = c("Neuroticism", "Depressed affect", "SESA", "Worry"))
psoriasis_outc_plot = psoriasis_outcome_res$Method <- factor(psoriasis_outcome_res$Method, 
                                                             levels = c("MR-PRESSO", "Weighted median", "MR-Egger/SIMEX", "IVW"))
psoriasis_outcome_plot = ggplot(psoriasis_outcome_res,
                                aes(x = OR,
                                    y = reorder(Exposure, desc(Exposure)),
                                    color = Method)) +
  geom_point(size = 3,
             position = position_dodge(0.5)) +
  geom_errorbar(aes(xmin = LO_CI,
                    xmax = UP_CI),
                position = position_dodge(0.5),
                width = .2) +
  geom_vline(xintercept = 1,
             linetype = "dotted",
             color = "black") +
  scale_color_manual(values = c("IVW" = "#990438",
                                "MR-Egger/SIMEX" = "#FFC107",
                                "Weighted median" = "#004D40",
                                "MR-PRESSO" = "#1E88E5"),
                     breaks = c("IVW", "Weighted median", "MR-Egger/SIMEX", "MR-PRESSO")) +
  theme_classic(base_size = 16) +
  labs(x = "Odds Ratio (95% Confidence Interval)",
       y = "",
       title = "Psoriasis as outcome") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(colour = "white", fill = "white")
  ) +
  coord_cartesian(xlim = c(0.5, 3))
pdf("D:/data/skin_neuroticism/clusters/Fig1.pdf",
width = 12,
height = 5
)
cowplot::plot_grid(psoriasis_exposure_plot, psoriasis_outcome_plot, rel_widths = c(1, 1.2), scale = 1, labels = "auto")
dev.off()
###eczema
eczema_exposure_res = list()
for(i in names(diseases_as_exposure[["eczema"]])){
  eczema_exposure_res[["ivw"]] = rbind(eczema_exposure_res[["ivw"]], data.frame(
    Exposure = diseases_as_exposure[["eczema"]][[i]]$results$id.exposure[3],
    Outcome = diseases_as_exposure[["eczema"]][[i]]$results$id.outcome[3],
    Method = "IVW",
    NSNPs = diseases_as_exposure[["eczema"]][[i]]$results$nsnp[3],
    Beta = diseases_as_exposure[["eczema"]][[i]]$results$b[3],
    SE = diseases_as_exposure[["eczema"]][[i]]$results$se[3],
    Pval = diseases_as_exposure[["eczema"]][[i]]$results$pval[3]
  ))
  eczema_exposure_res[["mregger"]] = rbind(eczema_exposure_res[["mregger"]], data.frame(
    Exposure = diseases_as_exposure[["eczema"]][[i]]$results$id.exposure[1],
    Outcome = diseases_as_exposure[["eczema"]][[i]]$results$id.outcome[1],
    Method = "MR-Egger",
    NSNPs = diseases_as_exposure[["eczema"]][[i]]$results$nsnp[1],
    Beta = diseases_as_exposure[["eczema"]][[i]]$results$b[1],
    SE = diseases_as_exposure[["eczema"]][[i]]$results$se[1],
    Pval = diseases_as_exposure[["eczema"]][[i]]$results$pval[1]
  ))
  eczema_exposure_res[["median"]] = rbind(eczema_exposure_res[["median"]], data.frame(
    Exposure = diseases_as_exposure[["eczema"]][[i]]$results$id.exposure[2],
    Outcome = diseases_as_exposure[["eczema"]][[i]]$results$id.outcome[2],
    Method = "Weighted median",
    NSNPs = diseases_as_exposure[["eczema"]][[i]]$results$nsnp[2],
    Beta = diseases_as_exposure[["eczema"]][[i]]$results$b[2],
    SE = diseases_as_exposure[["eczema"]][[i]]$results$se[2],
    Pval = diseases_as_exposure[["eczema"]][[i]]$results$pval[2]
  ))
  eczema_exposure_res[["presso"]] = rbind(eczema_exposure_res[["presso"]] , data.frame(
    Exposure = diseases_as_exposure[["eczema"]][[i]]$results$id.exposure[2],
    Outcome = diseases_as_exposure[["eczema"]][[i]]$results$id.outcome[2],
    Method = "MR-PRESSO",
    NSNPs = diseases_as_exposure[["eczema"]][[i]]$results$nsnp[3] - length(diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
    Beta = ifelse(is.na(diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[2]),
                  diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[1],
                  diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[2]),
    SE = ifelse(is.na(diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`Main MR results`$Sd[2]),
                diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`Main MR results`$Sd[1],
                diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`Main MR results`$Sd[2]),
    Pval = ifelse(is.na(diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[2]),
                  diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[1],
                  diseases_as_exposure[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[2])   
  ))
}
eczema_exposure_res = do.call(rbind, eczema_exposure_res)
eczema_exposure_res = get_labels_and_metrics(eczema_exposure_res)
eczema_exposure_res[eczema_exposure_res$Exposure == "eczema", ]$Exposure = "eczema"
eczema_exposure_res[eczema_exposure_res$Outcome == "neuro", ]$Outcome = "Neuroticism"
eczema_exposure_res[eczema_exposure_res$Outcome == "worry_cluster", ]$Outcome = "Worry"
eczema_exposure_res[eczema_exposure_res$Outcome == "depressed_cluster", ]$Outcome = "Depressed affect"
eczema_exposure_res[eczema_exposure_res$Outcome == "SESA", ]$Outcome = "SESA"
eczema_exposure_res$Outcome = factor(eczema_exposure_res$Outcome, levels = c("Neuroticism", "Depressed affect", "SESA", "Worry"))
library(ggplot2)
library(dplyr)
eczema_exp_plot = eczema_exposure_res$Method <- factor(eczema_exposure_res$Method, 
                                                             levels = c("MR-PRESSO", "Weighted median", "MR-Egger", "IVW"))
eczema_exposure_plot = ggplot(eczema_exposure_res,
                                 aes(x = OR,
                                     y = reorder(Outcome, desc(Outcome)),
                                     color = Method)) +
  geom_point(size = 3,
             position = position_dodge(0.5)) +
  geom_errorbar(aes(xmin = LO_CI,
                    xmax = UP_CI),
                position = position_dodge(0.5),
                width = .2) +
  geom_vline(xintercept = 1,
             linetype = "dotted",
             color = "black") +
  scale_color_manual(values = c("IVW" = "#990438",
                                "MR-Egger" = "#FFC107",
                                "Weighted median" = "#004D40",
                                "MR-PRESSO" = "#1E88E5"),
                     breaks = c("IVW", "Weighted median", "MR-Egger", "MR-PRESSO")) +
  theme_classic(base_size = 16) +
  labs(x = "Odds Ratio (95% Confidence Interval)",
       y = "",
       title = "Atopic dermatitis as exposure") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(colour = "white", fill = "white"),
    legend.position = "none"
  )
eczema_outcome_res = list()
for(i in names(diseases_as_outcome[["eczema"]])){
  eczema_outcome_res[["ivw"]] = rbind(eczema_outcome_res[["ivw"]], data.frame(
    Exposure = diseases_as_outcome[["eczema"]][[i]]$results$id.exposure[3],
    Outcome = diseases_as_outcome[["eczema"]][[i]]$results$id.outcome[3],
    Method = "IVW",
    NSNPs = diseases_as_outcome[["eczema"]][[i]]$results$nsnp[3],
    Beta = diseases_as_outcome[["eczema"]][[i]]$results$b[3],
    SE = diseases_as_outcome[["eczema"]][[i]]$results$se[3],
    Pval = diseases_as_outcome[["eczema"]][[i]]$results$pval[3]
  ))
  eczema_outcome_res[["simex"]] = rbind(eczema_outcome_res[["simex"]], data.frame(
    Exposure = diseases_as_outcome[["eczema"]][[i]]$results$id.exposure[1],
    Outcome = diseases_as_outcome[["eczema"]][[i]]$results$id.outcome[1],
    Method = "MR-Egger/SIMEX",
    NSNPs = diseases_as_outcome[["eczema"]][[i]]$results$nsnp[1],
    Beta = diseases_as_outcome[["eczema"]][[i]]$SIMEX$beta[1],
    SE = diseases_as_outcome[["eczema"]][[i]]$SIMEX$se[1],
    Pval = diseases_as_outcome[["eczema"]][[i]]$SIMEX$pval[1]
  ))
  eczema_outcome_res[["median"]] = rbind(eczema_outcome_res[["median"]], data.frame(
    Exposure = diseases_as_outcome[["eczema"]][[i]]$results$id.exposure[2],
    Outcome = diseases_as_outcome[["eczema"]][[i]]$results$id.outcome[2],
    Method = "Weighted median",
    NSNPs = diseases_as_outcome[["eczema"]][[i]]$results$nsnp[2],
    Beta = diseases_as_outcome[["eczema"]][[i]]$results$b[2],
    SE = diseases_as_outcome[["eczema"]][[i]]$results$se[2],
    Pval = diseases_as_outcome[["eczema"]][[i]]$results$pval[2]
  ))
  eczema_outcome_res[["presso"]] = rbind(eczema_outcome_res[["presso"]] , data.frame(
    Exposure = diseases_as_outcome[["eczema"]][[i]]$results$id.exposure[2],
    Outcome = diseases_as_outcome[["eczema"]][[i]]$results$id.outcome[2],
    Method = "MR-PRESSO",
    NSNPs = diseases_as_outcome[["eczema"]][[i]]$results$nsnp[3] - length(diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
    Beta = ifelse(is.na(diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[2]),
                  diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[1],
                  diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`Causal Estimate`[2]),
    SE = ifelse(is.na(diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`Main MR results`$Sd[2]),
                diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`Main MR results`$Sd[1],
                diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`Main MR results`$Sd[2]),
    Pval = ifelse(is.na(diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[2]),
                  diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[1],
                  diseases_as_outcome[["eczema"]][[i]]$presso[[1]]$`Main MR results`$`P-value`[2])   
  ))
}
eczema_outcome_res = do.call(rbind, eczema_outcome_res)
eczema_outcome_res = get_labels_and_metrics(eczema_outcome_res)
eczema_outcome_res[eczema_outcome_res$Outcome == "eczema", ]$Outcome = "eczema"
eczema_outcome_res[eczema_outcome_res$Exposure == "neuro", ]$Exposure = "Neuroticism"
eczema_outcome_res[eczema_outcome_res$Exposure == "worry_cluster", ]$Exposure = "Worry"
eczema_outcome_res[eczema_outcome_res$Exposure == "depressed_cluster", ]$Exposure = "Depressed affect"
eczema_outcome_res[eczema_outcome_res$Exposure == "SESA", ]$Exposure = "SESA"
eczema_outcome_res$Exposure = factor(eczema_outcome_res$Exposure, levels = c("Neuroticism", "Depressed affect", "SESA", "Worry"))
eczema_outc_plot = eczema_outcome_res$Method <- factor(eczema_outcome_res$Method, 
                                                             levels = c("MR-PRESSO", "Weighted median", "MR-Egger/SIMEX", "IVW"))
eczema_outcome_plot = ggplot(eczema_outcome_res,
                                aes(x = OR,
                                    y = reorder(Exposure, desc(Exposure)),
                                    color = Method)) +
  geom_point(size = 3,
             position = position_dodge(0.5)) +
  geom_errorbar(aes(xmin = LO_CI,
                    xmax = UP_CI),
                position = position_dodge(0.5),
                width = .2) +
  geom_vline(xintercept = 1,
             linetype = "dotted",
             color = "black") +
  scale_color_manual(values = c("IVW" = "#990438",
                                "MR-Egger/SIMEX" = "#FFC107",
                                "Weighted median" = "#004D40",
                                "MR-PRESSO" = "#1E88E5"),
                     breaks = c("IVW", "Weighted median", "MR-Egger/SIMEX", "MR-PRESSO")) +
  theme_classic(base_size = 16) +
  labs(x = "Odds Ratio (95% Confidence Interval)",
       y = "",
       title = "Atopic dermatitis as outcome") +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks = element_line(color = "black"),
    strip.background = element_rect(colour = "white", fill = "white")
  ) +
  coord_cartesian(xlim = c(0.5, 3))
pdf("D:/data/skin_neuroticism/clusters/Fig2.pdf",
width = 12,
height = 5
)
cowplot::plot_grid(eczema_exposure_plot, eczema_outcome_plot, rel_widths = c(1, 1.2), scale = 0.9, labels = "auto")
dev.off()
