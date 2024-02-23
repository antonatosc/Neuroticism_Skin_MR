#### Custom R script to perform MR
#### Summary statistics should be on a single file, with setwd() pointing to the directory list.
#### Cols should be 
# SNP for rsID, A1 for effect allele, A2 for alternate allele, BETA for effect estimate, STDERR for Standard Error, Freq for Allele frequency, P for p-value and N for sample size
# The code further tests for the IV strength for MR-Egger. If IGX2 is < 0.9, then SIMEX is performed.
# The final output is a list containing all relevant results (accompanied by SIMEX results if IGX < 0.9).
# Clumping is performed with a local EUR reference data to avoid IEU GWAS database server overload. The LD reference dataset was downloaded from: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz

#Load packages
#install.packages("devtools")
library(devtools)
#install_github("MRCIEU/TwoSampleMR")
library(TwoSampleMR)
#install.packages("ggplot2")
library(ggplot2)
#install.packages("data.table")

uvmr = function(exp, outc, path_to_bfile, path_to_plink_bin){
  print(paste0("Reading ", exp, ' instruments.'))
  exposure = fread(exp)
  print(paste0("Formatting ", exp, ' instruments.'))
  exposure[["P"]] = as.numeric(exposure[["P"]]) 
  exposure = exposure[exposure[["P"]] < 5e-8, ]
  exposure = as.data.frame(exposure)
  exposure = format_data(exposure,
                         type = "exposure",
                         header = TRUE,
                         snp_col = "SNP",
                         beta_col = "BETA",
                         se_col = "STDERR",
                         eaf_col = "Freq",
                         effect_allele_col = "A1",
                         other_allele_col = "A2",
                         pval_col = "P",
                         samplesize_col = "N"
  )
  print(paste0("Clumping ", exp, ' instruments.'))
  library(ieugwasr)
  bfile = path_to_bfile
  plink_bin = path_to_plink_bin
  exposure_clumped = ld_clump(dplyr::tibble(rsid=exposure$SNP,
                                            pval=exposure$pval.exposure,
                                            id=exposure$id.exposure),
                              bfile = bfile,
                              plink_bin = plink_bin)
  exposure = exposure[exposure$SNP %in% exposure_clumped$rsid, ]
  exposure = as.data.frame(exposure)
  print(paste0("Reading ", outc, ' instruments.'))
  outcome =read_outcome_data(
    snps = exposure$SNP,
    filename = outc,
    sep = "\t",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "STDERR",
    eaf_col = "Freq",
    effect_allele_col = "A1",
    other_allele_col = "A2",
    pval_col = "P",
    samplesize_col = "N"
  )
  print(paste0("Harmonizing: ", exp, " --> ", outc, "."))
  harmonized = harmonise_data(
    exposure_dat = exposure, 
    outcome_dat = outcome,
    action = 2
  )
  harmonized$var_exp = varexp(harmonized)
  harmonized$f_stat = fstat(harmonized)
  harmonized$id.exposure = exp
  harmonized$id.outcome = outc
  print(paste0("Performing MR: ", exp, " --> ", outc, "."))
  res = mr(harmonized)
  res$exposure = exp
  print(paste0("Performing MR-Egger: ", exp, " --> ", outc, "."))
  pleio = mr_pleiotropy_test(harmonized)
  print(paste0("Performing Heterogeneity test: ", exp, " --> ", outc, "."))
  het = mr_heterogeneity(harmonized)
  print(paste0("Performing single SNP analysis: ", exp, " --> ", outc, "."))
  single = mr_singlesnp(harmonized)
  print(paste0("Performing LOO test: ", exp, " --> ", outc, "."))
  loo = mr_leaveoneout(harmonized)
  print(paste0("Performing MR-PRESSO: ", exp, " --> ", outc, "."))
  mrpresso = tryCatch({
    run_mr_presso(harmonized, NbDistribution = 1000, SignifThreshold = 0.05) },
    error=function(e){print("ERROR: Not enough IVs to run MRPRESSO.")})
  all = list(
    "data" = harmonized,
    "results" = res,
    "pleiotropy" = pleio,
    "heterogeneity" = het,
    "single-snp" = single,
    "leave-one-out" = loo,
    "presso" = mrpresso
  )
  all[["Variance explained"]] = sum(harmonized[harmonized$mr_keep == T, ]$var_exp)
  all[["F-statistic"]] = mean(harmonized[harmonized$mr_keep == T, ]$f_stat)
  all[["IGX2"]] <- Isq(abs(harmonized[harmonized$mr_keep == T, ]$beta.exposure), harmonized[harmonized$mr_keep == T, ]$se.exposure)
  all[["I^2 heterogeneity"]] <- paste0(round((100*(het$Q[2]-het$Q_df[2])/het$Q),2), "%")
  if(all[["IGX2"]] > 0.9) {
    print(paste0("I-squared statistic for IV strength for MR-Egger: ", round(all[["IGX2"]], 2), ". Analysis complete."))
    print(paste0("Final results for: ", exp, " --> ", outc, "."))
    print(res)
    print(paste0("Mean F-statistic: ", exp, " --> ", outc, "."))
    print(all[["F-statistic"]])
    print(paste0("Variance explained: ", exp, " --> ", outc, "."))
    print(all[["Variance explained"]])
    print(paste0("Heterogeneity: ", exp, " --> ", outc, "."))
    print(all[["I^2 heterogeneity"]])
    print(paste0("MR-Egger: ", exp, " --> ", outc, "."))
    print(pleio)
    print(paste0("I-squared statistic for IV strength for MR-Egger: ", exp, " --> ", outc, "."))
    print(all[["IGX2"]])
    print(paste0("MR-PRESSO: ", exp, " --> ", outc, "."))
    print(all[["presso"]][[1]]$`Main MR results`)
  } else {
    print(paste0("I-squared statistic for IV strength for MR-Egger: ", round(all[["IGX2"]], 2), ", less than 0.9. Performing SIMEX."))
    #install.packages("simex")
    library(simex)
    df <- harmonized[harmonized$mr_keep == TRUE, ]
    beta.exp = df$beta.exposure
    beta.outc = df$beta.outcome
    se.exp = df$se.exposure
    se.outc = df$se.outcome
    BXG = abs(beta.exp)
    BYG = beta.outc*sign(beta.exp)
    fit = lm(BYG~BXG, x=T, y=T)
    print("Simulation extrapolation.")
    simex_fit = simex(fit,
                      B=1000,
                      measurement.error = se.exp,
                      SIMEXvariable="BXG",
                      fitting.method ="quad",
                      asymptotic="FALSE")
    simex_fit2 = summary(simex_fit)
    SIMEX_egger = data.frame(
      id.exposure = exp,
      id.outcome = outc,
      outcome = 'outcome',
      exposure = 'exposure',
      beta = simex_fit2$coefficients$jackknife[2,1],
      se = simex_fit2$coefficients$jackknife[2,2],
      pval = simex_fit2$coefficients$jackknife[2,4]
    )
    all[["SIMEX"]] = SIMEX_egger
    print(paste0("Final results for: ", exp, " --> ", outc, "."))
    print(res)
    print(paste0("Mean F-statistic: ", exp, " --> ", outc, "."))
    print(all[["F-statistic"]])
    print(paste0("Variance explained: ", exp, " --> ", outc, "."))
    print(all[["Variance explained"]])
    print(paste0("Heterogeneity: ", exp, " --> ", outc, "."))
    print(all[["I^2 heterogeneity"]])
    print(paste0("MR-Egger: ", exp, " --> ", outc, "."))
    print(pleio)
    print(paste0("I-squared statistic for IV strength for MR-Egger: ", exp, " --> ", outc, "."))
    print(all[["IGX2"]])
    print(paste0("SIMEX: ", exp, " --> ", outc, "."))
    print(paste0("MR-PRESSO: ", exp, " --> ", outc, "."))
    print(all[["SIMEX"]])
    print(all[["presso"]][[1]]$`Main MR results`)
  }
  return(all)
}
varexp <- function(data){
  pve = (2*(data$beta.exposure^2)*data$eaf.exposure*(1 - data$eaf.exposure))/
    ((2*(data$beta.exposure^2)*data$eaf.exposure*(1 - data$eaf.exposure)) + ((data$se.exposure^2)*2*data$samplesize.exposure*data$eaf.exposure*(1 - data$eaf.exposure)))
  return(pve)
}
fstat <- function(data){
  F_stat = ((data$samplesize.exposure - 1 - 1)/1)*(data$var_exp/(1 - data$var_exp))
  return(F_stat)
}
