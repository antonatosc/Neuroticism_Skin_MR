#### Custom R script to perform MR
#### Summary statistics should be on a single file, with setwd() pointing to the directory list.
#### Cols should be 
# SNP for rsID
# A1 for effect allele
# A2 for alternate allele
# BETA for effect estimate
# STDERR for Standard Error
# Freq for Allele frequency
# P for p-value
# N for sample size
# Clumping is performed with a local EUR reference data to avoid IEU GWAS database server overload. The LD reference dataset was downloaded from: http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz

uvmr = function(exp, outc, path_to_bfile, path_to_plink_bin){
  library(TwoSampleMR)
  library(data.table)
  print(paste0("Reading ", exp, ' instruments.\n'))
  exposure = fread(exp)
  print(paste0("Formatting ", exp, ' instruments.\n'))
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
  print(paste0("Clumping ", exp, ' instruments.\n'))
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
  print(paste0("Reading ", outc, ' instruments.\n'))
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
  print(paste0("Harmonizing: ", exp, " --> ", outc, ".\n"))
  harmonized = harmonise_data(
    exposure_dat = exposure, 
    outcome_dat = outcome,
    action = 2
  )
  harmonized$var_exp = varexp(harmonized)
  harmonized$f_stat = fstat(harmonized)
  harmonized$id.exposure = exp
  harmonized$id.outcome = outc
  print(paste0("Performing MR: ", exp, " --> ", outc, ".\n"))
  res = mr(harmonized)
  res$exposure = exp
  print(paste0("Performing MR-Egger: ", exp, " --> ", outc, ".\n"))
  pleio = mr_pleiotropy_test(harmonized)
  print(paste0("Performing Heterogeneity test: ", exp, " --> ", outc, ".\n"))
  het = mr_heterogeneity(harmonized)
  print(paste0("Performing single SNP analysis: ", exp, " --> ", outc, ".\n"))
  single = mr_singlesnp(harmonized)
  print(paste0("Performing LOO test: ", exp, " --> ", outc, ".\n"))
  loo = mr_leaveoneout(harmonized)
  print(paste0("Performing MR-PRESSO: ", exp, " --> ", outc, ".\n"))
  mrpresso = tryprintch({
    run_mr_presso(harmonized, NbDistribution = 1000, SignifThreshold = 0.05) },
    error=function(e){print("ERROR: Not enough IVs to run MRPRESSO.\n")})
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
  print(paste0("Final results for: ", exp, " --> ", outc, ".\n"))
  print(res)
  print(paste0("Mean F-statistic: ", exp, " --> ", outc, ".\n"))
  print(all[["F-statistic"]])
  print(paste0("Variance explained: ", exp, " --> ", outc, ".\n"))
  print(all[["Variance explained"]])
  print(paste0("Heterogeneity: ", exp, " --> ", outc, ".\n"))
  print(all[["I^2 heterogeneity"]])
  print(paste0("MR-Egger: ", exp, " --> ", outc, ".\n"))
  print(pleio)
  print(paste0("MR-PRESSO: ", exp, " --> ", outc, ".\n"))
  print(all[["presso"]][[1]]$`Main MR results`)
  return(all)
}
