setwd("D:/data/skin_neuroticism/clusters")
source("uvmr.r")
exposure = c("neuro", "depressed_cluster", "worry_cluster", "SESA")
disease = c("psoriasis", "eczema")
diseases_as_exposure <- list()
diseases_as_outcome <- list()
for(i in disease){
  for(j in exposure){
    diseases_as_exposure[[i]][[j]] = uvmr(i, j,
                                          "E:/MR/ld/EUR",
                                          "E:/MR/plink/plink.exe")
    diseases_as_outcome[[i]][[j]] = uvmr(j, i,
                                         "E:/MR/ld/EUR",
                                         "E:/MR/plink/plink.exe")
  }
}
variants_ecz_as_exp <- list()
variants_ecz_as_outc <- list()
for(i in exposure){
  variants_ecz_as_exp[[i]] = diseases_as_exposure[["eczema"]][[i]]$data[diseases_as_exposure[["eczema"]][[i]]$data$mr_keep == T, ]
  variants_ecz_as_outc[[i]] = diseases_as_outcome[["eczema"]][[i]]$data[diseases_as_outcome[["eczema"]][[i]]$data$mr_keep == T, ]
}
summary_data <- list()
for(i in c("eczema", exposure)){
  summary_data[[i]] = fread(i)
  summary_data[[i]] = summary_data[[i]][ , c("CHR", "POS", "SNP", "A1", "A2", "Freq", "BETA", "STDERR", "P", "N")]
  summary_data[[i]]$Z = summary_data[[i]]$BETA/summary_data[[i]]$STDERR
  names(summary_data[[i]]) = c("chr", "pos", "rsid", "ref", "alt", "freq", "beta", "se", "p", "N", "z")
}
library(MRlap)
for(i in exposure){
  diseases_as_exposure[["eczema"]][[i]][["mrlap"]] <- MRlap(exposure = summary_data[["eczema"]],
                                                                       exposure_name = "eczema",
                                                                       outcome = summary_data[[i]],
                                                                       outcome_name = i,
                                                                       ld = "D:/data/eur_w_ld_chr",
                                                                       hm3 = "D:/data/eur_w_ld_chr/w_hm3.snplist",
                                                                       do_pruning = F,
                                                                       user_SNPsToKeep = variants_ecz_as_exp[[i]]$SNP)
  diseases_as_outcome[["eczema"]][[i]][["mrlap"]] <- MRlap(exposure = summary_data[[i]],
                                                                      exposure_name = i,
                                                                      outcome = summary_data[["eczema"]],
                                                                      outcome_name = "eczema",
                                                                      ld = "D:/data/eur_w_ld_chr",
                                                                      hm3 = "D:/data/eur_w_ld_chr/w_hm3.snplist",
                                                                      do_pruning = F,
                                                                      user_SNPsToKeep = variants_ecz_as_outc[[i]]$SNP)
}
rm(summary_data)
save.image("res.RData")
