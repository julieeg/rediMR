# postprocess rediMR

library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
covarset=args[1] #confounders1/basic

pairset="pairset1"
tag="vCole"
dataDir="../data/processed/rediMR/vCole"

diet_trait_pairs=c("alch_alt","alch_cir", "breadtype_cvd", "breadtype_ldl", "oilyfish_cvd","oilyfish_tg")
pairs.l <- list(
  oilyfish_tg = list(
    exposure="oilyfish_QT",
    outcome="tg"),
  oilyfish_cvd = list(
    exposure="oilyfish_QT",
    outcome="cvd"),
  bread_ldl = list(
    exposure="bread_type_BIN",
    outcome="ldl"),
  bread_cvd = list(
    exposure="bread_type_BIN",
    outcome="cvd"),
  alch_cir = list(
    exposure="alch_glasspermonth_QT",
    outcome="cir"),
  alch_alt = list(
    exposure="alch_glasspermonth_QT",
    outcome="alt")
)
  

# beta-change ---------------------
bchange_files <- grep(
  covarset, grep(
    paste0("_bchange_full_",tag,".csv"), 
    list.files(dataDir, recursive=T), value=T),
  value=T
) ; bchange_files


bchange_all <- do.call(rbind.data.frame, lapply(1:length(bchange_files), function(i) {
  exp=pairs.l[[i]]$exposure ; out=pairs.l[[i]]$outcome
  f <- grep(exp, grep(out, bchange_files, value=T), value=T)
  read.csv(paste0(dataDir, "/", f)) %>%
    mutate(pair=names(pairs.l[i]), exposure=exp, outcome=out) 
}))


## mr results.csv ---------------------
mr_results_files <- grep(
  covarset, grep(
    paste0("_mr_results_",tag,".csv"), 
    list.files(dataDir, recursive=T), value=T),
  value=T
) ; mr_results_files


mr_results_all <- do.call(rbind.data.frame, lapply(1:length(mr_results_files), function(i) {
  exp=pairs.l[[i]]$exposure ; out=pairs.l[[i]]$outcome
  res_file <- grep(exp, grep(out, mr_results_files, value=T), value=T)
  read.csv(paste0(dataDir, "/",res_file)) %>%
    mutate(pair=names(pairs.l[i]), exposure=exp, outcome=out)
  }))


# mr_output.rda ----------------------
mr_output_files <- grep(
  covarset, grep(paste0("_mr_output_",tag,".rda"), 
    list.files(dataDir, recursive=T), value=T), value=T
  ) ; mr_output_files

mr_output_all <- lapply(1:length(pairs.l), function(i) {
  exp=pairs.l[[i]]$exposure ; out=pairs.l[[i]]$outcome
  output_file <- grep(exp, grep(out, mr_output_files, value=T), value=T)
  readRDS(paste0(dataDir, "/",output_file)) 
}) ; names(mr_output_all) = names(pairs.l)


## mr_instruments.csv ---------
mr_instrum_files <- grep(
  covarset, grep(
    paste0("_mr_instrument_",tag,".csv"), 
    list.files(dataDir, recursive=T), value=T),
  value=T
) ; mr_instrum_files


mr_instrum_all <- do.call(rbind.data.frame, lapply(1:length(mr_instrum_files), function(i) {
  exp=pairs.l[[i]]$exposure ; out=pairs.l[[i]]$outcome
  instrum_file <- grep(exp, grep(out, mr_instrum_files, value=T), value=T) 
  read.csv(paste0(dataDir, "/",instrum_file)) %>%
    mutate(pair=names(pairs.l[i]), exposure=exp, outcome=out, .before=nSNPs)
}))


## save combined files based on covariate set (?)

bchange_all %>% write.csv(paste0(dataDir,"/", pairset, "_", covarset, "_",tag,"_bchange.csv"), row.names=F)
mr_results_all %>% write.csv(paste0(dataDir,"/", pairset, "_", covarset, "_",tag,"_mr_result.csv"), row.names=F)
mr_output_all %>% saveRDS(paste0(dataDir,"/", pairset, "_", covarset, "_", tag, "_mr_output.rda"))
mr_instrum_all %>% write.csv(paste0(dataDir,"/", pairset, "_", covarset, "_",tag,"_mr_instrument.csv"), row.names=F)


##EOF





