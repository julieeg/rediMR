# load required packages
library(dplyr)
library(data.table)


# command args
args <- commandArgs(trailingOnly = TRUE)
pheno <- args[1] #{pheno}
tag <- args[2]

pheno_tag=paste0(pheno, "_", tag)
ss_path <- paste0("../data/processed/gwas/", pheno_tag)


### Define functions
calc_lambda <- function(x, p=0.5){
  # Calculate genomic inflation lambda value
  x = x[!is.na(x)]
  x.quantile <- quantile(x, p)
  round(qchisq(1 - x.quantile, 1) / qchisq(p, 1), 2)
}

make_qq <- function(data, pval_col, main=""){
  # Make a quantile-quantile plot
  data <- filter(data, data[[pval_col]] > 0)  # In case extremely low p-values are stored as zero
  
  # Process p-values
  y_vals <- sort(-log10(data[[pval_col]]))
  x_vals <- -log10(rev(ppoints(length(y_vals))))  # ppoints generates a uniform probability distribution
  
  # Trim points at higher p-values (credit to RaMWAS package for code snippet)
  levels = as.integer((x_vals - x_vals[1]) / (tail(x_vals, 1) - x_vals[1]) * 2000)
  keep = c(TRUE, diff(levels) != 0)
  levels = as.integer((y_vals - y_vals[1])/(tail(y_vals, 1) - y_vals[1]) * 2000)
  keep = keep | c(TRUE, diff(levels) != 0)
  keep = which(keep)
  
  par(ps = 18)
  plot(x = x_vals[keep], y = y_vals[keep], 
       xlab = expression(-log[10](italic(p)) * " (Expected)"), 
       ylab = expression(-log[10](italic(p)) * " (Observed)"),
       main = main, cex = 0.8, 
       cex.lab = 0.8, cex.main = 0.9, 
       pch = 16, ylim = c(0, ceiling(max(y_vals))))
  abline(0, 1, lty = 2)
  legend(x = 'topleft', y = 'topleft',
         bquote(lambda == .(calc_lambda(data[[pval_col]]))), 
         cex = 0.9, bty = "n")
}


### Read in summary stats and subset to columns of interest
ss <- fread(paste0(ss_path, ".gwas")) %>% rename(CHR="#CHROM") %>%
  mutate(across(c(CHR, P), ~ as.numeric(.))) %>%
  select(CHR, POS, ID, REF, ALT, A1, TEST, OBS_CT, BETA, SE, T_STAT, P, ERRCODE)

# Save ssInput file with selected columns
readr::write_tsv(ss %>% rename("#CHROM"=CHR), paste0(ss_path, ".gwas"))

# File storage
plot_dir <- paste0(dirname(ss_path), "/gwas_plots")
system(paste0("mkdir -p ", plot_dir))


### Create Manhattan Plot

pdf(paste0(plot_dir, "/", pheno_tag, "_manhattan.pdf"), height = 5, width = 9)
qqman::manhattan(x=ss %>% filter(P<0.05), chr="CHR", bp="POS", p="P", snp="ID")
dev.off()


### Create Q-Q plot

write(calc_lambda(ss$P), paste0(plot_dir, "/", pheno, "_lambda"))
pdf(paste0(plot_dir, "/", pheno_tag, "_qq.pdf"))
make_qq(ss, "P")
dev.off()

##EOF

