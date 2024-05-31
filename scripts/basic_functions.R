## Basic functions


# load packages
library(tidyverse) ; library(table1)



#####################################################
## ~~~~  Data display for descriptive tables  ~~~~ ##
#####################################################

# ======================================
## Print continuous vars as mean +- SD 
# ======================================

mean_sd<-function(x, d=2) {
  sprintf("%s \u00B1 %s", round(mean(x, na.rm=T), digits = d), 
          round(sd(x, na.rm=T), digits = d))
}


# ==================================
## Print categorical vars as n (%)
# ==================================

n_pct <- function(x, level=F) {
  if(level==F) {
    sapply(as.list(names(table(x))), function(lvl) {
      paste0(lvl, ", ", sum(x == lvl, na.rm=T), " (", round(sum(x == lvl, na.rm=T)/n()*100,1), "%)") }) } 
  else{paste0(sum(x == level, na.rm=T), " (", round(sum(x == level, na.rm=T)/n()*100,1), "%)")}
}


# =====================================================
## Print continuous vars as median [25th, 75th %tile]
# =====================================================

median_25to75<-function(x, d=2) {
  qs<-round(quantile(x, breaks=seq(0,1,0.25), na.rm=T), d)
  sprintf("%s [%s, %s]", round(median(x, na.rm=T), digits = d), 
          qs[[2]], qs[[4]])
}



################################################
##  ~~~~ Data wrangling & transformation ~~~~ ##
################################################

# ========================
## Remove outliers by SD
# ========================

remove_outliers.fun <- function(x, SDs=5) {
  bounds <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x <- ifelse(x>bounds[1] & x<bounds[2], x, NA)
  x
}


# ========================
## Median impute for negative or missing values
# ========================

median_imp.fun <- function(x) {
  x.new <- ifelse(x == -1 | x == -3 | x == -9 | is.na(x) == T, median(x, na.rm=T), x)
  return(x.new)
}


# ===================
## Calculate zscore
# ===================
zscore.fun <- function(x) {
  z<-((x - mean(x, na.rm=T)) / sd(x, na.rm=T))
  return(z)
}



#########################################################################
##  ~~~~  Basic descriptive plots of continuous/categorical vars  ~~~~ ##
#########################################################################

# ============================================
## pre-built color palettes & ggplot themes
# ============================================

library("paletteer") ; library("RColorBrewer")
palettes <- list(NatComms= paletteer_d("ggsci::nrc_npg", n=10),
                 greens5 = paletteer_dynamic("cartography::green.pal", 5),
                 greens = rev(paletteer_dynamic("cartography::green.pal", 10)),
                 pugr = c("#9C50CE", "#9C50CE95","#30900195", "#309001"),
                 #oranges5 = rev(paletteer_dynamic("cartography::orange.pal", 5)),
                 oranges = rev(paletteer_dynamic("cartography::orange.pal", 10)),
                 #blues5 = rev(paletteer_dynamic("cartography::blue.pal", 5)),
                 blues = rev(paletteer_dynamic("cartography::blue.pal", 10)),
                 purples=rev(paletteer_dynamic("cartography::purple.pal", 10))[3:10],
                 CovarGroups=c(brewer.pal(9,"Oranges")[4:6], brewer.pal(9, "Greens")[5:4], 
                               brewer.pal(9, "Blues")[5:6],  brewer.pal(9, "Purples")[c(8:4)], 
                               brewer.pal(9, "PiYG")[1:3], brewer.pal(9, "PuRd")[5:6]))


# =========================================
## Plot continuous x categorical
# =========================================

# Functions 
plot_contXcat <- function(cat_var, cont_var=pheno, data=dat, ...) {
  data %>% select(cat=all_of(cat_var), cont=all_of(cont_var)) %>% 
    filter(!is.na(cont)) %>%
    group_by(cat) %>% summarise(m=mean(cont, na.rm=T), sd=sd(cont, na.rm=T)) %>%
    ggplot(aes(x=cat, y=m, ymin=m-sd, ymax=m+sd, fill=cat)) + ggthemeF + 
    theme(axis.text.x=element_text(angle=35, hjust=1, color="black")) +
    geom_bar(stat="identity") + geom_errorbar(width=0.15) +
    scale_fill_manual(...) +
    labs(title=paste("Bar plot of", cont_var, "by", cat_var), x=" ", y="mean (SD)")
}

plot_contXcont <- function(cont_var_x, cont_var_y=pheno, data=dat, ...) {
  data %>% select(cont_x=all_of(cont_var_x), cont_y=all_of(cont_var_y)) %>% filter(complete.cases(.)) %>%
    ggplot(aes(x=cont_x, y=cont_y)) + ggthemeF +
    geom_point(size=1.25, position=position_jitter(0.15), alpha=0.75, ...) + 
    geom_smooth(method="lm", color = "black") + 
    labs(title=paste("Scatter plot of", cont_var_y, "by", cont_var_x), y=cont_var_y, x=cont_var_x)
}


plot_catxcat_pct <- function(cat_var_group, cat_var_x="phenoQ", data=dat, ...) {
  prop.table(table(data %>% select(cat_x=all_of(cat_var_x), cat_group=all_of(cat_var_group))),1) %>% 
    as.data.frame() %>%
    ggplot(aes(x=cat_x, y=Freq*100, group=cat_group, fill=cat_group)) +
    geom_bar(stat="identity", position=position_dodge()) +
    labs(title=paste(cat_var_x, "by", cat_var_group), y="Percent", x=" ") + 
    scale_fill_manual(...) +
    theme(plot.title = element_text(face="bold", size=10),
          axis.text.x = element_text(angle=25, vjust=0.95)) +
    ggthemeF 
}


#EOF







