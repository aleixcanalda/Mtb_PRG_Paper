
#libraries
library(dplyr)
library(ggplot2)
library(vroom)
library(data.table)
library(tidyverse)
library(ggpubr)
library(ggExtra)
library(plyr)

#long read assemblies

#get list of files from the quast output, downloaded from spartan quast output
files = list.files(path='long/',pattern = ".*tsv")
file.paths <<- as.vector(paste("long", files, sep = "/"))
#read each file
tables <- lapply(X = file.paths, FUN = read.table, header = TRUE,sep = "\t",row.names=1)
#bind all dataframes into one
quast = do.call(cbind.data.frame,tables)
#transpose
t.quast = as.data.frame(t(quast))
#make numeric
t.quast = t.quast %>% mutate_all(as.numeric)

#add NGS information, ont and hifi data from spartan
t.quast$NGS = "PB"
ont = scan("~/Downloads/oxford_all.txt", character(), quote = "")
pb_hifi = scan("~/Downloads/pacbio_hifi_list.txt", character(), quote = "")
t.quast = t.quast %>% 
  rownames_to_column('ID') %>%
  mutate(NGS = case_when(ID %in% ont ~ "ONT",
                         ID %in% pb_hifi~ "PB-Hifi",
                         TRUE ~ "PB"))
#ncbi assemblies
filesn = list.files(path='ncbi/',pattern = ".*tsv")
file.pathsn <<- as.vector(paste("ncbi", filesn, sep = "/"))
tablesn <- lapply(X = file.pathsn, FUN = read.table, header = TRUE,sep = "\t",row.names=1)
quastn = do.call(cbind.data.frame,tablesn)
t.quastn = as.data.frame(t(quastn))
t.quastn = t.quastn %>% mutate_all(as.numeric)
t.quastn$NGS = "NCBI"
t.quastn = t.quastn %>% 
  rownames_to_column('ID')
#put long read and ncbi data together
t.quast.all = rbind.fill(t.quast,t.quastn)

## FILTERING
t.quast.gen_fac.filt = t.quast.all[t.quast.all$Genome_fraction>95,]
t.quast.gen_fac.n50.filt = t.quast.gen_fac.filt[t.quast.gen_fac.filt$N50 > 10606,]
t.quast.gen_fac.n50.nopb.nomiss.filt = t.quast.gen_fac.n50.filt[t.quast.gen_fac.n50.filt$misassemblies<100 & t.quast.gen_fac.n50.filt$NGS!="PB",]

## PLOTS

# Plot Genome Fraction, FIRST filter
p_gf = ggplot(t.quast.all, aes(Genome_fraction, fill=NGS,alpha=0.05)) +
  geom_density(color="black") + 
  theme_pubr() +
  labs(x="Genome fraction (%)",y="Density") +
  scale_fill_manual(values=c("gray","yellow3","turquoise","lightblue"),name="") +
  geom_segment(aes(x=95,y=0,xend=95,yend=8),color="black",linetype=2) +
  scale_alpha(guide = "none")

# Plot of assembly length + N50, AFTER genome fraction filter
p_n50len = ggplot(t.quast.gen_fac.filt, aes(x=N50,y=Total_length, color=NGS,alpha=0.05)) +
  geom_point() + 
  theme_pubr() +
  labs(y="Length of the assembly") +
  scale_color_manual(values=c("gray","yellow3","turquoise","lightblue"),name="") +
  scale_alpha(guide = "none") +
  geom_segment(aes(x=10606,y=4e6,xend=10606,yend=1.4e7),colour="black",linetype=2) +
  ylim(4e6,NA)
  #geom_text(aes(x=4411532,y=1e5,label="H37Rv"),colour="black")
p_n50lenmarg = ggMarginal(p_n50len, type = "density",groupFill = TRUE)


# Plot of misassemblies, AFTER n50 filter
p_miss = ggplot(t.quast.gen_fac.n50.filt, aes(misassemblies, fill=NGS,alpha=0.05)) +
  geom_density(color="black") + 
  theme_pubr() +
  labs(x="Misassemblies",y="Density") +
  scale_fill_manual(values=c("gray","yellow3","turquoise","lightblue"),name="") +
  scale_alpha(guide = "none") +
  geom_segment(aes(x=100,y=0,xend=100,yend=0.55),color="black",linetype=2)

p_filt = ggplot(t.quast.gen_fac.n50.nopb.nomiss.filt, aes(x=N50,y=Total_length, color=NGS,alpha=0.05)) +
  geom_point() + 
  theme_pubr() +
  labs(y="Length of the assembly") +
  scale_color_manual(values=c("gray","yellow3","lightblue"),name="") +
  scale_alpha(guide = "none") +
  ylim(4.2e6,NA)
#geom_text(aes(x=4411532,y=1e5,label="H37Rv"),colour="black")
p_filtmarg = ggMarginal(p_filt, type = "density",groupFill = TRUE)

#final filtered data plot
png("~/Documents/PHD/SV/SV_analysis/figures/quast_figures/quast_filtered_95_24.png",width=750,height=1000)
ggarrange(p_gf,p_n50lenmarg,p_miss,p_filtmarg, ncol = 2, labels = c("a","b","c", "d"),nrow = 2)
dev.off()
write_csv(t.quast.gen_fac.n50.nopb.nomiss.filt,"~/Documents/PHD/SV/SV_analysis/quast/assemblies_filtered_final.csv")

