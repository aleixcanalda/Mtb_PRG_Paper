library(rworldmap)
library(dplyr)
library(ggplot2)
library(vroom)
library(data.table)
library(tidyverse)
library(ggpubr)
library(ggExtra)
library(plyr)
library(RColorBrewer)
library(tidyr)
library(gridExtra)
library(cowplot)

lin = read.table("~/Documents/PHD/SV/SV_analysis/temp/final_lineages.txt",sep=" ")

colnames(lin) = c("ID","Lineage")
lin$L1 = NA
lin$L2 = NA
lin$L3 = NA
lin$L4 = NA
lin$L5 = NA


for (i in 1:nrow(lin)) {
  if (length(str_split(lin$Lineage,"[.]")[[i]]) == 1) {
      lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
  }
  else if (length(str_split(lin$Lineage,"[.]")[[i]]) == 2) {
    lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
    lin[i,]$L2 = str_split(lin$Lineage,"[.]")[[i]][2]
  }
  
  else if (length(str_split(lin$Lineage,"[.]")[[i]]) == 3) {
    lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
    lin[i,]$L2 = str_split(lin$Lineage,"[.]")[[i]][2]
    lin[i,]$L3 = str_split(lin$Lineage,"[.]")[[i]][3]
  }
  
  else if (length(str_split(lin$Lineage,"[.]")[[i]]) == 4) {
    lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
    lin[i,]$L2 = str_split(lin$Lineage,"[.]")[[i]][2]
    lin[i,]$L3 = str_split(lin$Lineage,"[.]")[[i]][3]
    lin[i,]$L4 = str_split(lin$Lineage,"[.]")[[i]][4]
  }
  
  else if (length(str_split(lin$Lineage,"[.]")[[i]]) == 5) {
    lin[i,]$L1 = str_split(lin$Lineage,"[.]")[[i]][1]
    lin[i,]$L2 = str_split(lin$Lineage,"[.]")[[i]][2]
    lin[i,]$L3 = str_split(lin$Lineage,"[.]")[[i]][3]
    lin[i,]$L4 = str_split(lin$Lineage,"[.]")[[i]][4]
    lin[i,]$L5 = str_split(lin$Lineage,"[.]")[[i]][5]
  }
}    

lin_color = c("#DDCC77","#117733","#88CCEE","#CC6677","#882255","#AA4499","#332288","#44AA99","#DDAA55")

lin_count = as.data.frame(table(lin$Lineage))

lin_count$L1 = NA
for (i in 1:nrow(lin_count)) {
  lin_count[i,]$L1 = str_split(lin_count$Var1,"[.]")[[i]][1]
}

linfig = lin_count %>% arrange(L1,Freq) %>%
  mutate(Var1=factor(Var1,Var1)) %>%
  ggplot(aes(x=Var1,y=Freq,color=L1)) +
    geom_segment( aes(x=Var1 ,xend=Var1, y=0, yend=Freq),size=2, color="black") +
    geom_point(size=7) +
    theme_pubr() +
    theme(axis.text.y = element_text(size = 40),axis.text.x=element_text(angle=55,hjust = 1,size=40),legend.text = element_text(size = 40),
          axis.title = element_text(size = 40),legend.title = element_text(size = 40)) +
    scale_color_manual(values=lin_color) +
    labs(x="Lineages in PRG",color="Lineage") +
    ylab("Number of isolates") +
    scale_y_continuous(expand=c(0.01,0.01),limits=c(0,NA))

png("~/Documents/PHD/SV/SV_analysis/figures/r_figures/maplin.png",height = 2000,width = 2500)
plot_grid(map_grob,linfig,labels = c('a','b'),nrow = 2,ncol = 1,label_size = 40)
dev.off()