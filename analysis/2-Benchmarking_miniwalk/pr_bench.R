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
library(cowplot)

#READ DATA
data_hist = read_csv("~/Documents/PHD/SV/SV_analysis/bench/pr_svimasm_NN.csv",col_names=c("TYPE","SIZE","VERDICT"))
pr_bench95 = read_tsv("~/Documents/PHD/SV/SV_analysis/bench/svimasm_longass_NN.tsv",col_names = c("Values","SV_Caller","Stat"))
short = read_tsv("~/Documents/PHD/SV/SV_analysis/bench/short_read_bench_NN.tsv",col_names = c("Values","SV_Caller","Stat"))
short_svimasm = read_tsv("~/Documents/PHD/SV/SV_analysis/bench/short_read_bench_svimasm_NN.tsv",col_names = c("Values","SV_Caller","Stat"))
data_hist_manta = read_csv("~/Documents/PHD/SV/SV_analysis/bench/pr_manta_histogram_NN.csv",col_names=c("TYPE","SIZE","VERDICT"))
data_hist_graph = read_csv("~/Documents/PHD/SV/SV_analysis/bench/pr_graph_histogram_NN.csv",col_names=c("TYPE","SIZE","VERDICT"))
data_hist_manta_svimasm = read_csv("~/Documents/PHD/SV/SV_analysis/bench/pr_manta_histogram_svimasm_NN.csv",col_names=c("TYPE","SIZE","VERDICT"))
data_hist_graph_svimasm = read_csv("~/Documents/PHD/SV/SV_analysis/bench/pr_graph_histogram_svimasm_NN.csv",col_names=c("TYPE","SIZE","VERDICT"))
cov = read_tsv("~/Documents/PHD/SV/SV_analysis/bench/short_read_bench_svimasm_changedprcode2_cov.tsv",col_names = c("Values","SV_Caller","Stat","Coverage"))
###

#FUNCTIONS
comparison_plot = function(short,title){
  short_manta = short[short$SV_Caller=="manta",]
  colnames(short_manta) = c("Values_Manta","SV_manta","Stat_manta")
  short_mg = cbind(short_manta,short[short$SV_Caller=="miniwalk-short-read",])
  comparison = ggplot(short_mg,aes(x=Values),fill=) + 
    geom_density(aes(x=Values*100,y=..density..,fill= "miniwalk-short-read")) +
    geom_density(aes(x=Values_Manta*100,y=-..density..,fill="manta")) +
    theme_pubr() +
    coord_flip() +
    xlim(0,100) +
    facet_grid(~Stat,switch = "x") +
    xlab("Values (%)") +
    scale_fill_manual(values = c("#404080","#69b3a2")) + 
    theme(strip.background = element_rect(fill = "white",color="white"),strip.text.x = element_text(size=40),axis.text.x = element_blank(),axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),axis.title.x = element_blank(),legend.title = element_blank(),plot.title = element_text(size = 40),
          legend.text = element_text(size=40),axis.text = element_text(size = 40),axis.title = element_text(size = 40))+
    scale_x_continuous(limits = c(0,100),breaks = seq(0,100,by=10)) +
    ggtitle(title)
  return(comparison)
}
#gets precision and recall distributed into different SV sizes and separated into INS and DEL
prepare_data <- function(data, size_col = "SIZE") {
  data2 = data %>%
    filter(TYPE != "INV") %>%
    mutate(Range = case_when(
      -3000 > SIZE ~ -3000,
      -3000 < SIZE & SIZE <= -2500 ~ -2500,
      -2500 < SIZE & SIZE <= -2000 ~ -2000,
      -2000 < SIZE & SIZE <= -1500 ~ -1500,
      -1500 < SIZE & SIZE <= -1000 ~ -1000,
      -1000 < SIZE & SIZE <= -500 ~ -500,
      -500 < SIZE & SIZE <= -100 ~ -100,
      -100 < SIZE & SIZE <= -50 ~ -50,
      100 > SIZE & SIZE >= 50 ~ 50,
      500 > SIZE & SIZE >= 100 ~ 100,
      1000 > SIZE & SIZE >= 500 ~ 500,
      1500 > SIZE & SIZE >= 1000 ~ 1000,
      2000 > SIZE & SIZE >= 1500 ~ 1500,
      2500 > SIZE & SIZE >= 2000 ~ 2000,
      3000 > SIZE & SIZE >= 2500 ~ 2500,
      SIZE >= 3000 ~ 3000
    )) %>%
    group_by(Range, VERDICT) %>%
    tally() %>%
    pivot_wider(names_from = VERDICT, id_cols = Range, values_from = n) %>%
    replace(is.na(.),0) %>%
    mutate(
      Precision = TP / (TP + FP),
      Recall = TP / (TP + FN),
      Range_end = case_when(
        Range < 0 & Range == -3000 ~ Range,
        Range < 0 & Range != -3000 & Range != -50 & Range != -100 ~ Range - 500,
        Range > 0 & Range != 3000 & Range != 50 & Range != 100 ~ Range + 500,
        Range > 0 & Range == 3000 ~ Range,
        Range > 0 & Range == 50 ~ Range + 50,
        Range > 0 & Range == 100 ~ Range + 400,
        Range < 0 & Range == -100 ~ Range - 400,
        Range < 0 & Range == -50 ~ Range - 50
      )) %>%
    replace(is.na(.),0)
  pm = mean(data2$Precision)
  pr = mean(data2$Recall)
  data2 %>%  mutate(sd_p = sqrt((Precision - pm)^2 / (FP + TP)),
      sd_r = sqrt((Recall - pr)^2 / (FN + TP))
    )
}
#Plot SV histogram length + PR
plot_histogram = function(data_hist_noinv,data_hist_segment,breaks=250,title){
  ggplot(data_hist_noinv%>%filter(SIZE > 49 | SIZE < -49),aes(x=SIZE)) +
    geom_segment(data=data_hist_segment%>%filter(Range != 0),aes(x=Range,y=Precision*500,xend=Range_end,yend=case_when(Range<0 & Range!=-3000~ lag(Precision),
                                                                                                  Range<0 & Range ==-3000 ~ Precision,
                                                                                                  Range>0 & Range != 3000 ~ lead(Precision),
                                                                                                  Range>0 & Range == 3000 ~ Precision)*500,colour="Precision"),size=2) +
    geom_segment(data=data_hist_segment%>%filter(Range != 0),aes(x=Range,y=Recall*500,xend=Range_end,yend=case_when(Range<0 & Range!=-3000~ lag(Recall),
                                                                                               Range<0 & Range ==-3000 ~ Recall,
                                                                                               Range>0 & Range != 3000 ~ lead(Recall),
                                                                                               Range>0 & Range == 3000 ~ Recall)*500,colour="Recall"),size=2) +
    geom_histogram(breaks=seq(-3000,3000,by=50)) +
    geom_ribbon(data=data_hist_segment%>%filter(Range != 0),aes(x=Range,ymin = Precision*500 - 1.96 * sd_p*500,
                                           ymax = Precision*500 + 1.96 * sd_p*500,fill="Precision"),alpha=0.2) +
    geom_ribbon(data=data_hist_segment%>%filter(Range != 0),aes(x=Range,ymin = Recall*500 - 1.96 * sd_r*500,
                                           ymax = Recall*500 + 1.96 * sd_r*500,fill="Recall"),alpha=0.2)+
    xlab("SV size (bp)") +
    theme_pubr()+
    scale_y_continuous("Counts",sec.axis = sec_axis(~ . / 5, name = "Values (%)")) +
    scale_x_continuous(breaks = seq(-3000,3000,by=breaks))+
    theme(axis.text.x=element_text(angle=45,hjust = 1,size = 40),legend.title=element_blank(),legend.key.width = unit(3,"cm"),
          legend.text = element_text(size = 40),axis.text = element_text(size = 40),axis.title = element_text(size = 40),axis.ticks = element_line(linewidth = 3),
          plot.title = element_text(size = 40))+
    guides(fill=FALSE) + scale_color_manual(values=c("lightblue","yellow3"))+
    scale_fill_manual(values = c("lightblue","yellow3")) +
    coord_cartesian(ylim = c(0,500),expand = FALSE) +
    ggtitle(title)
}
calculate_precision <- function(file_path,data) {
  # Read the .tsv file
  df <- read.csv(file_path, sep = " ", stringsAsFactors = FALSE,header = FALSE)
  colnames(df) = c("position","SVTYPE","SVSIZE","score")
  # Remove duplicates based on position, keeping the first occurrence
  df <- df %>%
    mutate(is_duplicate = position == lag(position)) %>%
    filter(is.na(is_duplicate) | !is_duplicate) %>%
    select(-is_duplicate)
  # Create a vector of thresholds from 0.001 to 1, incrementing by 0.001
  thresholds <- seq(0.001, 1, by = 0.001)
  
  # Initialize a data frame to store the results
  results <- data.frame(threshold = numeric(), precision = numeric())
  
  # Loop through each threshold
  for (threshold in thresholds) {
    # Filter the dataframe based on the threshold
    filtered_df <- df %>% filter(score >= threshold)
    
    # Calculate true positives and false positives
    true_positives <- nrow(filtered_df)
    false_positives <- nrow(df) - true_positives
    
    # Calculate precision
    precision <- true_positives / (true_positives + false_positives)
    
    # Append the result to the results data frame
    results <- rbind(results, data.frame(threshold = threshold, precision = precision,type=data))
  }
  
  # Return the results dataframe
  return(results)
}
###

#First benchmark: 
data_hist_segment = prepare_data(data_hist)
plot_hist_longass = plot_histogram(data_hist%>%filter(TYPE!="INV"),data_hist_segment,title="miniwalk-long-read against SVIM-ASM truth SV")

#Violin plot svimasm-assemblies PR
violin_longass = ggplot(pr_bench95,aes(x=Stat,y=Values*100)) +
  geom_violin(fill="#99CCFF") +
  geom_boxplot(fill="#99CCFF",width=0.1) +
  theme_pubr() +
  ylim(70,100) +
  ylab("Values (%)") +
  xlab("") +
  theme(axis.text = element_text(size = 40),axis.title = element_text(size = 40),plot.title = element_text(size = 40)) + ggtitle("SVIM-ASM truth")

##SHORT READS BENCHMARK
comparison_svimasm = comparison_plot(short_svimasm,title = "SVIM-ASM truth SV")
comparison = comparison_plot(short,title = "miniwalk-long-read truth SV")

#Graph Standard
#GRAPH
data_hist_segment_graph = prepare_data(data_hist_graph)
plot_graph = plot_histogram(data_hist_graph%>%filter(TYPE!="INV"),data_hist_segment_graph,breaks=500,title="miniwalk-short-read against miniwalk-long-read truth SV")

#MANTA
data_hist_segment_manta = prepare_data(data_hist_manta[abs(data_hist_manta$SIZE)>49,])
plot_manta = plot_histogram(data_hist_manta%>%filter(TYPE!="INV"),data_hist_segment_manta,breaks=500,title="manta against miniwalk-long-read truth SV")

#SVIM-ASM Standard
#GRAPH
data_hist_segment_graph_svimasm = prepare_data(data_hist_graph_svimasm)
plot_graph_svimasm = plot_histogram(data_hist_graph_svimasm%>%filter(TYPE!="INV"),data_hist_segment_graph_svimasm,breaks=500,title="miniwalk-short-read against SVIM-ASM truth SV")

#MANTA
data_hist_segment_manta_svimasm = prepare_data(data_hist_manta_svimasm[abs(data_hist_manta_svimasm$SIZE)>49,])
plot_manta_svimasm = plot_histogram(data_hist_manta_svimasm%>%filter(TYPE!="INV"),data_hist_segment_manta_svimasm,breaks=500,title="manta against SVIM-ASM truth SV")

#vcfdist ct score plot
vcfdist_long_p = calculate_precision("~/Documents/PHD/SV/SV_analysis/bench/vcfdist_longass_precision_NN.tsv",data="miniwalk-long-read")
vcfdist_long_r = calculate_precision("~/Documents/PHD/SV/SV_analysis/bench/vcfdist_longass_recall_NN.tsv",data="miniwalk-long-read")
colnames(vcfdist_long_r) = c("threshold","recall","type")
vcfdist = merge(vcfdist_long_p,vcfdist_long_r,by=c("threshold","type"))

vcfdist_graph_p = calculate_precision("~/Documents/PHD/SV/SV_analysis/bench/vcfdist_shortass_precision_NN.tsv",data="miniwalk-short-read")
vcfdist_graph_r = calculate_precision("~/Documents/PHD/SV/SV_analysis/bench/vcfdist_shortass_recall_NN.tsv",data="miniwalk-short-read")
colnames(vcfdist_graph_r) = c("threshold","recall","type")
vcfdist2 = merge(vcfdist_graph_p,vcfdist_graph_r,by=c("threshold","type"))
vcfdist = rbind(vcfdist,vcfdist2)

vcfdist_manta_p = calculate_precision("~/Documents/PHD/SV/SV_analysis/bench/vcfdist_manta_precision.tsv",data="manta")
vcfdist_manta_r = calculate_precision("~/Documents/PHD/SV/SV_analysis/bench/vcfdist_manta_recall.tsv",data="manta")
colnames(vcfdist_manta_r) = c("threshold","recall","type")
vcfdist2 = merge(vcfdist_manta_p,vcfdist_manta_r,by=c("threshold","type"))
vcfdist = rbind(vcfdist,vcfdist2)

vcfdist_p = ggplot(vcfdist,aes(threshold,precision,color=type)) +
  geom_line(size=3) +
  theme_pubr() +
  scale_color_manual(values = c("#0072B2","#E69F00","#66CC99")) +
  coord_cartesian(ylim = c(0,1.01),xlim = c(0,1.01),expand = FALSE) +
  scale_y_continuous(breaks = seq(0,1,by=0.1)) +
  xlab("Vcfdist partial credit threshold") +
  ylab("Precision") +
  theme(axis.text.x=element_text(angle=45,hjust = 1,size = 40),legend.title=element_blank(),legend.key.width = unit(3,"cm"),
        legend.text = element_text(size = 40),axis.text = element_text(size = 40),axis.title = element_text(size = 40),axis.ticks = element_line(linewidth = 3))
  
vcfdist_r = ggplot(vcfdist,aes(threshold,recall,color=type)) +
  geom_line(size=3) +
  theme_pubr() +
  scale_color_manual(values = c("#0072B2","#E69F00","#66CC99")) +
  coord_cartesian(ylim = c(0,1.01),xlim = c(0,1.01),expand = FALSE) +
  scale_y_continuous(breaks = seq(0,1,by=0.1)) +
  xlab("Vcfdist partial credit threshold") +
  ylab("Recall") +
  theme(axis.text.x=element_text(angle=45,hjust = 1,size = 40),legend.title=element_blank(),legend.key.width = unit(3,"cm"),
        legend.text = element_text(size = 40),axis.text = element_text(size = 40),axis.title = element_text(size = 40),axis.ticks = element_line(linewidth = 3))
blank <- grid::nullGrob()
#All plots in one grid!
row_2_3 = plot_grid(comparison_svimasm,plot_manta_svimasm,plot_graph_svimasm,blank,blank,blank,comparison,plot_manta,plot_graph,labels = c('c','d','e','','','','f','g','h'),rel_heights = c(1,0.1,1),nrow = 3,ncol = 3,label_size = 40)
row_1 = plot_grid(violin_longass,plot_hist_longass,labels = c('a','b'),ncol=2,nrow=1,label_size = 40)
row_4 = plot_grid(row_1,blank,row_2_3,ncol = 1,nrow = 3,rel_heights = c(1,0.1,1.75),rel_widths = c(1,1.5))
row_5 = plot_grid(vcfdist_p,vcfdist_r,labels=c('i','j'),ncol = 2,nrow = 1,label_size = 40)
final_benchmark_plot = plot_grid(row_4,blank,row_5,ncol = 1,nrow = 3,rel_heights = c(2,0.1,1))


png("~/Documents/PHD/SV/SV_analysis/figures/r_figures/benchmarking.png",height = 3000,width = 3300)
final_benchmark_plot
dev.off()

cov_pr = short_svimasm[short_svimasm$SV_Caller=="miniwalk-short-read" & short_svimasm$Stat=="Precision",]
cov_re = short_svimasm[short_svimasm$SV_Caller=="miniwalk-short-read" & short_svimasm$Stat=="Recall",]
cov_pr$Coverage = c(38,122,27,110,151,135,29,23,131,93,43,154,215,134,54,27)
cov_re$Coverage = c(38,122,27,110,151,135,29,23,131,93,43,154,215,134,54,27)
#Supplementary coverage analysis
cov_p = ggplot(cov_pr,aes(x=Coverage,y=Values)) +
  geom_point()+
  theme_pubr() +
  ylab("Precision") +
  geom_smooth(method = "lm")
cov_r = ggplot(cov_re,aes(x=Coverage,y=Values)) +
  geom_point()+
  theme_pubr()+
  ylab("Recall")+
  geom_smooth(method = "lm")

png("~/Documents/PHD/SV/SV_analysis/figures/r_figures/coverage_bench.png",height = 350,width = 570)
plot_grid(cov_p,cov_r,labels = c("a","b"),nrow=1,ncol=2,label_size = 20)
dev.off()


#Short read t tests
shapiro.test(short_svimasm[short_svimasm$Stat=="Precision" & short_svimasm$SV_Caller=="miniwalk-short-read",]$Values) #ALL NORMAL
t.test(short_svimasm[short_svimasm$Stat=="Precision" & short_svimasm$SV_Caller=="miniwalk-short-read",]$Values,short_manta[short_manta$Stat_manta=="Precision",]$Values_Manta)
t.test(short[short$Stat=="Recall" & short$SV_Caller=="Mtb-PRG",]$Values,short_manta[short_manta$Stat_manta=="Recall",]$Values_Manta)
t.test(short[short$Stat=="Precision" & short$SV_Caller=="Mtb-PRG",]$Values,short_manta[short_manta$Stat_manta=="Precision",]$Values_Manta)


#Clustering distance benchmark
clust_size = data.frame(Values=c(0.84,0.83,0.83,0.84,0.83,0.83,),
                        Stat=c("P","R","F","P","R","F","P","R","F","P","R","F","P","R","F"),
                        Size=c(0,0,0,50,50,50,100,100,100,150,150,150,200,200,200))
ggplot(clust_size,aes(x=Size,y=Values,color=Stat))+
  theme_pubr() +
  geom_line()
