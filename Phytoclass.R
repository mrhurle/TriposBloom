library(phytoclass) 
library(ggplot2)     
library(ggforce)     
library(dplyr)       
library(lubridate)   
library(scales)

#Phytoclass HPLC/pigment data
pigment<-read.csv("BiologicalData/pigment.csv")
pig1<-pigment[-c(1,2)]
pig.metadata<-read.csv("BiologicalData/pigmetadata.csv")
Fmat<-read.csv("BiologicalData/fm.csv")
rownames(Fmat)<-Fmat[,1]
Fmat<-Fmat[-c(1)]
ratios<-read.csv("BiologicalData/values.csv")
pigment.results<-simulated_annealing(pig1,Fmat,ratios)
pigabund<-pigment.results$Figure$data
pigabund<-cbind(pig.metadata,pigabund)
pigabund<-pigabund[-c(3)]
write.csv(pigabund,file="pigabund.csv")

pigabund$Date <- factor(pigabund$Date, levels = c("23-Jun-23","5-Jul-23","11-Jul-23","24-Jul-23","1-Aug-23","7-Aug-23","7-Sep-23"),labels = c("June 23", "July 05", "July 11", "July 24", "August 01", "August 07", "September 07"))
pigabund$Depth<-factor(pigabund$Depth, levels = c("1.5", "2", "4", "5", "6", "10", "11", "13", "14", "15", "20", "30", "67"),labels = c("1.5m", "2m", "4m", "5m", "6m", "10m", "11m", "13m", "14m", "15m", "20m", "30m", "67m"))
pigabund$names<-factor(pigabund$names, levels = rev(c("Dinoflagellates", "Diatoms","Pras","Haptophytes", "Cryptophytes","Pelagophytes","Chlorophytes","Synechococcus")), labels = rev(c("Dinoflagellates", "Diatoms","Prasinophytes","Haptophytes", "Cryptophytes","Pelagophytes","Chlorophytes","Synechococcus")))

sum<-sum(pigabund$vals)
pigabund$relabund<-pigabund$vals/sum
pigabund$relabund<-pigabund$relabund*100

pigment<-ggplot(data=pigabund, aes(x=relabund, y="", group=names, fill=names))

pigment+geom_bar(aes(), stat = "identity", position = "fill", width = 2) +geom_hline(yintercept = 0) +scale_fill_manual(values = rev(c("#44AA99", "#F5793A", "#882255", "#332288","#117733", "#88CCEE", "#DDCC77", "#AA4499"))) +scale_x_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +facet_nested(Date + Depth ~ .) +labs(x = "Relative Abundance (%)", y = "", fill = "Taxa") +guides(fill = guide_legend(nrow = 9, ncol = 1)) +theme_bw(base_size = 13) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),strip.text.y = element_text(size = 8, angle = 0),strip.background.y = element_rect(fill = "grey85", color = "black", linewidth = 0.5),legend.position = "right",panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),panel.background = element_rect(fill = "white", color = NA),plot.background = element_rect(fill = "white", color = NA))
