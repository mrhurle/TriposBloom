library(tidytext)
library(tidyverse)
library(scales)
library(oce)
library(MBA)
library(lubridate)
library(reshape2)
library(ggplot2)

tripos<-read.csv("BiologicalData/masterdata.csv")

live<-tripos[c(1:3,6,9)]
avg.l<-tripos %>% #specify data frame
  group_by(Date,Depth,category) %>% #specify group indicator
  summarise_at(vars(Tripos.per.L), list(Tripos.p.L=mean))

dead<-tripos[c(1:3,20,23)]
colnames(dead)[4]<-"Tripos.per.L"
colnames(dead)[5]<-"category"

avg.d<-tripos%>%
  group_by(Date,Depth, category.4)%>%
  summarise_at(vars(DBP.per.L), list(Tripos.p.L=mean))
colnames(avg.d)[4]<-"Tripos.p.L"
colnames(avg.d)[3]<-"category"

tripos<-rbind(live,dead)
avg.tripos<-rbind(avg.l,avg.d)

tripos$Tripos.per.L<-tripos$Tripos.per.L+1
avg.tripos$Tripos.p.L<-avg.tripos$Tripos.p.L+1

tripos$Date <- factor(tripos$Date, levels = c("7-May-23","16-Jun-23","17-Jun-23","23-Jun-23","5-Jul-23","11-Jul-23","24-Jul-23","1-Aug-23","7-Aug-23","21-Aug-23","7-Sep-23"),labels = c("May 07", "June 16", "June 17","June 23","July 05","July 11","July 24","August 01","August 07","August 21","September 07"))
avg.tripos$Date <- factor(avg.tripos$Date, levels = c("7-May-23","16-Jun-23","17-Jun-23","23-Jun-23","5-Jul-23","11-Jul-23","24-Jul-23","1-Aug-23","7-Aug-23","21-Aug-23","7-Sep-23"),labels = c("May 07", "June 16", "June 17","June 23","July 05","July 11","July 24","August 01","August 07","August 21","September 07"))

tripos$category<-factor(tripos$category, levels = c("live", "dead"), labels = c("Live Cells","Dead Particles"))
avg.tripos$category<-factor(avg.tripos$category, levels = c("live", "dead"), labels = c("Live Cells","Dead Particles"))

t1<-ggplot(tripos,aes(x=Tripos.per.L, y=Depth, group=category, color=category))
t1 + geom_point(aes(fill = category), shape = 21, size = 1.3, stroke = 0.6, color = "black") + geom_path(data = avg.tripos, aes(x = Tripos.p.L, y = Depth, group = category, color = category)) + scale_fill_manual(values = c("#44AA99", "#757575"), name = "Category") + scale_color_manual(values = c("#44AA99", "#757575"), name = "Category") + scale_y_reverse() + scale_x_log10(labels = label_log(digits = 2), position = "top", breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5), limits = c(1e0, 1e5)) + facet_wrap(~Date, scales = "fixed") + labs(x = expression("Log Tripos per L"^-1), y = "Depth (m)") + theme_classic(base_size = 13) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.5), legend.position = c(0.88, 0.18), legend.title = element_blank(), legend.text = element_text(size = 10)) + guides(fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 2, stroke = 0.5)), color = guide_legend(override.aes = list(linetype = 1, size = 1)))


plank<-read.csv("BiologicalData/masterdata.csv")
pico<-plank[c(1:3,10,12)]
colnames(pico)[4]<-"Plank.per.mL"
colnames(pico)[5]<-"category"
avg.pico<- plank%>% #specify data frame
  group_by(Date,Depth,category.1) %>% #specify group indicator
  summarise_at(vars(Pico), list(Plankton.p.mL=mean))
avg.pico<-avg.pico[-c(51:54),]
colnames(avg.pico)[3]<-"category"

nano<-plank[c(1:3,13,15)]
colnames(nano)[4]<-"Plank.per.mL"
colnames(nano)[5]<-"category"
avg.nano<- plank%>% #specify data frame
  group_by(Date,Depth,category.2) %>% #specify group indicator
  summarise_at(vars(Nano), list(Plankton.p.mL=mean))
avg.nano<-avg.nano[-c(51:54),]
colnames(avg.nano)[3]<-"category"

syn<-plank[c(1:3,16,18)]
colnames(syn)[4]<-"Plank.per.mL"
colnames(syn)[5]<-"category"
avg.syn<- plank%>% #specify data frame
  group_by(Date,Depth,category.3) %>% #specify group indicator
  summarise_at(vars(Syn), list(Plankton.p.mL=mean))
avg.syn<-avg.syn[-c(51:54),]
colnames(avg.syn)[3]<-"category"

plank<-rbind(pico,nano,syn)
avg.plank<-rbind(avg.pico,avg.nano,avg.syn)

plank$Plank.per.mL<-plank$Plank.per.mL+1
avg.plank$Plankton.p.mL<-avg.plank$Plankton.p.mL+1

plank$Date <- factor(plank$Date, levels = c("7-May-23","16-Jun-23","17-Jun-23","23-Jun-23","5-Jul-23","11-Jul-23","24-Jul-23","1-Aug-23","7-Aug-23","21-Aug-23","7-Sep-23"),labels = c("May 07","June 16", "June 17","June 24","July 05","July 11","July 24","August 01","August 07","August 21","September 07"))
avg.plank$Date <- factor(avg.plank$Date, levels = c("7-May-23","16-Jun-23","17-Jun-23","23-Jun-23","5-Jul-23","11-Jul-23","24-Jul-23","1-Aug-23","7-Aug-23","21-Aug-23","7-Sep-23"),labels = c("May 07","June 16", "June 17","June 24","July 05","July 11","July 24","August 01","August 07","August 21","September 07"))

plank$category<-factor(plank$category, levels = c("Pico", "Nano", "Syn"), labels = c("Pico", "Nano", "Synechococcus"))
avg.plank$category<-factor(avg.plank$category, levels = c("Pico", "Nano", "Syn"), labels = c("Pico", "Nano", "Synechococcus"))

t2<-ggplot(plank,aes(x=Plank.per.mL, y=Depth, group=category, color=category))

t2 + geom_point(aes(fill = category), shape = 21, size = 1.3, stroke = 0.6, color = "black") + geom_path(data = avg.plank, aes(x = Plankton.p.mL, y = Depth, group = category, color = category))+ scale_y_reverse() + scale_x_log10(labels = label_log(digits = 2), position = "top", breaks = c(1e0, 1e1, 1e2, 1e3, 1e4, 1e5), limits = c(1e0, 1e5)) + facet_wrap(~Date, scales = "fixed") + labs(x = expression("Log Plankton per mL"^-1), y = "Depth (m)") + theme_classic(base_size = 13) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), panel.background = element_rect(fill = "white", color = NA), plot.background = element_rect(fill = "white", color = NA), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.5), legend.position = c(0.88, 0.18), legend.title = element_blank(), legend.text = element_text(size = 10)) + guides(fill = guide_legend(override.aes = list(shape = 21, color = "black", size = 2, stroke = 0.5)), color = guide_legend(override.aes = list(linetype = 1, size = 1)))


#cellratios
cellratio<-read.csv("BiologicalData/CellRatios.csv")
cellratio<-cellratio %>% #specify data frame
  group_by(Date,Depth) %>% #specify group indicator
  summarise_at(vars(Ratio), list(Avg.Ratio=mean))
logratio<-log10(cellratio$Avg.Ratio)
logratio<-as.data.frame(logratio)
cellratio<-cbind(cellratio,logratio)

cellratio$Date<-factor(cellratio$Date, levels = c("7-May-23","16-Jun-23","17-Jun-23","23-Jun-23","5-Jul-23","11-Jul-23","24-Jul-23","1-Aug-23","7-Aug-23","21-Aug-23","7-Sep-23") ,labels = c("May 07","June 16", "June 17","Jun 23","July 05","July 11","July 24","Aug 01","Aug 07","Aug 21","Sep 07"))
ratio<-ggplot(cellratio,aes(x=Date, y=Depth))

ratio + geom_point(aes(fill = logratio), shape = 21, color = "black", size = 3, stroke = 0.6) + scale_y_reverse() + labs(x = "", y = "Depth (m)", fill = "") + scale_fill_viridis_c() + scale_x_discrete(position = "top") + theme_classic(base_size = 13) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12), axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))


#nutrient plots
nuts<-read.csv("Oceanographic/nutrients.csv")
nuts$Date<-factor(nuts$Date, levels = c("23-Jun-23", "5-Jul-23", "11-Jul-23", "24-Jul-23", "1-Aug-23", "7-Aug-23", "7-Sep-23"), labels = c("June 23", "July 05", "July 11", "July 24", "Aug 01", "Aug 07", "Sep 07"))
ODV_colors <- c("#feb483", "#d31f2a", "#ffc000", "#27ab19", "#0db5e6", "#7139fe", "#d16cfa")

#nuts$Date<-as.Date(nuts$Date, format = "%d-%b-%Y")%>%as.Date(nuts$Date, format = "%d-%m-%Y")

#inorganic P
P<-nuts[,-c(2,4,6,7,8,9)]
P<-ggplot(P,aes(x=Date, y=Depth))

P + theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_point(aes(fill = Inorg.P), shape = 21, color = "black", size = 3, stroke = 0.6) + scale_y_reverse() + labs(x = "", y = "Depth (m)", fill = "") + scale_fill_viridis_c(limits = c(0.04, 0.17), breaks = c(0.05, 0.10, 0.15, 0.20)) + scale_x_discrete(position = "top")


#DOP
dop<-nuts[,-c(2,4,5,7,8,9)]
dop<-dop[-c(17),] #remove NA rows
dop<-ggplot(dop,aes(x=Date, y=Depth))


#NO3
no3<-nuts[,-c(2,4,5,6,7,9)]
no3<-ggplot(no3,aes(x=Date, y=Depth))
no3 + theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_point(aes(fill = NO3), shape = 21, color = "black", size = 3, stroke = 0.6) + scale_y_reverse() + labs(x = "", y = "Depth (m)", fill = "") + scale_fill_viridis_c() + scale_x_discrete(position = "top")



#DON
don<-nuts[,-c(2,4,5,6,7,8)]
don<-don[-c(9,15,29,37),] #remove NA rows
don<-ggplot(don,aes(x=Date, y=Depth))
don + theme_classic() + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12)) + geom_point(aes(fill = DON), shape = 21, color = "black", size = 3, stroke = 0.6) + scale_y_reverse() + labs(x = "", y = "Depth (m)", fill = "") + scale_fill_viridis_c() + scale_x_discrete(position = "top")


#import CTD data
CTD<-read.csv("Oceanographic/CTD.csv")
CTD$date<-as.Date(CTD$date, format = "%d-%b-%Y")%>%as.Date(CTD$date, format = "%d-%m-%Y")

#function to calculate MLD based on density
compute_mld <- function(depth, density, threshold = 0.01, max_ref_depth = 5) {if (length(depth) != length(density)) { stop("Depth and density vectors must be the same length.")}
  ord <- order(depth)
  depth <- depth[ord]
  density <- density[ord]
  ref_index <- which(depth <= max_ref_depth)
  if (length(ref_index) == 0) {stop("No depth measurements shallower than max_ref_depth.")}
  z_ref <- depth[max(ref_index)]
  rho_ref <- density[max(ref_index)]
  delta_rho <- density - rho_ref
  exceed_indices <- which(delta_rho >= threshold)
  if (length(exceed_indices) == 0) {return(max(depth))} else {return(depth[exceed_indices[1]])}}

den<-CTD[-c(383:433),] #remove 7/24 NA rows
den<-den[,-c(3,4,5,7,8,9,10)] #remove other physical paramerters 

mld_results <- den %>%
  group_by(date) %>%
  nest() %>%
  mutate(
    mld = map_dbl(data, ~ compute_mld(.x$depth, .x$density))
  ) %>%
  select(date, mld)

print(mld_results)

mld<-read.csv("Oceanographic/MLD.csv")
mld$Date<-as.Date(mld$Date, format = "%d-%b-%Y")%>%as.Date(mld$Date, format = "%d-%m-%Y")
#start_date <- as.Date("0023-06-23")
#end_date <- max(nuts$Date)  # Assuming temp_mba$date is your date range
#date_breaks <- seq.Date(from = start_date, to = end_date,length.out = 11)  # First, get weekly breaks
#date_breaks <- date_breaks[seq(1, length(date_breaks), by = 1)]  # Adjust breaks

#temperature ODV plot
temp<-CTD[,-c(4,5,6,7,8,9,10)]
temp$date<-as.numeric(temp$date) 
temp_mba<-mba.surf(temp,no.X = 800, no.Y = 800, extend = T) 
dimnames(temp_mba$xyz.est$z)<-list(temp_mba$xyz.est$x,temp_mba$xyz.est$y) 
temp_mba<-melt(temp_mba$xyz.est$z, varnames = c('date','depth'),value.name='temperature') 
temp_mba$date<-as.Date(temp_mba$date,"1970-01-01")

ggplot(data = temp_mba, aes(x = date, y = depth)) +geom_raster(aes(fill = temperature)) +geom_contour(aes(z = temperature), binwidth = 1.5, color = "black", alpha = 0.1) +scale_y_reverse(breaks = seq(10, 60, by = 10),expand = expansion(add = c(1, 1))) +scale_fill_gradientn(colors = rev(ODV_colors),breaks = c(8, 12, 16, 20)) +geom_line(data = mld, aes(x = Date, y = MLD), linetype = "dotted", color = "black") +geom_point(data = temp, aes(x = date, y = depth),color = "black", fill = "white",size = 1, shape = 21, stroke = 0.5) +labs(y = "Depth (m)", x = "", fill = "")+theme_classic()+theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white", color = NA),plot.background = element_rect(fill = "white", color = NA))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))+scale_x_date(labels = date_format("%b %d"),breaks = date_breaks("2 week"),position = "top",expand = expansion(add = c(3, 3)))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))

#salinity ODV plot
sal<-CTD[,-c(3,5,6,7,8,9,10)]
sal$date<-as.numeric(sal$date)
sal_mba<-mba.surf(sal,no.X = 800, no.Y = 800, extend = T)
dimnames(sal_mba$xyz.est$z)<-list(sal_mba$xyz.est$x,sal_mba$xyz.est$y)
sal_mba<-melt(sal_mba$xyz.est$z, varnames = c('date','depth'),value.name='salinity')
sal_mba$date<-as.Date(sal_mba$date,"1970-01-01")

ggplot(data = sal_mba, aes(x = date, y = depth)) +geom_raster(aes(fill = salinity)) +geom_contour(aes(z = salinity), binwidth = 1.5, color = "black", alpha = 0.1) +scale_y_reverse(breaks = seq(10, 60, by = 10),expand = expansion(add = c(1, 1))) +scale_fill_gradientn(colors = rev(ODV_colors)) +geom_line(data = mld, aes(x = Date, y = MLD), linetype = "dotted", color = "black") +geom_point(data = temp, aes(x = date, y = depth),color = "black", fill = "white",size = 1, shape = 21, stroke = 0.5) +labs(y = "Depth (m)", x = "", fill = "") +theme_classic() +theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white", color = NA),plot.background = element_rect(fill = "white", color = NA))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))+scale_x_date(labels = date_format("%b %d"),breaks = date_breaks("2 week"),position = "top",expand = expansion(add = c(3, 3)))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))

#fluorescence ODV plot
fluor<-CTD[,-c(3,4,5,6,8,9,10)]
fluor<-fluor[-c(383:433),] #remove 7/24 NA rows
fluor$fluorescence<-fluor$fluorescence+1
fluor$date<-as.numeric(fluor$date)
fluorlog<-log10(fluor$fluorescence)
fluorlog<-as.data.frame(fluorlog)
fluor<-cbind(fluor,fluorlog)
fluor<-fluor[,-c(3)]
fluor_mba<-mba.surf(fluor,no.X = 800, no.Y = 800, extend = T)
dimnames(fluor_mba$xyz.est$z)<-list(fluor_mba$xyz.est$x,fluor_mba$xyz.est$y)
fluor_mba<-melt(fluor_mba$xyz.est$z, varnames = c('date','depth'),value.name='fluorlog')
fluor_mba$date<-as.Date(fluor_mba$date,"1970-01-01")

ggplot(data = fluor_mba, aes(x = date, y = depth)) +geom_raster(aes(fill = fluorlog)) +geom_contour(aes(z = fluorlog), binwidth = 1.5, color = "black", alpha = 0.1) +scale_y_reverse(breaks = seq(10, 60, by = 10),expand = expansion(add = c(1, 1))) +scale_fill_gradientn(colors = rev(ODV_colors),breaks = c(0.25, 0.50, 0.75, 1)) +geom_line(data = mld, aes(x = Date, y = MLD), linetype = "dotted", color = "black") +geom_point(data = temp, aes(x = date, y = depth),color = "black", fill = "white",size = 1, shape = 21, stroke = 0.5) +labs(y = "Depth (m)", x = "", fill = "") +theme_classic() +theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white", color = NA),plot.background = element_rect(fill = "white", color = NA))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))+scale_x_date(labels = date_format("%b %d"),breaks = date_breaks("2 week"),position = "top",expand = expansion(add = c(3, 3)))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))

#par
#par<-CTD[,-c(3,4,5,6,7,8,9)]
#par<-par[-c(383:433),] #remove 7/24 NA rows
#par$PAR<-par$PAR+1
#par$date<-as.numeric(par$date)
#parlog<-log10(par$PAR)
#parlog<-as.data.frame(parlog)
#par<-cbind(par,parlog)
#par<-par[,-c(3)]
#par$date<-as.numeric(par$date)
#par_mba<-mba.surf(par,no.X = 800, no.Y = 800, extend = T)
#dimnames(par_mba$xyz.est$z)<-list(par_mba$xyz.est$x,par_mba$xyz.est$y)
#par_mba<-melt(par_mba$xyz.est$z, varnames = c('date','depth.m'),value.name='PAR')
#par_mba$date<-as.Date(par_mba$date,"1970-01-01")
#ggplot(data = par_mba,aes(x=date,y=depth.m))+geom_raster(aes(fill = PAR))+scale_fill_gradientn(colors=rev(ODV_colors))+geom_contour(aes(z = PAR), binwidth = 20, color = "black", alpha =0.2)+scale_y_reverse(breaks = seq(10, 60, by = 10))+scale_x_date(labels=date_format("%b-%d"), breaks=date_breaks("2 week"), position = "top") +coord_cartesian(expand = 0)+labs(y = "Depth (m)", x = "", fill = "")+theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12))+theme(axis.title.y=element_text(size=12))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))

#beam attenuation ODV
beamat<-CTD[,-c(3,4,5,6,7,9,10)]
beamat<-beamat[-c(383:433),] #remove 7/24 NA rows
beamat$date<-as.numeric(beamat$date)
beam_mba<-mba.surf(beamat,no.X = 800, no.Y = 800, extend = T)
dimnames(beam_mba$xyz.est$z)<-list(beam_mba$xyz.est$x,beam_mba$xyz.est$y)
beam_mba<-melt(beam_mba$xyz.est$z, varnames = c('date','depth'),value.name='beamAttenuation')
beam_mba$date<-as.Date(beam_mba$date,"1970-01-01")

ggplot(data = beam_mba, aes(x = date, y = depth)) +geom_raster(aes(fill = beamAttenuation)) +geom_contour(aes(z = beamAttenuation), binwidth = 1.5, color = "black", alpha = 0.1)+scale_y_reverse(breaks = seq(10, 60, by = 10),expand = expansion(add = c(1, 1))) +scale_fill_gradientn(colors = rev(ODV_colors),breaks=c(0.5,1.0,1.5,2)) +geom_line(data = mld, aes(x = Date, y = MLD), linetype = "dotted", color = "black") +geom_point(data = temp, aes(x = date, y = depth),color = "black", fill = "white",size = 1, shape = 21, stroke = 0.5) +labs(y = "Depth (m)", x = "", fill = "") +theme_classic() +theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white", color = NA),plot.background = element_rect(fill = "white", color = NA))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))+scale_x_date(labels = date_format("%b %d"),breaks = date_breaks("2 week"),position = "top",expand = expansion(add = c(3, 3)))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))

#beam transmission
beamtr<-CTD[,-c(3,4,5,6,7,8,10)]
beamtr<-beamtr[-c(383:433),] #remove 7/24 NA rows
beamtr$date<-as.numeric(beamtr$date)
beamt_mba<-mba.surf(beamtr,no.X = 800, no.Y = 800, extend = T)
dimnames(beamt_mba$xyz.est$z)<-list(beamt_mba$xyz.est$x,beamt_mba$xyz.est$y)
beamt_mba<-melt(beamt_mba$xyz.est$z, varnames = c('date','depth'),value.name='beamTransmission')
beamt_mba$date<-as.Date(beamt_mba$date,"1970-01-01")

ggplot(data = beamt_mba, aes(x = date, y = depth)) +geom_raster(aes(fill = beamtr)) +geom_contour(aes(z = salinity), binwidth = 1.5, color = "black", alpha = 0.1) +scale_y_reverse(breaks = seq(10, 60, by = 10),expand = expansion(add = c(1, 1))) +scale_fill_gradientn(colors = rev(ODV_colors),breaks = c(8, 12, 16, 20)) +geom_line(data = mld, aes(x = Date, y = MLD), linetype = "dotted", color = "white") +geom_point(data = temp, aes(x = date, y = depth),color = "black", fill = "white",size = 1, shape = 21, stroke = 0.5) +labs(y = "Depth (m)", x = "", fill = "") +theme_classic() +theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12),panel.background = element_rect(fill = "white", color = NA),plot.background = element_rect(fill = "white", color = NA))+guides(fill = guide_colorbar(barwidth = 1.5, barheight = 10))+scale_x_date(labels = date_format("%b %d"),breaks = date_breaks("2 week"),position = "top",expand = expansion(add = c(3, 3)))


#phytoclass -- pigment data
library(phytoclass)
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

library(ggplot2)
library(ggforce)

pigment+geom_bar(aes(), stat = "identity", position = "fill", width = 2) +geom_hline(yintercept = 0) +scale_fill_manual(values = rev(c("#44AA99", "#F5793A", "#882255", "#332288","#117733", "#88CCEE", "#DDCC77", "#AA4499"))) +scale_x_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = "")) +facet_nested(Date + Depth ~ .) +labs(x = "Relative Abundance (%)", y = "", fill = "Taxa") +guides(fill = guide_legend(nrow = 9, ncol = 1)) +theme_bw(base_size = 13) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),strip.text.y = element_text(size = 8, angle = 0),strip.background.y = element_rect(fill = "grey85", color = "black", linewidth = 0.5),legend.position = "right",panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),panel.background = element_rect(fill = "white", color = NA),plot.background = element_rect(fill = "white", color = NA))


avgpig <-avgpig %>% 
  group_by(names) %>% 
  summarise(vals = sum(vals, na.rm = TRUE))

sum<-sum(avgpig$vals)
avgpig$tot<-(avgpig$vals*100)/sum


pigabund$Date<-as.Date(pigabund$Date, format = "%d-%b-%y")%>%as.Date(pig.abund$Date, format = "%d-%m-%Y")
avgpig<-subset(pigabund,Date>=as.Date("2023-06-23")&Date<=as.Date("2023-07-11"))

