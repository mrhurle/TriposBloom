library(phyloseq)
library(phyloseq.extended)
library(qiime2R)
library(vegan)
library(tidyverse)
library(fantaxtic)
library(ggh4x)
library(pairwiseAdonis)
library(dplyr)
library(tidyr)
library(tibble)

#load 18s taxonomy file - assigned with PR2 database
taxonomy<-read_qza("18S-rRNA/tax/asv-taxo.qza")
tax_tab<-taxonomy$data %>% #convert to data frame, tab separate and rename taxa levels, remove row with confidence values
  as.data.frame() %>%
  separate(Taxon, sep = ";", c("Domain","Supergroup", "Division", "Subdivision", "Class", "Order", "Family","Genus", "Species")) %>%
  column_to_rownames("Feature.ID") %>%
  dplyr::select(-Confidence)

#load ASV count table
table<-read_qza("18S-rRNA/dada2/table.qza")
count_tab<-table$data %>% as.data.frame() #convert to data frame

#load metadata file
metadata<-read.table("18S-rRNA/metadata/metadata18s.tsv", sep="\t", header=TRUE)

rownames(metadata)<-metadata$sample.id

#merge into phyloseq
ps<-phyloseq(tax_table(as.matrix(tax_tab)), otu_table(count_tab,taxa_are_rows = T),sample_data(metadata))
sample_sums(ps)

#remove unwanted groups 
ps_new=subset_taxa(ps, Supergroup!="Obazoa" |is.na(Supergroup))
sample_sums(ps_new)
ps_new=subset_taxa(ps_new, Division!="Proteobacteria" |is.na(Division))
sample_sums(ps_new)
ps_new=subset_taxa(ps_new, Domain!="Unassigned", Prune=T)

ps_new=name_na_taxa(ps_new) #adds an unassigned label to better identify lowest taxonomic level

sample_sums(ps_new)

#remove singletons (ASVs present once)
ps_filt_18S = filter_taxa(ps_new, function (x) {sum(x) > 1}, prune=TRUE)

#remove low counts (ASVs less than 5)
ps_filt_18S = filter_taxa(ps_filt_18S, function (x) {sum(x) > 5}, prune=TRUE)

#estimate minimum mean, and maximum reads counts
ps_min<-min(sample_sums(ps_filt_18S))
ps_mean<-mean(sample_sums(ps_filt_18S))
ps_max<-max(sample_sums(ps_filt_18S))
ps_sum<-sum(sample_sums(ps_filt_18S))

rare_18S <- ggrare(ps_filt_18S, step = 100, plot = FALSE, parallel = FALSE, se = FALSE)
rare_18S + theme_bw()

#rarefy to even sampling depth
ps_rare<-rarefy_even_depth(ps_filt_18S, sample.size = min(sample_sums(ps_filt_18S)), rngseed = 711, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)
sample_sums(ps_rare)

#rename ASVs in seqential order
taxa_names(ps_rare)<-paste0("ASV", seq(ntaxa(ps_rare)))

#export ASV information 
otu_filt=as(otu_table(ps_rare), "matrix")
tax_filt=as(tax_table(ps_rare), "matrix")
merge_filt<-cbind(otu_filt, tax_filt)
write.csv(merge_filt, file="Table_18sMerge.csv")

#GGPlot2
#relative abundance at class level~top10
ps_class<-tax_glom(ps_rare, taxrank = "Class") #agglomerate to class level
ps_class_trans<-transform_sample_counts(ps_class, function(x) 100*x/sum(x)) #transform to relative abundance
plot<- ps_class_trans %>%
  psmelt() %>% #melt data
  group_by(CollectionDate, Class, Depth) %>% #group by collection date and depth
  summarise_at("Abundance", .funs = mean)
focus <- c("Polycystinea", "Spirotrichea", "Embryophyceae", "Bacillariophyceae", "Cryptophyceae", "Prymnesiophyceae", "Syndiniales", "Mamiellophyceae", "Mediophyceae", "Dinophyceae") #focus on top 10 classes
plot$Class <- ifelse(plot$Class %in% focus, plot$Class, "Others") # Others category
plot_18S = plot # Rename data frame
plot_18S$Class<- as.character(plot_18S$Class) # Convert to character
p <- ggplot(data=plot_18S, aes(x=Abundance, y="", group=Class, fill=Class, color=Class))
p$data$Class <- factor(p$data$Class, levels = rev(c("Dinophyceae","Mediophyceae","Mamiellophyceae","Syndiniales", "Prymnesiophyceae","Cryptophyceae","Bacillariophyceae", "Embryophyceae","Spirotrichea","Polycystinea","Others"))) # Set order of groups in the plot
p$data$CollectionDate <- factor(p$data$CollectionDate, levels = c("23-Jun-23","5-Jul-23", "11-Jul-23","24-Jul-23", "1-Aug-23", "7-Aug-23"),labels = c("June 23","July 05", "July 11","July 24", "August 01", "August 07")) # Set order of transects
p$data$Depth <-factor(p$data$Depth,levels = c("1.5","2","4","5","6","10","13","14","15"), labels = c("1.5m","2m","4m","5m","6m","10m","13m","14m","15m")) #change facet labels

p +geom_bar(aes(), stat = "identity", position = "fill", width = 1.5) +scale_x_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""))+scale_fill_manual(values = rev(c("#44AA99", "#F5793A", "#882255", "#88CCEE", "#332288","#117733", "#DDCC77", "#AA4499", "#F7CDA4", "#A5CFCC", "#757575"))) +scale_color_manual(values = rev(c("#44AA99", "#F5793A", "#882255", "#88CCEE", "#332288","#117733", "#DDCC77", "#AA4499", "#F7CDA4", "#A5CFCC", "#757575"))) +facet_nested(CollectionDate + Depth ~ .) +labs(x = "Relative Abundance (%)", y = "") +guides(fill = guide_legend(nrow = 11, ncol = 1)) +theme_bw(base_size = 13) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),strip.text.y = element_text(size = 9, angle = 0),legend.position = "right",panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

#relative abundance plot within class Dinophyceae~by family
ps_dino=subset_taxa(ps_rare, Class=="Dinophyceae", prune = T) #subset class dinophyceae
dino_fam<-tax_glom(ps_dino, taxrank = "Family") #agglomerate to family level
ps_dino_trans<-transform_sample_counts(dino_fam, function(x) 100*x/sum(x)) #transform to relative abundance
plot1<- ps_dino_trans %>%
  psmelt() %>% #melt data
  group_by(CollectionDate, Family,Depth) %>% #group by collection date and depth
  summarise_at("Abundance", .funs = mean)
focus1 <- c("Ceratiaceae","Suessiaceae","Gymnodiniaceae","Dinophysiaceae","Pyrocystaceae","Kareniaceae","Prorocentraceae","Heterocapsaceae") #focus on top 8 Family
plot1$Family <- ifelse(plot1$Family %in% focus1, plot1$Family, "Others") # Others category
plot_18S_Dino = plot1 # Rename data frame
plot_18S_Dino$Family<- as.character(plot_18S_Dino$Family) # Convert to character
p1 <- ggplot(data=plot_18S_Dino, aes(x=Abundance, y="", group=Family, fill=Family, color=Family))
p1$data$Family <- factor(p1$data$Family, levels = rev(c("Ceratiaceae","Suessiaceae","Gymnodiniaceae","Dinophysiaceae","Pyrocystaceae","Kareniaceae","Prorocentraceae","Heterocapsaceae","Others"))) # Set order of groups in the plot
p1$data$CollectionDate <- factor(p1$data$CollectionDate,levels = c("23-Jun-23","5-Jul-23", "11-Jul-23","24-Jul-23", "1-Aug-23", "7-Aug-23"),labels = c("June 23","July 05", "July 11","July 24", "August 01", "August 07")) # Set order of transects
p1$data$Depth <-factor(p1$data$Depth,levels = c("1.5","2","4","5","6","10","13","14","15"), labels = c("1.5m","2m","4m","5m","6m","10m","13m","14m","15m")) #change facet labels
p1+geom_bar(aes(), stat = "identity", position = "fill", width = 1.5) +scale_x_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""))+scale_fill_manual(values = rev(c("#44AA99", "#F5793A", "#882255", "#88CCEE", "#332288","#117733", "#DDCC77", "#AA4499", "#757575"))) +scale_color_manual(values = rev(c("#44AA99", "#F5793A", "#882255", "#88CCEE", "#332288","#117733", "#DDCC77", "#AA4499", "#757575"))) +facet_nested(CollectionDate + Depth ~ .) +labs(x = "Relative Abundance (%)", y = "") +guides(fill = guide_legend(nrow = 11, ncol = 1)) +theme_bw(base_size = 13) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),strip.text.y = element_text(size = 9, angle = 0),legend.position = "right",panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

#relative abundance within class Dinophyceae~genus~top5
ps_genus<-tax_glom(ps_dino, taxrank = "Genus") #agglomerate to genus level
ps_genus_trans<-transform_sample_counts(ps_genus, function(x) 100*x/sum(x)) #transform to relative abundance
plot2<- ps_genus_trans %>%
  psmelt() %>% #melt data
  group_by(CollectionDate, Genus, Depth) %>% #group by collection date and depth
  summarise_at("Abundance", .funs = mean)
focus2 <- c("Tripos","Biecheleriopsis","Pelagodinium","Biecheleria","Dinophysis") #focus on top 5 genus
plot2$Genus <- ifelse(plot2$Genus %in% focus2, plot2$Genus, "Others") # Others category
plot_18S_Genus = plot2 # Rename data frame
plot_18S_Genus$Genus<- as.character(plot_18S_Genus$Genus) # Convert to character
p2 <- ggplot(data=plot_18S_Genus, aes(x=Abundance, y="", group=Genus, fill=Genus, color=Genus))
p2$data$Genus <- factor(p2$data$Genus, levels = rev(c("Tripos","Biecheleriopsis","Pelagodinium","Biecheleria","Dinophysis","Others"))) # Set order of groups in the plot
p2$data$CollectionDate <- factor(p2$data$CollectionDate,levels = c("23-Jun-23","5-Jul-23", "11-Jul-23","24-Jul-23", "1-Aug-23", "7-Aug-23"),labels = c("June 23","July 05", "July 11","July 24", "August 01", "August 07")) # Set order of transects
p2$data$Depth <-factor(p2$data$Depth,levels = c("1.5","2","4","5","6","10","13","14","15"), labels = c("1.5m","2m","4m","5m","6m","10m","13m","14m","15m")) #change facet labels
p2+geom_bar(aes(), stat = "identity", position = "fill", width = 1.5) +scale_x_continuous(labels = scales::label_percent(scale = 100, prefix = "", suffix = ""))+scale_fill_manual(values = rev(c("#44AA99", "#F5793A", "#882255", "#88CCEE", "#332288", "#757575"))) +scale_color_manual(values = rev(c("#44AA99", "#F5793A", "#882255", "#88CCEE", "#332288", "#757575"))) +facet_nested(CollectionDate + Depth ~ .) +labs(x = "Relative Abundance (%)", y = "") +guides(fill = guide_legend(nrow = 11, ncol = 1)) +theme_bw(base_size = 13) +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1),axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),strip.text.y = element_text(size = 9, angle = 0),legend.position = "right",panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))

#Beta diversity for whole community grouped by date/depth
ordu=ordinate(ps_rare, "PCoA", "bray") #PCoA ordination based on Bray-Curtis
pord=plot_ordination(ps_rare, ordu, color = "CollectionDate")
pord$data$CollectionDate<-as.factor(pord$data$CollectionDate) #convert collection.date to factor
#pord$data$WaterCol<-as.factor(pord$data$WaterCol) #set as factor
pord$data$CollectionDate <- factor(pord$data$CollectionDate,levels = c("23-Jun-23","5-Jul-23", "11-Jul-23","24-Jul-23", "1-Aug-23", "7-Aug-23"),labels = c("June 23","July 05", "July 11","July 24", "August 01", "August 07"))#set order of collection dates
pord$data$WaterCol <-factor(pord$data$WaterCol,levels = c("surface", "mid", "deep"), labels = c("Surface", "Mid", "Deep")) #set order of water column depths

pord+geom_point(aes(color = CollectionDate, fill = CollectionDate, shape = WaterCol),size = 2.5,stroke = 0.8,alpha = 0.9,color = "black") +scale_shape_manual(values = c("Surface" = 21, "Mid" = 22, "Deep" = 23)) +theme_bw() +theme(legend.position = "right",panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),,panel.border = element_rect(color = "black", fill = NA, size = 1),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +labs(color = "Date", shape = "Depth", fill = "Date") +guides(color = guide_legend(order = 1),fill = "none")

#Pairwise PERMANOVA ~ test significance of collection date (or depth)
metadata<-as(sample_data(ps_rare), "data.frame")
metadata$CollectionDate=as.factor(metadata$CollectionDate)
bray=phyloseq::distance(ps_rare, method = "bray")
adonis2(phyloseq::distance(ps_rare, method = "bray")~CollectionDate,data = metadata, permutations = 9999)
pairwise.adonis(as.dist(bray), as.factor(metadata$CollectionDate),p.adjust.m = 'bonferroni', perm = 9999)

#alpha diversity between collection dates (or depth)
ps_rich<-estimate_richness(ps_rare, measures = c("Observed", "Shannon")) #estimate richness and Shannon diversity index
Date=sample_data(ps_rare)$CollectionDate #define different collection dates
Depth=sample_data(ps_rare)$WaterCol
rich_all<-data.frame(ps_rich,Date,Depth) #merge diversity values with collection date
rich_all$Date=as.factor(rich_all$Date)
ad<-ggplot(rich_all, aes(x=Date, y=Shannon, color=Date, shape =Depth ))
ad$data$Date<-factor(ad$data$Date,levels = c("23-Jun-23","5-Jul-23", "11-Jul-23","24-Jul-23", "1-Aug-23", "7-Aug-23"),labels = c("June 23","July 05", "July 11","July 24", "August 01", "August 07")) #order dates
ad$data$Depth<-factor(ad$data$Depth,levels = c("surface", "mid", "deep"), labels = c("Surface", "Mid", "Deep"))
#point plot ~ by date and water column depth

ad +ylab("Shannon Diversity Index") +theme(legend.position = "right",axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),panel.border = element_rect(color = "black", fill = NA, size = 1),plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +geom_point(aes(fill = Date, shape = Depth),color = "black",size = 2.5,stroke = 0.8) +scale_shape_manual(values = c("Surface" = 21, "Mid" = 22, "Deep" = 23))+guides(fill = guide_legend(title = "Date",override.aes = list(color = "black", fill = scales::hue_pal()(6), shape = 21)),shape = guide_legend(title = "Depth"))

#ANOVA ~ check if variable of date affects Shannon diversity
anova.sh=aov(ps_rich$Shannon ~ sample_data(ps_rare)$CollectionDate)
summary(anova.sh)
TukeyHSD(anova.sh)

# ---------------------------
# 1. Bray-Curtis distance
# ---------------------------
ps_rel <- transform_sample_counts(ps_rare, function(x) 100 * x / sum(x))

bray_dist <- phyloseq::distance(ps_rel, method = "bray")


# ---------------------------
# 2. Environmental data (aligned ONLY samples)
# ---------------------------
metadata <- as(sample_data(ps_rel), "data.frame")

env <- metadata[, c("NO3", "PO4", "DON", "DOP",
                    "Salinity", "Temp", "Beam.Attenuation")]

# IMPORTANT: remove rows with NA BEFORE analysis
keep <- complete.cases(env)

env_clean <- env[keep, ]
bray_clean <- as.matrix(bray_dist)[keep, keep]


# ---------------------------
# 3. Z-score standardization (REVIEWER REQUIREMENT)
# ---------------------------
env_clean <- log1p(env_clean)
env_clean <- as.data.frame(scale(env_clean))


# ---------------------------
# 4. dbRDA model (NO stepwise selection)
# ---------------------------
dbrda_final <- dbrda(bray_clean ~ NO3 + PO4 + DON + DOP +
                       Salinity + Temp,
                     data = env_clean)


# ---------------------------
# 5. Model significance
# ---------------------------
anova(dbrda_final, permutations = 9999)
anova(dbrda_final, by = "term", permutations = 9999)

#adjusted R2 for full dbRDA model
RsquareAdj(dbrda_final)

# ---------------------------
# 6. VIF check (must be reported in manuscript)
# ---------------------------
vif_vals <- vif.cca(dbrda_final)
print(vif_vals)


# ---------------------------
# 7. Variance explained
# ---------------------------
eig_vals <- summary(dbrda_final)$cont$importance
var1 <- round(eig_vals[2, 1] * 100, 1)
var2 <- round(eig_vals[2, 2] * 100, 1)


# ---------------------------
# 8. Site scores (FIXED ALIGNMENT)
# ---------------------------
site_scores <- scores(dbrda_final, display = "sites")
site_scores <- as.data.frame(site_scores)

# align metadata AFTER filtering
meta_clean <- metadata[keep, ]

site_scores$CollectionDate <- meta_clean$CollectionDate
site_scores$Depth <- meta_clean$WaterCol


# ---------------------------
# 9. Environmental vectors
# ---------------------------
env_scores <- scores(dbrda_final, display = "bp")
arrowdf <- as.data.frame(env_scores)
arrowdf$labels <- rownames(arrowdf)

# Remove DOP vector and label
arrowdf <- subset(arrowdf, labels != "DOP")

arrow_map <- aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2)
label_map <- aes(x = dbRDA1 * 1.15, y = dbRDA2 * 1.15, label = labels)


# ---------------------------
# 10. PLOT (keep your visualization, but correct stats interpretation)
# ---------------------------

site_scores$CollectionDate <- factor(
  site_scores$CollectionDate,
  levels = c("23-Jun-23","5-Jul-23","11-Jul-23","24-Jul-23","1-Aug-23","7-Aug-23"),
  labels = c("June 23","July 05","July 11","July 24","August 01","August 07")
)

site_scores$Depth <- factor(
  site_scores$Depth,
  levels = c("surface","mid","deep"),
  labels = c("Surface","Mid","Deep")
)


ggplot(site_scores, aes(x = dbRDA1, y = dbRDA2)) +
  geom_point(aes(color = CollectionDate, shape = Depth),
             size = 2.5, stroke = 0.8) +
  geom_vline(xintercept = 0, color = "grey70", linetype = 2) +
  geom_hline(yintercept = 0, color = "grey70", linetype = 2) +
  geom_segment(data = arrowdf, mapping = arrow_map,
               size = 0.3,
               arrow = arrow(length = unit(0.01, "npc"))) +
  geom_text(data = arrowdf, mapping = label_map, size = 3) +
  theme_bw() +
  labs(
    color = "Date",
    shape = "Depth",
    x = paste0("dbRDA1 (", var1, "%)"),
    y = paste0("dbRDA2 (", var2, "%)")
  )+ theme_classic()

