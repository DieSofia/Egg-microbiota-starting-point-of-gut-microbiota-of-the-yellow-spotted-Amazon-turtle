#########################################################################################################################################################
R-script used for the publication "Egg microbiota is the starting point of hatchling gut microbiota in the endangered yellow-spotted Amazon river turtle" in the Molecular Ecology Journal by Carranco et al.,2022
################################################################################################################## Amazon freshwater turtle Project ############################################################################################################################################################################################################
#Created by Ana Sofia Carranco, Department of Evolutionary Ecology and Conservation Genomics, Ulm University, Germany
# For questions or requests contact me to: ana.carranco-narvaez@uni-ulm.de
########################################################################################
# Sabías que estes donde estes, el aire que respiramos y el agua que tomamos viene en gran parte de los bosques Amazónicos? Increíble verdad? Mi deseo es que siempre tengamos ese privilegio.
########################################################################################


#R version 3.6.3 (2020-02-29)
#load all the necessary packages

library('BiocManager')
library('devtools')
library('phyloseq')
library('microbiome')
library('gtable')
library('gridExtra')
library('nloptr')
library('ggpubr')
library('MASS')
library('MuMIn')
library('R2admb')
library('glmmADMB')
library('glmmTMB')
library('lme4.0')
library('tidyverse')
library('mgcv')
library('gratia')
library('R2admb')
library('multcomp')
library('vegan')
library('pairwiseAdonis')
library('cluster')
library('cowplot')
library('devtools')
library('MASS')
library('microbiomeutilities')
library('stringr')
library('phyloseq')
library('devtools')
library('stringr')
library('data.table')
library('ggplot2')
library('plyr')
library('dplyr')
library('data.table') 
library('ggsci')
library('grid')
library('emmeans')

#import the phyloseq which contains all the data used for this study, these data are already filtered for chimeras, chloroplast ASVs and mitochondria, and are clean of contaminants that may come from the controls.

#Read the RDS file
punifilis_pure <- readRDS("punifilis_pure.rds")


#Here we subset our Blank samples in a different object becasue after this step we remove samples with less than 10 000 reads (and our Blanks have very low number of reads but we want to keep them for the analyses)

punifilis_pure_control =subset_samples(punifilis_pure, SampleForm=="Blank")


punifilis_pure_control
#We will merge this object with all the turtle data after a couple of steps

### Here we remove samples with less than 10 000 reads from our data set
punifilis_pure =subset_samples(punifilis_pure, sample_sums(punifilis_pure)>10000)# check the threshhold with 5000 and add the environmentla samples 

## Here we have our object, which basically contains all the data from eggs, turtles' cloaca and environmental samples for analysis (N=195)
summarize_phyloseq(punifilis_pure)

# However we want to add Blank samples to our analysis. So here we merge the two datasets and create an object with the biological dataset and the blanks
punifilis_control<-merge_phyloseq(punifilis_pure_control,punifilis_pure)

## Here we have our object, with all the data (eggs, turtles' cloaca, environment, and blank samples) for analysis (N=203)
punifilis_control


#################################################################

						      ###### ALPHA DIVERSITY ######	

#################################################################


Purenifilis_alpha <- alpha(punifilis_control, index = c("shannon","observed","simpson"))
data.table(Purenifilis_alpha)

Purenifilis_alpha_meta<- meta(punifilis_control)
Purenifilis_alpha_meta$Featureid2<-rownames(Purenifilis_alpha_meta)
Purenifilis_alpha$Featureid2<-rownames(Purenifilis_alpha)
Purenifilis_alfadiv<- merge(Purenifilis_alpha_meta,Purenifilis_alpha, by="Featureid2")

colnames(Purenifilis_alfadiv)

head(sample_data(Purenifilis_alpha, 100))

### Some color pallets that I used for the plots
ana_co <- c("#15542e","#ffe734","#abc100","#ffa951","#01a39c","#cc73e7","#668578","#0083fb","#00e098","#5dc7ff","#34ea62","#cedfff","#51426b","#fff0c0","#8cffd6")

ana_col_alot <- c("#250034","#51b222","#2234b2","#55d876","#8a3ec3","#00b357","#b934b7","#368500","#cf7aff","#00670a","#006cee","#f68117","#00a3fd","#c5590c","#0181d0","#fdb046","#005caa","#8d8e00","#65006b","#c2c55c","#3d0050","#d9be51","#003572","#b57600","#93a4ff","#5a7800","#e275d8","#727e00","#984b91","#76d1a5","#0e001e","#02d2cb","#834d00","#8abcff","#4f3500","#58ceee","#322000","#cdb4f7","#002403","#e89362","#01a4d7","#ad6f4b","#003962","#deba85","#00171b","#c3c28d","#1c1500","#96cbab","#71436c","#00693b","#a37c9e","#1c2600","#9dc6d1","#004b3f","#bdc0c0","#004456","#01947b","#806251","#017788","#016c51")

ana_co5<- c("#0e5200","#bd8c10","#6164bd","#009f22","#0189fb","#a47e91","#f7671f","#003e8f","#c17cc3","#01a0b8")

ana_co5.2<- c("#6164bd","#bd8c10","#0e5200","#009f22","#0189fb","#a47e91","#f7671f","#003e8f","#c17cc3","#01a0b8")

ana_co5.3<- c("#405461","#aaaa37","#b90f73","#c17cc3","#00a436","#cf8935","#6b15b0","#01fac4","#f7671f","#4c2d23","#6479cc","#000000","#a37c9e")


			#############################################
			########## ALPHA DIVERSITY PLOTS ############
			#############################################


# Here we create a new column called Stage, in this column we put together the cloaca samples from day 30 and over day 30 of development 

levels(Purenifilis_alfadiv$Age)

Purenifilis_alfadiv$Stage<-factor(Purenifilis_alfadiv$Age,levels=c("Embryo","Egg","0","5","10","15","20","25","30","42","44","51","61","62","EnvSand","EnvWater","ExtBlank", "FieldBlank"),labels=c("UnhatchedEgg","HatchedEgg", "day0", "day5", "day10", "day15", "day20", "day25", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater","NestEnvironment","RiverEnvironment","Blank","Blank"))

# Here we create a new column to be used for the labels in the graphs
levels(Purenifilis_alfadiv$Stage)

Purenifilis_alfadiv$SampleFrom<-factor((Purenifilis_alfadiv$Stage), levels = c("UnhatchedEgg","HatchedEgg", "day0", "day5", "day10", "day15", "day20", "day25", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater","NestEnvironment","RiverEnvironment","Blank"), labels= c("Inner eggshell","Inner eggshell","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Env","Env","Blank"))






				###################################################
				############# PLOT NO.1 OBSERVED DIVERSITY ########
				#############      (Species Richness)      ########
				###################################################


ggplot(Purenifilis_alfadiv,aes(x = Stage, y = observed,fill= SampleForm))+geom_boxplot(outlier.shape=NA)+theme_linedraw()+geom_jitter(colour= "black", position=position_dodge(0.8), size=1.5)+xlab(" ")+ylab("Species Richness\n")+guides(fill=guide_legend(title = element_blank(), position = "right", reverse = TRUE,keywith = 1, keyheight = 1),shape=guide_legend(title = waiver(), position = "right", reverse = TRUE, keywith = 1, keyheight = 1))+scale_fill_manual(values= ana_co5,breaks=c("Egg","Turtle","Env"),labels=c("Egg","Turtle","Environment"))+facet_grid(~SampleFrom, scales = "free", space = "free")+theme(strip.text.x = element_text(size = 16))+theme(axis.text=element_text(size = 14))+theme(axis.title=element_text(size = 16))+scale_x_discrete(labels = function(x) str_wrap(x, width = 18))+ theme(legend.text = element_text(size = 16))+theme(panel.grid.minor = element_line(size = 0.05), panel.grid.major = element_line(size = 0.05))+scale_x_discrete(labels=c("UnhatchedEgg" = "L.D.","HatchedEgg" = "Hatched","day0" = "Day 0", "day5" = "Day 5", "day10" = "Day 10", "day15" = "Day 15", "day20" = "Day 20", "day25" = "Day 25", "day30andGreater" = "Day \u226530","NestEnvironment" = "Nest", "RiverEnvironment" = "River", "Blank" = "Blank"))

		
				###################################################
				############# PLOT NO.2 Shannon Index #############
				###################################################



ggplot(Purenifilis_alfadiv,aes(x = Stage, y = diversity_shannon,fill= SampleForm))+geom_boxplot(outlier.shape=NA)+theme_linedraw()+geom_jitter(colour= "black", position=position_dodge(0.8), size=1.5)+xlab(" ")+ylab("Shannon diversity\n")+guides(fill=guide_legend(title = element_blank(), position = "right", reverse = TRUE,keywith = 1, keyheight = 1),shape=guide_legend(title = waiver(), position = "right", reverse = TRUE, keywith = 1, keyheight = 1))+scale_fill_manual(values= ana_co5,breaks=c("Egg","Turtle","Env"),labels=c("Egg","Turtle","Environment"))+facet_grid(~SampleFrom, scales = "free", space = "free")+theme(strip.text.x = element_text(size = 16))+theme(axis.text=element_text(size = 14))+theme(axis.title=element_text(size = 16))+scale_x_discrete(labels = function(x) str_wrap(x, width = 18))+ theme(legend.text = element_text(size = 16))+theme(panel.grid.minor = element_line(size = 0.05), panel.grid.major = element_line(size = 0.05))+ scale_x_discrete(labels=c("UnhatchedEgg" = "L.D.","HatchedEgg" = "Hatched","day0" = "Day 0", "day5" = "Day 5", "day10" = "Day 10", "day15" = "Day 15", "day20" = "Day 20", "day25" = "Day 25", "day30andGreater" = "Day \u226530","NestEnvironment" = "Nest", "RiverEnvironment" = "River", "Blank" = "Blank"))



##############################   ALPHA DIVERSITY ANALYSES  ###################################



# We calculate Sequencing depth
SeqDepth<-as.data.frame(sample_sums(punifilis_control))
SeqDepth<-setDT(SeqDepth, keep.rownames = TRUE)[]
	head(SeqDepth)
	names(SeqDepth)[2] <-"SeqDepth"
	names(SeqDepth)[1] <-"Featureid2"
	head(SeqDepth)

# We merge the table holding the Sequencing depth information to our main data table
Purenifilis_alfadiv<-merge(Purenifilis_alfadiv,SeqDepth,by="Featureid2")

DevAlphaEnv<-Purenifilis_alfadiv

########## Observed diverstiy (Species Richness) ##########

M1Dev<-gam(observed~Stage+s(scale(SeqDepth)),data=DevAlphaEnv,na.action=na.fail,family=nb,method="ML")

#Lets have a look to the results
summary(M1Dev)

#Lets have a look to our model selection
dredge(M1Dev, extra="adjR^2")


###### Multiple comparison post-hoc test for Species Richness
posthocOptimal<-summary(glht(M1Dev))
posthocOptimal$test$sigma
posthocOptimal$test

#Here we have a table with all the multiple comparision calculations and the meta data. Just load the table for the plot
posthocOptimalDF<-read.table("posthocOptimalDF.txt", header=TRUE, sep="")

## For our plot we add the significance to the p-values
posthocOptimalDF$Sign<-ifelse(posthocOptimalDF$Pvalue < 0.001, "***", ifelse(posthocOptimalDF$Pvalue < 0.01, "**",ifelse(posthocOptimalDF$Pvalue < 0.05, "*","")))

ggplot(posthocOptimalDF, aes(reorder(StageFig, Coefficient), Coefficient,label = Sign)) +geom_hline(yintercept=0, colour="#8C2318", size=1) +   geom_pointrange(aes(ymin=lwr, ymax=upr)) +  labs(x="Sample type comparison", y="Estimate", title="Species diversity") + coord_flip() + theme_bw()+geom_text(nudge_y=0.35,nudge_x=0.18, size= 5)+theme(axis.text=element_text(size = 10))+theme(axis.title=element_text(size = 14))+theme(title=element_text(size = 10))

########## Shannon diversity ##########

M1DevSh<-glm(diversity_shannon~Stage+scale(SeqDepth), family=Gamma(log),data=DevAlphaEnv,na.action=na.fail)

#Lets have a look to the results
summary(M1DevSh)

dredge(M1DevSh, extra="adjR^2")

#Multiple comparison post-hoc Shannon diversity

posthocOptimal<-summary(glht(M1DevSh))
posthocOptimal$test$sigma
posthocOptimal$test

## Read table
read.table(posthocOptimalShDF,"/home/student-gen-a/Dropbox/SharedFolder/MolEcolSecondSub/Tables/posthocOptimalShDF.txt")

# we add the significance to the p-values
posthocOptimalShDF$Sign<-ifelse(posthocOptimalShDF$Pvalue < 0.001, "***", ifelse(posthocOptimalShDF$Pvalue < 0.01, "**", ifelse(posthocOptimalShDF$Pvalue < 0.05, "*",""))) 

ggplot(posthocOptimalShDF, aes(reorder(StageFigSh, Coefficient), Coefficient,label = Sign)) +geom_hline(yintercept=0, colour="#8C2318", size=1) +   geom_pointrange(aes(ymin=lwr, ymax=upr)) +  labs(x="Sample type comparison", y="Estimate", title="Shannon diversity") + coord_flip() + theme_bw()+geom_text(nudge_y=0.14,nudge_x=0.10, size= 5)+theme(axis.text=element_text(size = 10))+theme(axis.title=element_text(size = 14))+theme(title=element_text(size = 10))



##############################################
##########  BETA DIVERSITY ANALISYS ########## 
##############################################

###############Compositional##############

############### Object with Blanks ##############

#Taxa filtering
punifilis_control_comp <-  phyloseq::filter_taxa(punifilis_control, function(x) sum(x>0) > 30, TRUE)

# load OTU table
tableclrc<-t(otu_table(punifilis_control_comp))

# call the clr function
clr_coordinates <- function(X) {
  lX = log(X)
  lX - apply(lX, 1, mean)
} 


######### for data with blanks #######
## integrate the clr function to the OTU table and load the necessary data
tableclrc <- clr_coordinates(tableclrc+1)
data2c<-sample_data(punifilis_control)
taxaclrc<-tax_table(punifilis_control)
treeclrc<-phy_tree(punifilis_control)

## we merge our data in a single phyloseq
clrc <- merge_phyloseq(tableclrc,data2c,taxaclrc,treeclrc)

## check abundances and load meta data
otuc<- abundances(clrc)
metac<- meta(clrc)

## as for Alpha diversity we create a new factor merging the samples from cloaca at day 30 of development and after day 30 together.
metac$Stage<-factor(metac$Age,levels=c("Embryo","Egg","0","5","10","15","20","25","30","42","44","51","61","62","EnvSand","EnvWater","ExtBlank","FieldBlank"),labels=c("UnhatchedEgg","HatchedEgg", "day0", "day5", "day10", "day15", "day20", "day25", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater","NestEnvironment","RiverEnvironment","Blank","Blank"))

#Permanovas
adonis(t(otuc)~Stage,data=metac, permutations=99999, method= "euclidean")

## We calculate euclidean distances for the Pairwise analysis
distc<- vegdist(t(otuc),method="euclidean")

#Pairwise table
pairwise.adonis(distc,metac$Stage)


#For the PCA
gp.dist <- dist(tableclrc, method="euclidean")
gp.pcoa <- ordinate(punifilis_control_comp, 'PCoA', distance=gp.dist)

plot_ordination(punifilis_control_comp, gp.pcoa, color="Stage1", shape="SampleForm2") + geom_point(size=2.5)+ stat_ellipse(type = "norm", linetype = 2)+geom_point(size=2.5)+theme_light()+scale_color_manual("Stage1", values = ana_co5.3, breaks=c("UnhatchedEgg","HatchedEgg","day0","day5","day10","day15","day20","day25","day30andGreater","NestEnvironment","RiverEnvironment","Blank"),labels=c("L.D.","Hatched Eggs","Day 0","Day 5","Day 10","Day 15","Day 20","Day 25","Day \u226530","Nest\nEnvironment","River\nEnvironment", "Blank")) + xlab("PCA 1 [26.4%]")+ ylab("PCA 2 [9.4%]")+theme(axis.text=element_text(size = 14))+theme(axis.title=element_text(size = 14))+ theme(legend.text = element_text(size = 11))+ scale_shape(name =" ", breaks=c("Egg","Egg","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Environment","Environment", "Blank"),labels=c("Egg","Egg","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Environment","Environment", "Blank"))+ theme(legend.position=c(1,1), legend.justification=c(1,1), legend.direction= "horizontal")



###Object without Blanks to use for Beta diversity dispersion###

# Taxa filtering
punifilis_pure_comp <-  phyloseq::filter_taxa(punifilis_pure, function(x) sum(x>0) > 30, TRUE)

#OTU table of the filtered data set
tableclr<-t(otu_table(punifilis_pure_comp))

#Our new table with CLR trasnformed data
tableclr <- clr_coordinates(tableclr+1)

#We extract data from our main data set
data2<-sample_data(punifilis_pure)
taxaclr<-tax_table(punifilis_pure)
treeclr<-phy_tree(punifilis_pure)

#We merge our data into a new phyloseq for our beta dispersion analysis
clr <- merge_phyloseq(tableclr,data2,taxaclr,treeclr)

otu<- abundances(clr)
meta<- meta(clr)

meta$Stage<-factor(meta$Age,levels=c("Embryo","Egg","0","5","10","15","20","25","30","42","44","51","61","62","EnvSand","EnvWater"),labels=c("UnhatchedEgg","HatchedEgg", "day0", "day5", "day10", "day15", "day20", "day25", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater","NestEnvironment","RiverEnvironment"))

# We calculate euclidean distances for beta diversity dispersion analysis
dist<- vegdist(t(otu),method="euclidean")
betadist<-betadisper(dist,meta$Stage, type = "median")$distance

meta$betadist<-betadist

## For the beta diversity dispersion plot we seaparate the environmental samples, since we are more interested in the biological result from the time of development of turtles
metaEnv1<-meta[which(!meta$Stage=="RiverEnvironment"),]
metaEnv<-meta[which(!metaEnv1$Stage=="NestEnvironment"),]
metaEnv

levels(metaEnv$Stage)

# We create a new factor with nicer labels for our plot
metaEnv$SampleFrom<-factor((metaEnv$Stage), levels = c("UnhatchedEgg","HatchedEgg", "day0", "day5", "day10", "day15", "day20", "day25", "day30andGreater"), labels= c("Inner eggshell","Inner eggshell","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca"))

levels(metaEnv$SampleFrom)

### Plot for beta diversity dispersion 
ggplot(metaEnv,aes(x = Stage, y = betadist,fill= SampleForm))+geom_boxplot(outlier.shape=NA)+theme_linedraw()+geom_jitter(colour= "black", position=position_dodge(0.8), size=1.5)+xlab(" ")+ylab("Beta Diversity Dispersion")+guides(fill=guide_legend(title = element_blank(), position = "right", reverse = TRUE, keywith = 0.5, keyheight = 0.5),shape=guide_legend(title = waiver(), position = "right", reverse = TRUE, keywith = 1, keyheight = 1))+scale_fill_manual(values= ana_co5,breaks=c("Egg","Turtle"),labels=c("Egg","Turtle"))+facet_grid(~SampleFrom, scales = "free", space = "free")+theme(strip.text.x = element_text(size = 14))+theme(axis.text=element_text(size = 12))+theme(axis.title=element_text(size = 14))+scale_x_discrete(labels = function(x) str_wrap(x, width = 18))+ theme(legend.text = element_text(size = 14))+ scale_x_discrete(labels=c("UnhatchedEgg" = "L.D.","HatchedEgg" = "Hatched","day0" = "Day 0", "day5" = "Day 5", "day10" = "Day 10", "day15" = "Day 15", "day20" = "Day 20", "day25" = "Day 25", "day30andGreater" = "Day \u226530"))



## GLS to control for microbial composition heterogeneity during turtles' development

Mdist2<- gls(betadist~Stage,na.action=na.fail,data=metaEnv, weights = varIdent(form = ~ 1|Stage),method="ML")

# Lets have a look to the results of our model
summary(Mdist2)

dredge(Mdist2, extra="adjR^2")

#Multiple comparison post-hoc beta diversity dispersion

posthocOptimalbeta<-summary(glht(Mdist2))
posthocOptimalbeta$test$sigma
posthocOptimalbeta$test

#read the table for plotting 
read.table(posthocOptimalBDF,"/home/student-gen-a/Dropbox/SharedFolder/MolEcolSecondSub/Tables/posthocOptimalBDF.txt")

posthocOptimalBDF$Sign<-ifelse(posthocOptimalBDF$Pvalue < 0.001, "***", ifelse(posthocOptimalBDF$Pvalue < 0.01, "**", ifelse(posthocOptimalBDF$Pvalue < 0.05, "*",""))) 

ggplot(posthocOptimalBDF, aes(reorder(StageFigbeta, Coefficient), Coefficient,label = Sign)) +geom_hline(yintercept=0, colour="#8C2318", size=1) +   geom_pointrange(aes(ymin=lwr, ymax=upr)) +  labs(x="Sample type comparison", y="Estimate", title="Beta diversity") + coord_flip() + theme_bw()+geom_text(nudge_y= 1.5, nudge_x=0.1, size = 5)+theme(axis.text=element_text(size = 10))+theme(axis.title=element_text(size = 14))+theme(title=element_text(size = 11))




############## Beta diversity Analysis withouth Filter for Supplementary Information ############## 


#### We create an object withouth applying the filtering step for betadiverstiy analysis. The data from this analysis is only to see if there is different in the results between filtered and non filtered data. The Graphics from this step are presented in the Supplementay Information of our manuscript.

############### Object with Blanks AND NO FILTER for PERMANOVA, Pairwise and PCoA graphics ##############

#### We do not apply any filter
punifilis_control_comp_nf <-  phyloseq::filter_taxa(punifilis_control, function(x) sum(x>0) > 0, TRUE)

tableclrc_nf<-t(otu_table(punifilis_control_comp_nf))

clr_coordinates <- function(X) {
  lX = log(X)
  lX - apply(lX, 1, mean)
} 


#Our new table with CLR trasnformed data
tableclrc_nf <- clr_coordinates(tableclrc_nf+1)

#We extract data from our main data set
data2c_nf<-sample_data(punifilis_control)
taxaclrc_nf<-tax_table(punifilis_control)
treeclrc_nf<-phy_tree(punifilis_control)

# we merge the data into a single phyloseq object
clrc_nf <- merge_phyloseq(tableclrc_nf,data2c_nf,taxaclrc_nf,treeclrc_nf)

otuc_nf<- abundances(clrc_nf)
metac_nf<- meta(clrc_nf)


metac_nf$Stage<-factor(metac_nf$Age,levels=c("Embryo","Egg","0","5","10","15","20","25","30","42","44","51","61","62","EnvSand","EnvWater","ExtBlank","FieldBlank"),labels=c("UnhatchedEgg","HatchedEgg", "day0", "day5", "day10", "day15", "day20", "day25", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater","NestEnvironment","RiverEnvironment","Blank","Blank"))

levels(metac_nf$Stage) 

permturtc_nf<-adonis(t(otuc_nf)~Stage,data=metac_nf, permutations=99999, method= "euclidean")

permturtc_nf

distc_nf<- vegdist(t(otuc_nf),method="euclidean")

pairwisec_nf<- pairwise.adonis(distc_nf,metac_nf$Stage)

pairwisec_nf

### PCoA withouth filtering
gp.dist <- dist(tableclrc_nf, method="euclidean")
gp.pcoa <- ordinate(punifilis_control_comp_nf, 'PCoA', distance=gp.dist)


plot_ordination(punifilis_control_comp_nf, gp.pcoa, color="Stage1", shape="SampleForm2") + geom_point(size=2.5)+ stat_ellipse(type = "norm", linetype = 2)+geom_point(size=2.5)+theme_light()+scale_color_manual("Stage1", values = ana_co5.3, breaks=c("UnhatchedEgg","HatchedEgg","day0","day5","day10","day15","day20","day25","day30andGreater","NestEnvironment","RiverEnvironment","Blank"),labels=c("L.D.","Hatched Eggs","Day 0","Day 5","Day 10","Day 15","Day 20","Day 25","Day \u226530","Nest\nEnvironment","River\nEnvironment", "Blank")) + xlab("PCA 1 [12.9%]")+ ylab("PCA 2 [5%]")+theme(axis.text=element_text(size = 14))+theme(axis.title=element_text(size = 14))+ theme(legend.text = element_text(size = 11))+ scale_shape(name =" ", breaks=c("Egg","Egg","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Environment","Environment", "Blank"),labels=c("Egg","Egg","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Turtle","Environment","Environment", "Blank"))+ theme(legend.position=c(0.28,0.92), legend.justification=c(1,1), legend.direction= "horizontal")


################# OBJECT WITHOUT BLANKS AND NO FILTER #########
################   FOR BETADIPS ANALYSIS  ##########
punifilis_pure_comp_nf <-  phyloseq::filter_taxa(punifilis_pure, function(x) sum(x>0) >0, TRUE)

tableclr_nf<-t(otu_table(punifilis_pure_comp_nf))

######## for data without blanks ########
tableclr_nf <- clr_coordinates(tableclr_nf+1)
data2_nf<-sample_data(punifilis_pure)
taxaclr_nf<-tax_table(punifilis_pure)
treeclr_nf<-phy_tree(punifilis_pure)

clr_nf <- merge_phyloseq(tableclr_nf,data2_nf,taxaclr_nf,treeclr_nf)

otu_nf<- abundances(clr_nf)
meta_nf<- meta(clr_nf)


meta_nf$Stage<-factor(meta_nf$Age,levels=c("Embryo","Egg","0","5","10","15","20","25","30","42","44","51","61","62","EnvSand","EnvWater"),labels=c("UnhatchedEgg","HatchedEgg", "day0", "day5", "day10", "day15", "day20", "day25", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater", "day30andGreater","NestEnvironment","RiverEnvironment"))

levels(meta_nf$Stage) 

dist_nf<- vegdist(t(otu_nf),method="euclidean")

betadist_nf<-betadisper(dist_nf,meta_nf$Stage, type = "median")$distance

meta_nf$betadist_nf<-betadist_nf

metaEnv1_nf<-meta_nf[which(!meta_nf$Stage=="RiverEnvironment"),]
metaEnv_nf<-meta_nf[which(!metaEnv1_nf$Stage=="NestEnvironment"),]
metaEnv_nf


metaEnv_nf$SampleFrom<-factor((metaEnv_nf$Stage), levels = c("UnhatchedEgg","HatchedEgg", "day0", "day5", "day10", "day15", "day20", "day25", "day30andGreater"), labels= c("Inner eggshell","Inner eggshell","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca","Hatchling cloaca"))

levels(metaEnv_nf$SampleFrom)

ggplot(metaEnv_nf,aes(x = Stage, y = betadist_nf,fill= SampleForm))+geom_boxplot(outlier.shape=NA)+theme_linedraw()+geom_jitter(colour= "black", position=position_dodge(0.8), size=1.5)+xlab(" ")+ylab("Beta Diversity Dispersion")+guides(fill=guide_legend(title = element_blank(), position = "right",reverse = TRUE, keywith = 0.5, keyheight = 0.5),shape=guide_legend(title = waiver(), position = "right", reverse = TRUE, keywith = 1, keyheight = 1))+scale_fill_manual(values= ana_co5.2,breaks=c("Env","Turtle","Egg"),labels=c("Environment","Turtle","Egg"))+facet_grid(~SampleFrom, scales = "free", space = "free")+theme(strip.text.x = element_text(size = 14))+theme(axis.text=element_text(size = 12))+theme(axis.title=element_text(size = 14))+scale_x_discrete(labels = function(x) str_wrap(x, width = 18))+ theme(legend.text = element_text(size = 14))+theme(legend.position=c(.6,.3), legend.direction="horizontal")+scale_x_discrete(labels=c("UnhatchedEgg" = "L.D.","HatchedEgg" = "Hatched","day0" = "Day 0", "day5" = "Day 5", "day10" = "Day 10", "day15" = "Day 15", "day20" = "Day 20", "day25" = "Day 25", "day30andGreater" = "Day \u226530"))


## Beta diversity dispersion model selection

Mdist2_nf<- gls(betadist_nf~Stage,na.action=na.fail,data=metaEnv_nf, weights = varIdent(form = ~ 1|Stage),method="ML")

summary(Mdist2_nf)

dredge(Mdist2_nf, extra="adjR^2")

#read table for plotting 
posthocOptimalBDF_nf<-read.table("/home/student-gen-a/Dropbox/SharedFolder/MolEcolSecondSub/Tables/posthocOptimalBDF_Sup.txt")

posthocOptimalBDF_nf$Sign<-ifelse(posthocOptimalBDF_nf$Pvalue < 0.001, "***", ifelse(posthocOptimalBDF_nf$Pvalue < 0.01, "**", ifelse(posthocOptimalBDF_nf$Pvalue < 0.05, "*",""))) 

ggplot(posthocOpt imalBDF_nf, aes(reorder(StageFigbeta, Coefficient), Coefficient,label = Sign)) +geom_hline(yintercept=0, colour="#8C2318", size=1) +   geom_pointrange(aes(ymin=lwr, ymax=upr)) +  labs(x="Sample type comparison", y="Estimate", title="Beta diversity") + coord_flip() + theme_bw()+geom_text(nudge_y= 2, nudge_x=0.1, size = 5)+theme(axis.text=element_text(size = 10))+theme(axis.title=element_text(size = 14))+theme(title=element_text(size = 11))





					#################################
					######### SOURCETRACKER #########
					#################################


## Sourcetracker data comparing the nest environment with egg, cloaca and river environment samples. Non normalized data.

# We load our metadata from the text file
meta<-read.table('MetaSourTracker.txt', sep='\t',h=T,row.names=1,check=F,comment='')

# We load our otu table
otus<-read.table('OTUTable.txt', sep='\t',h=T,row.names=1,check=F,comment='')

# We transform our otu table to a matrix table,
otus <- t(as.matrix(otut))

# We organize our tables by sample IDs
common.sample.ids <- intersect(rownames(meta), rownames(otus))
otus <- otus[common.sample.ids,]
meta <- meta[common.sample.ids,]

# We check if both of our tables, meta and otu table, are aligened and have the same names (by sample ID)
if(length(common.sample.ids) <= 1) {
    message <- paste(sprintf('Error: there are %d sample ids in common '),
                    'between the metadata file and data table')
    stop(message)
}

# We create a new factor where we give new labels to our levels. This is important since here we say which group of samples is the source and which are the sink.
meta$SourceSinkNest<-factor(meta$Stage1, levels=c("UnhatchedEgg","HatchedEgg","day0","day5","day10","day15","day20","day25", "day30andGreater","NestEnvironment","RiverEnvironment"), labels=c("sink","sink","sink","sink","sink","sink","sink","sink","sink","source","sink"))

#We create a new meta data four the next analysis
metaNestAll<-meta[which(meta$Stage1=="NestEnvironment" | meta$Stage1=="HatchedEgg" | meta$Stage1=="UnhatchedEgg" | meta$Stage1=="day0" | meta$Stage1=="day5" | meta$Stage1=="day10" | meta$Stage1=="day15" | meta$Stage1=="day20" | meta$Stage1=="day25" | meta$Stage1=="day30andGreater" | meta$Stage1=="RiverEnvironment"),]

# Now we do our Souerce Tracker analysis. Run the next steps only once.
# We repeat this for all the "Source" group that we want to compare.

#Let's call our source tracker function
source('SourceTracker.r')

# Here we train our data
train.ix <- which(metaNestAll$SourceSinkNest=='source')
test.ix <- which(metaNestAll$SourceSinkNest=='sink')
envs <- metaNestAll$SourceSinkNest
if(is.element('Description',colnames(metaNestAll))) desc <- metaNestAll$Description
alpha1 <- alpha2 <- 0.001

# Here we train our SourceTracker object on our previous trained data
str <- sourcetracker(otus[train.ix,], envs[train.ix], rarefaction_depth = 10000)

# Here we estimate source proportions in data
results_NestAll <- predict(st,otus[test.ix,], alpha1=alpha1, alpha2=alpha2, rarefaction_depth = 10000)

# We create a data frame with the calculated proportions
sourcetracker1_NestAll <-data.frame(results_NestAll$proportions)

# We save the our results into a txt. file, since these analyses take long time, you just have to repeat them once 
sourcetrackerNest<-write.table(sourcetracker1_NestAll,file = "SourcetrackerNest.txt", sep = "\t", row.names = TRUE, col.names = NA)

# We read the table and we are ready to plot and perform further analysis
metaNestAllSour<-read.table("SourceNestTable.txt", header=TRUE, sep="")


p6SouNest<-ggplot(metaNestAllSour, aes(x = Stage1, y = source,fill= SampleForm))+geom_boxplot(outlier.shape=NA)+theme_linedraw()+geom_jitter(colour= "black", position=position_dodge(0.8), size=1.5)+xlab(" ")+ylab("Source Nest")+ guides(fill=guide_legend(title = element_blank(), position = "right", reverse = FALSE, keywith = 0.5, keyheight = 0.5), shape=guide_legend(title = waiver(), position = "right", reverse = FALSE, keywith = 0.5, keyheight = 0.5))+ scale_fill_manual(values= ana_co5,breaks=c("Egg","Turtle","Env"),labels=c("Egg","Turtle","Environment"))+facet_grid(~SampleFrom, scales = "free", space = "free")+theme(strip.text.x = element_text(size = 14))+theme(axis.text=element_text(size = 12))+theme(axis.title=element_text(size = 14))+scale_x_discrete(labels = function(x) str_wrap(x, width = 18))+ theme(legend.text = element_text(size = 14))+ theme(legend.position=c(.59,.95), legend.direction="horizontal")+ scale_x_discrete(labels=c("UnhatchedEgg" = "L.D.","HatchedEgg" = "Hatched","day0" = "Day 0", "day5" = "Day 5", "day10" = "Day 10", "day15" = "Day 15", "day20" = "Day 20", "day25" = "Day 25", "day30andGreater" = "Day \u226530", "RiverEnvironment" = "River"))

p6SouNest


###### STATS

# we omit the nest environment label in our factor
metaNestAllSour$Stage1=="NestEnvironment"

M1<-glmmTMB(pos~Stage1,data=metaNestAllSour, family=nbinom2, na.action="na.fail")

# Lets have a look to the results of our model
summary(M1)
dredge(M1)

#call the system file for the next analysis 
source(system.file("other_methods","lsmeans_methods.R",package="glmmTMB"))

#For the multiple comparison calculations

ConfSource<-confint(glht(M1, linfct = mcp(Stage1 = "Tukey")))

### Here we have a table ready with all mul-comparison data. Load the table for plotting
posthocOptimalsourceDF<-read.table("posthocOptimalsourceDF.txt", header=TRUE, sep="")

## for significance of p-values
posthocOptimalsourceDF$Sign<-ifelse(posthocOptimalsourceDF$Pvalue < 0.001, "***", ifelse(posthocOptimalsourceDF$Pvalue < 0.01, "**", ifelse(posthocOptimalsourceDF$Pvalue < 0.05, "*",""))) 

ggplot(posthocOptimalsourcenestDF, aes(reorder(StageFigNest, Coefficient), Coefficient,label = Sign)) +geom_hline(yintercept=0, colour="#8C2318", size=1) +   geom_pointrange(aes(ymin=lwr, ymax=upr)) +  labs(x="Sample type comparison", y="Estimate", title="Source Nest") + coord_flip() + theme_bw()+geom_text(nudge_y=0.9, nudge_x=0.15, size = 5) +theme(axis.text=element_text(size = 10))+theme(axis.title=element_text(size = 14))+theme(title=element_text(size = 11))



###### with river samples #######
# we repeat the same analysis as above, we only change the source group, in this case the river. 

meta$SourceSinkRiver<-factor(meta$Stage1, levels=c("UnhatchedEgg","HatchedEgg","day0","day5","day10","day15","day20","day25", "day30andGreater","NestEnvironment","RiverEnvironment"), labels=c("sink","sink","sink","sink","sink","sink","sink","sink","sink","sink","source"))

metaRiverAll<-meta[which(meta$Stage1=="RiverEnvironment" | meta$Stage1=="HatchedEgg" | meta$Stage1=="UnhatchedEgg" | meta$Stage1=="day0" | meta$Stage1=="day5" | meta$Stage1=="day10" | meta$Stage1=="day15" | meta$Stage1=="day20" | meta$Stage1=="day25" | meta$Stage1=="day30andGreater" | meta$Stage1=="NestEnvironment"),]

#### call the source tracker function and perform the analysis as above.

## Load the table with source tracker analysis for ploting and further analysis
metaRiverAllSour<-read.table("SourceRiverTable.txt", header=TRUE, sep="")


ggplot(metaRiverAllSour, aes(x = Stage1, y = source,fill= SampleForm))+geom_boxplot(outlier.shape=NA)+theme_linedraw()+geom_jitter(colour= "black", position=position_dodge(0.8), size=1.5)+xlab(" ")+ylab("Source River")+facet_grid(~SampleFrom, scales = "free", space = "free")+ scale_fill_manual(values= ana_co5,breaks=c("Env","Turtle","Egg"),labels=c("Environment","Turtle","Egg")) + guides(fill=guide_legend(title = element_blank(), position = "right", reverse = TRUE, keywith = 0.5, keyheight = 0.5), shape=guide_legend(title = waiver(), position = "right", reverse = TRUE, keywith = 0.5, keyheight = 0.5))+theme(strip.text.x = element_text(size = 14))+theme(axis.text=element_text(size = 12))+theme(axis.title=element_text(size = 14))+scale_x_discrete(labels = function(x) str_wrap(x, width = 18))+ theme(legend.text = element_text(size = 14))+ theme(legend.position=c(.59,.95), legend.direction="horizontal")+ scale_x_discrete(labels=c("UnhatchedEgg" = "L.D.","HatchedEgg" = "Hatched","day0" = "Day 0", "day5" = "Day 5", "day10" = "Day 10", "day15" = "Day 15", "day20" = "Day 20", "day25" = "Day 25", "day30andGreater" = "Day \u226530", "NestEnvironment" = "Nest"))


### stats

M2<-glmmTMB(posR~Stage1,data=metaRiverAllSour, family=nbinom2, na.action="na.fail")

# Lets have a look to the results of our model
summary(M2)
dredge(M2)

##For the multiple comparison calculations
ConfSourceRiver<-confint(glht(M2, linfct = mcp(Stage1 = "Tukey")))

# Read table for plotting results of multiple comparison post-hoc
posthocOptimalsourceDFR<-read.table("posthocOptimalsourceDFR.txt", header=TRUE, sep="")

## for significance of p-values
posthocOptimalsourceDFR$Sign<-ifelse(posthocOptimalsourceDFR$Pvalue < 0.001, "***", ifelse(posthocOptimalsourceDFR$Pvalue < 0.01, "**", ifelse(posthocOptimalsourceDFR$Pvalue < 0.05, "*",""))) 

ggplot(posthocOptimalsourceriverDFR, aes(reorder(StageFigRiver, Coefficient), Coefficient,label = Sign)) +geom_hline(yintercept=0, colour="#8C2318", size=1) +   geom_pointrange(aes(ymin=lwr, ymax=upr)) +  labs(x="Sample type comparison", y="Estimate", title="Source River") + coord_flip() + theme_bw()+geom_text(nudge_y=1, nudge_x=0.15, size = 5) + theme(axis.text=element_text(size = 10))+theme(axis.title=element_text(size = 14))+theme(title=element_text(size = 11)) 



####### SOURCE TRACKER WITH HATCHED EGGS AS SOURCE without environment #########

metanoEnv<-meta[which(!meta$Stage1=="NestEnvironment"),]

metanoEnv<-metanoEnv[which(!metanoEnv$Stage1=="RiverEnvironment"),]

levels(metanoEnv$Stage1)

metanoEnv$SourceHatched<-factor(metanoEnv$Stage1, levels=c("UnhatchedEgg","HatchedEgg","day0","day5","day10","day15","day20","day25", "day30andGreater"), labels=c("sink","source","sink","sink","sink","sink","sink","sink","sink"))

metaHatchednoEnv<-metanoEnv[which(metanoEnv$Stage1=="UnhatchedEgg" | metanoEnv$Stage1=="HatchedEgg" | metanoEnv$Stage1=="day0" | metanoEnv$Stage1=="day5" | metanoEnv$Stage1=="day10" | metanoEnv$Stage1=="day15" | metanoEnv$Stage1=="day20" | metanoEnv$Stage1=="day25" | metanoEnv$Stage1=="day30andGreater"),]

## read the table with source tracker analysis data
metaHatchedSour<-read.table("SourceHatchedTable.txt", header=TRUE, sep="")

## the color used for this plot
ana_co5.4<- c("#bd8c10")

ggplot(metaHatchedSour, aes(x = Stage1, y = source,fill= SampleForm))+geom_boxplot(outlier.shape=NA)+theme_linedraw()+geom_jitter(colour= "black", position=position_dodge(0.8), size=1.5)+xlab(" ")+ylab("Source Hatched eggs")+facet_grid(~SampleFrom, scales = "free", space = "free")+ scale_fill_manual(values= ana_co5.4,breaks=c("Turtle"),labels=c("Turtle")) + guides(fill=guide_legend(title = element_blank(), position = "right", reverse = TRUE, keywith = 0.5, keyheight = 0.5), shape=guide_legend(title = waiver(), position = "right", reverse = TRUE, keywith = 0.5, keyheight = 0.5))+theme(strip.text.x = element_text(size = 14))+theme(axis.text=element_text(size = 12))+theme(axis.title=element_text(size = 14))+scale_x_discrete(labels = function(x) str_wrap(x, width = 18))+theme(legend.text = element_text(size = 14))+theme(legend.position="none", legend.direction="horizontal")+scale_x_discrete(labels=c("UnhatchedEgg" = "L.D.","HatchedEgg" = "Hatched","day0" = "Day 0", "day5" = "Day 5", "day10" = "Day 10", "day15" = "Day 15", "day20" = "Day 20", "day25" = "Day 25", "day30andGreater" = "Day \u226530"))


## Stats

M4<-glmmTMB(posR~Stage1,data=metaHatchedSour, family=nbinom2, na.action="na.fail")

# Lets have a look to the results of our model
summary(M4)
dredge(M4)

#### Read the table of multiple comparison post-hoc calculations 
posthocOptimalsourceHatchDFR<-read.table("posthocOptimalsourceHatchDFR.txt",header=TRUE, sep="")

## for significance of p-values
posthocOptimalsourceHatchDFR$Sign<-ifelse(posthocOptimalsourceHatchDFR$Pvalue < 0.001, "***", ifelse(posthocOptimalsourceHatchDFR$Pvalue < 0.01, "**", ifelse(posthocOptimalsourceHatchDFR$Pvalue < 0.05, "*",""))) 

ggplot(posthocOptimalsourceHatchDFR, aes(reorder(StageFigHatched, Coefficient), Coefficient,label = Sign)) +geom_hline(yintercept=0, colour="#8C2318", size=1) +   geom_pointrange(aes(ymin=lwr, ymax=upr)) +  labs(x="Sample type comparison", y="Estimate", title="Source Hatched") + coord_flip() + theme_bw()+geom_text(nudge_y=0.5, nudge_x=0.15, size = 5) + theme(axis.text=element_text(size = 10))+ theme(axis.title=element_text(size = 14))+theme(title=element_text(size = 11)) 


####### SOURCE TRACKER WITH EMBRYO EGGS AS SOURCE without environment. The graphics of this analysis are reported in the Supplementary infromation of the graphic #########

metaEmbryonoEnvSour<-read.table("SourceEmbryonoEnvTable.txt", header=TRUE, sep="")

# color pallet for this plot
ana_co5.5<- c("#bd8c10","#6164bd","#0e5200")

ggplot(metaEmbryonoEnvSour, aes(x = Stage1, y = source,fill= SampleForm))+geom_boxplot(outlier.shape=NA)+theme_linedraw()+geom_jitter(colour= "black", position=position_dodge(0.8), size=1.5)+xlab(" ")+ylab("Source Hatched eggs")+facet_grid(~SampleFrom, scales = "free", space = "free")+ scale_fill_manual(values= ana_co5.4,breaks=c("Turtle"),labels=c("Turtle")) + guides(fill=guide_legend(title = element_blank(), position = "right", reverse = TRUE, keywith = 0.5, keyheight = 0.5), shape=guide_legend(title = waiver(), position = "right", reverse = TRUE, keywith = 0.5, keyheight = 0.5))+theme(strip.text.x = element_text(size = 14))+theme(axis.text=element_text(size = 12))+theme(axis.title=element_text(size = 14))+scale_x_discrete(labels = function(x) str_wrap(x, width = 18))+theme(legend.text = element_text(size = 14))+theme(legend.position="none", legend.direction="horizontal")+scale_x_discrete(labels=c("UnhatchedEgg" = "L.D.","HatchedEgg" = "Hatched","day0" = "Day 0", "day5" = "Day 5", "day10" = "Day 10", "day15" = "Day 15", "day20" = "Day 20", "day25" = "Day 25", "day30andGreater" = "Day \u226530"))


### Stats 
M5<-glmmTMB(posR~Stage1,data=metaEmbryonoEnvSour, family=nbinom2, na.action="na.fail")

# Lets have a look to the results of our model
summary(M5)
dredge(M5)

#### Read the table for plotting results of multiple comparison post-hoc
posthocOptimalsourcEmbryoDFR<-read.table("posthocOptimalsourcEmbryoDFR.txt", header=TRUE, sep="")

## for significance of p-values
posthocOptimalsourcEmbryoDFR$Sign<-ifelse(posthocOptimalsourcEmbryoDFR$Pvalue < 0.001, "***", ifelse(posthocOptimalsourcEmbryoDFR$Pvalue < 0.01, "**", ifelse(posthocOptimalsourcEmbryoDFR$Pvalue < 0.05, "*",""))) 


ggplot(posthocOptimalsourcEmbryoDFR, aes(reorder(StageFigEmbryo, Coefficient), Coefficient,label = Sign)) +geom_hline(yintercept=0, colour="#8C2318", size=1) +   geom_pointrange(aes(ymin=lwr, ymax=upr)) +  labs(x="Sample type comparison", y="Estimate", title="Source Embryo") + coord_flip() + theme_bw()+geom_text(nudge_y=0.5, nudge_x=0.15) +theme(axis.text=element_text(size = 10)) +theme(axis.title=element_text(size = 14))+theme(title=element_text(size = 11))




					#################################
					############# ANCOM #############
					#################################


# Load the necessary packages
library(coin) #### Here I am using the package coin instead of "exactRankTests"
library("nlme")
library("compositions")

#Load the source code for ANCOM
source("ANCOM_updated_code.R")


### For this publication we compared the differences in ASV abundance between samples taken at different developmental times, starting with unhatched eggs (referred to as late embryonic eggs in the publication and in our dataset as unhatched eggs) up to day 30 or more than 30 days of development. We compared the different developmental stages where we observed a significant difference in ASV abundance.


#We start with L.D eggs ("UnhatchedEgg" in our data table") and hatched eggs

#we subset both sample groups and we merge them in a single object
EmbryoEggDF<-punifilis_pure %>% subset_samples(Stage1 == "UnhatchedEgg") 
HatchedEggDF<-punifilis_pure %>% subset_samples(Stage1 == "HatchedEgg") 
EmbryoEgg<-merge_phyloseq(EmbryoEggDF,HatchedEggDF)

#Filtering using the prevanlece threshhold as used for betadiversity
EmbryoEgg <-  phyloseq::filter_taxa(EmbryoEgg, function(x) sum(x>0) > length(sample_data(EmbryoEgg)$PCRControl)*0.3, TRUE)

#we create and object with an OTU table and their abundances, and we transform our object as data frame 
phylo_T <- t(abundances(EmbryoEgg))
phylo_B <- bind_cols(data.frame(Sample.ID = rownames(phylo_T)), data.frame(phylo_T))

#we create an object with the information of our table 
phylo_met <- meta(EmbryoEgg)

#we remove the X. from sample id names
colnames(phylo_met)[which(names(phylo_met) == "X.SampleID")] <- "Sample.ID"

#we check if our both meta files are consistent in the sample.id names
phylo_met$Sample.ID<-phylo_met$Featureid
phylo_met$Sample.ID==phylo_B$Sample.ID

### Here we perform ANCOM analysis
system.time({comparison_test_EmbryoEgg<- ANCOM.main(OTUdat=phylo_B, #OTU table
                           Vardat=phylo_met, #metadata file
                           adjusted=FALSE, # yes to include covariates
                           repeated=FALSE, # No repeated measurements
                           main.var="Stage1",
                           adj.formula=NULL,
                           longitudinal = FALSE, # No longitudinal data
                           random.formula=NULL, #follows the nlme syntax; isn't taken into consideration unless you have longitudinal data!!!
                           repeat.var=NULL, # No repeated measurements
                           multcorr=2, #1 is conservative, 2 is moderate and 3 is no correction
                           sig=0.05, #significance level
prev.cut=1) # OTUs with proportion of zeroes greater than prev.cut are not included in the analysis, so here we will not exclude any OTUs
})

## load the table where we have calculations from anova anlysis on this data set
fValue_dataframe_EmbryoEgg<-read.table("fValue_dataframe_EmbryoEgg.csv", header=TRUE, sep="")

#We merge our object containing the Ancom stats and we merge it together with our loaded table
volcano_df_EmbryoEgg <- merge(fValue_dataframe_EmbryoEgg, comparison_test_EmbryoEgg$W, by = "otu.names")
volcano_df_EmbryoEgg

### arrange the otu labels as characters
volcano_df_EmbryoEgg$otu.names <- as.character(volcano_df_EmbryoEgg$otu.names)

## we take the "X" from the otu names
volcano_df_EmbryoEgg$otu.names <- gsub('X','', volcano_df_EmbryoEgg$otu.names)

# we create an object with the names of the taxa in our data frame 
Tax<-data.frame(tax_table(EmbryoEgg))
Tax <- bind_cols(data.frame(ASV = rownames(Tax)), data.frame(Tax))
factor(Tax$ASV)
ASV <-ifelse(substring(Tax$ASV, 1,1) %in% c(0:9), paste0("X", Tax$ASV), paste(Tax$ASV))
Tax$ASV<-ASV
Tax$otu.names<-Tax$ASV
## we take the "X" from the otu names
Tax$otu.names <- gsub('X','', Tax$otu.names)

# we merge the two objects by otu.names
volcano_df_EmbryoEgg<-merge(volcano_df_EmbryoEgg,Tax,by="otu.names")

volcano_df_EmbryoEgg$St<-"EmbryoEgg"

volcano_df_EmbryoEgg

#we check our results and our data is ready to plot
volcano_df_EmbryoEgg$estimates

ggplot(volcano_df_EmbryoEgg) +
        geom_point(aes(x=estimates, y=W_stat, color= Order)) +
        scale_fill_gradient(limits=c(0, 50)) +
        ggtitle("Embryo - Hatched eggs") +
        xlab("clr estimates") + 
        ylab("W")



########### ANCOM cloaca samples day 0 - day 5 ###########


#we subset both sample groups and we merge them in a single object
day5DF<-punifilis_pure %>% subset_samples(Stage1 == "day5") 
day0_5<-merge_phyloseq(day5DF,day0DF)

#Filtering using the prevanlece threshhold as used for betadiversity
day0_5 <-  phyloseq::filter_taxa(day0_5, function(x) sum(x>0) > length(sample_data(day0_5)$PCRControl)*0.3, TRUE)

#we create and object with an OTU table and their abundances, and we transform our object as data frame 
phylo_T <- t(abundances(day0_5))#OTU table where row=sample; column =ASV 
phylo_B <- bind_cols(data.frame(Sample.ID = rownames(phylo_T)), data.frame(phylo_T)) #### Adding Sample.ID as a vector

phylo_met <- meta(day0_5)
str(phylo_met$Stage1)

colnames(phylo_met)[which(names(phylo_met) == "X.SampleID")] <- "Sample.ID"

phylo_met$Sample.ID<-phylo_met$Featureid
phylo_met$Sample.ID==phylo_B$Sample.ID


system.time({comparison_test_day0_5<- ANCOM.main(OTUdat=phylo_B, #OTU table
                           Vardat=phylo_met, #metadata file
                           adjusted=FALSE, # yes to include covariates
                           repeated=FALSE, # No repeated measurements
                           main.var="Stage1",
                           adj.formula=NULL,
                           longitudinal = FALSE, # No longitudinal data
                           random.formula=NULL, #follows the nlme syntax; isn't taken into consideration unless you have longitudinal data!!!
                           repeat.var=NULL, # No repeated measurements
                           multcorr=2, #1 is conservative, 2 is moderate and 3 is no correction
                           sig=0.05, #significance level
prev.cut=1) # OTUs with proportion of zeroes greater than prev.cut are not included in the analysis, so here we will not exclude any OTUs
})

#Table for Volcano plot
#here is the final table after merging our Ancom and our anova calculations. Read the table
volcano_df_day0_5<-read.table("Ancom_day0_5.txt", header=TRUE, sep="")

ggplot(volcano_df_day0_5) +
        geom_point(aes(x=estimates, y=W_stat,color= Order)) +
        scale_fill_gradient(limits=c(0, 50)) +
        ggtitle("Day 0 - Day 5") +
        xlab("clr estimates") + 
        ylab("W")


########### ANCOM cloaca samples day 5 - day 10 ###########


#we subset both sample groups and we merge them in a single object
day10DF<-punifilis_pure %>% subset_samples(Stage1 == "day10") 
day5_10<-merge_phyloseq(day5DF,day10DF)


#Filtering using the prevanlece threshhold as used for betadiversity
day5_10 <-  phyloseq::filter_taxa(day5_10, function(x) sum(x>0) > length(sample_data(day5_10)$PCRControl)*0.3, TRUE)


#we create and object with an OTU table and their abundances, and we transform our object as data frame 
phylo_T <- t(abundances(day5_10))#OTU table where row=sample; column =ASV 
phylo_B <- bind_cols(data.frame(Sample.ID = rownames(phylo_T)), data.frame(phylo_T)) #### Adding Sample.ID as a vector

phylo_met <- meta(day5_10)
str(phylo_met$Stage1)

colnames(phylo_met)[which(names(phylo_met) == "X.SampleID")] <- "Sample.ID"

phylo_met$Sample.ID<-phylo_met$Featureid
phylo_met$Sample.ID==phylo_B$Sample.ID

system.time({comparison_test_day5_10<- ANCOM.main(OTUdat=phylo_B, #OTU table
                           Vardat=phylo_met, #metadata file
                           adjusted=FALSE, # yes to include covariates
                           repeated=FALSE, # No repeated measurements
                           main.var="Stage1",
                           adj.formula=NULL,
                           longitudinal = FALSE, # No longitudinal data
                           random.formula=NULL, #follows the nlme syntax; isn't taken into consideration unless you have longitudinal data!!!
                           repeat.var=NULL, # No repeated measurements
                           multcorr=2, #1 is conservative, 2 is moderate and 3 is no correction
                           sig=0.05, #significance level
prev.cut=1) # OTUs with proportion of zeroes greater than prev.cut are not included in the analysis, so here we will not exclude any OTUs
})

#Table for Volcano plot
#here is the final table after merging our Ancom and our anova calculations. Read the table
volcano_df_day5_10<-read.table("Ancom_day5_10.txt", header=TRUE, sep="")

#check the data
volcano_df_day5_10

ggplot(volcano_df_day5_10) +
        geom_point(aes(x=estimates, y=W_stat,color=Class)) +
        scale_fill_gradient(limits=c(0, 50)) +
        ggtitle("Day 5 - Day 10") +
        xlab("clr estimates") + 
        ylab("W")



########### ANCOM cloaca samples day 25 - day 30 ###########

#we subset both sample groups and we merge them in a single object
day30DF<-punifilis_pure %>% subset_samples(Stage1 == "day30andGreater") 
day25_30<-merge_phyloseq(day25DF,day30DF)

#Filtering using the prevanlece threshhold as used for betadiversity
day25_30 <-  phyloseq::filter_taxa(day25_30, function(x) sum(x>0) > length(sample_data(day25_30)$PCRControl)*0.3, TRUE)


#we create and object with an OTU table and their abundances, and we transform our object as data frame 
phylo_T <- t(abundances(day25_30))#OTU table where row=sample; column =ASV 
phylo_B <- bind_cols(data.frame(Sample.ID = rownames(phylo_T)), data.frame(phylo_T)) #### Adding Sample.ID as a vector

phylo_met <- meta(day25_30)
colnames(phylo_met)[which(names(phylo_met) == "X.SampleID")] <- "Sample.ID"

phylo_met$Sample.ID<-phylo_met$Featureid
phylo_met$Sample.ID==phylo_B$Sample.ID


system.time({comparison_test_day25_30<- ANCOM.main(OTUdat=phylo_B, #OTU table
                           Vardat=phylo_met, #metadata file
                           adjusted=FALSE, # yes to include covariates
                           repeated=FALSE, # No repeated measurements
                           main.var="Stage1",
                           adj.formula=NULL,
                           longitudinal = FALSE, # No longitudinal data
                           random.formula=NULL, #follows the nlme syntax; isn't taken into consideration unless you have longitudinal data!!!
                           repeat.var=NULL, # No repeated measurements
                           multcorr=2, #1 is conservative, 2 is moderate and 3 is no correction
                           sig=0.05, #significance level
prev.cut=1) # OTUs with proportion of zeroes greater than prev.cut are not included in the analysis, so here we will not exclude any OTUs
})

#Table for Volcano plot
#here is the final table after merging our Ancom and our anova calculations. Read the table
volcano_df_day25_30<-read.table("Ancom_day25_30.txt", header=TRUE, sep="")

ggplot(volcano_df_day25_30) +
        geom_point(aes(x=estimates, y=W_stat,color=Class)) +
        scale_fill_gradient(limits=c(0, 50)) +
        ggtitle("Day 25- Day 30") +
        xlab("clr estimates") + 
        ylab("W")


########### ANCOM cloaca samples day 0 - day 30 ###########

#we subset both sample groups and we merge them in a single object
day30_0<-merge_phyloseq(day0DF,day30DF)

#Filtering using the prevanlece threshhold as used for betadiversity
day30_0 <-  phyloseq::filter_taxa(day30_0, function(x) sum(x>0) > length(sample_data(day30_0)$PCRControl)*0.3, TRUE)

#we create and object with an OTU table and their abundances, and we transform our object as data frame 
phylo_T <- t(abundances(day30_0))#OTU table where row=sample; column =ASV 
phylo_B <- bind_cols(data.frame(Sample.ID = rownames(phylo_T)), data.frame(phylo_T)) #### Adding Sample.ID as a vector

phylo_met <- meta(day30_0)
colnames(phylo_met)[which(names(phylo_met) == "X.SampleID")] <- "Sample.ID"

phylo_met$Sample.ID<-phylo_met$Featureid
phylo_met$Sample.ID==phylo_B$Sample.ID

system.time({comparison_test_day30_0<- ANCOM.main(OTUdat=phylo_B, #OTU table
                           Vardat=phylo_met, #metadata file
                           adjusted=FALSE, # yes to include covariates
                           repeated=FALSE, # No repeated measurements
                           main.var="Stage1",
                           adj.formula=NULL,
                           longitudinal = FALSE, # No longitudinal data
                           random.formula=NULL, #follows the nlme syntax; isn't taken into consideration unless you have longitudinal data!!!
                           repeat.var=NULL, # No repeated measurements
                           multcorr=2, #1 is conservative, 2 is moderate and 3 is no correction
                           sig=0.05, #significance level
prev.cut=1) # OTUs with proportion of zeroes greater than prev.cut are not included in the analysis, so here we will not exclude any OTUs
})

#Table for Volcano plot
#here is the final table after merging our Ancom and our anova calculations. Read the table
volcano_df_day30_0<-read.table("Ancom_day30_0.txt", header=TRUE, sep="")

ggplot(volcano_df_day30_0) +
        geom_point(aes(x=estimates, y=W_stat,color=Class)) +
        scale_fill_gradient(limits=c(0, 50)) +
        ggtitle("Day 25- Day 30") +
        xlab("clr estimates") + 
        ylab("W")



########################################################################################################################################################################################################################################################################

########## ANCOM SUP.INFORMATION

#We repeat the same for the next sample groups. Following we merge the samples groups that we want to compare and we repeat the same analysis as above
#############ANALYSIS EGGS DAY 0 ###############

HatchedEggDF<-punifilis_pure %>% subset_samples(Stage1 == "HatchedEgg") 
day0DF<-punifilis_pure %>% subset_samples(Stage1 == "day0") 
egg0<-merge_phyloseq(HatchedEggDF,day0DF)
egg0 <-  phyloseq::filter_taxa(egg0, function(x) sum(x>0) > length(sample_data(egg0)$PCRControl)*0.3, TRUE)

phylo_T <- t(abundances(egg0))#OTU table where row=sample; column =ASV 
phylo_B <- bind_cols(data.frame(Sample.ID = rownames(phylo_T)), data.frame(phylo_T)) #### Adding Sample.ID as a vector

phylo_met <- meta(egg0)
colnames(phylo_met)[which(names(phylo_met) == "X.SampleID")] <- "Sample.ID"

phylo_met$Sample.ID<-phylo_met$Featureid
phylo_met$Sample.ID==phylo_B$Sample.ID


system.time({comparison_test_egg0<- ANCOM.main(OTUdat=phylo_B, #OTU table
                           Vardat=phylo_met, #metadata file
                           adjusted=FALSE, # yes to include covariates
                           repeated=FALSE, # No repeated measurements
                           main.var="Stage1",
                           adj.formula=NULL,
                           longitudinal = FALSE, # No longitudinal data
                           random.formula=NULL, #follows the nlme syntax; isn't taken into consideration unless you have longitudinal data!!!
                           repeat.var=NULL, # No repeated measurements
                           multcorr=2, #1 is conservative, 2 is moderate and 3 is no correction
                           sig=0.05, #significance level
prev.cut=1) # OTUs with proportion of zeroes greater than prev.cut are not included in the analysis, so here we will not exclude any OTUs
})

## read table
volcano_df_egg0<-read.table("Ancom_egg_day0.txt", header=TRUE, sep="")

Wclr_egg0<-ggplot(volcano_df_egg0) +
        geom_point(aes(x=estimates, y=W_stat, color= Phylum)) +
        scale_fill_gradient(limits=c(0, 50)) +
        ggtitle("Eggs - Day 0") +
        xlab("clr estimates") + 
        ylab("W")



########### ANCOM cloaca samples day 10 - day 15 ###########

#day10DF
day15DF<-punifilis_pure %>% subset_samples(Stage1 == "day15") 
day10_15<-merge_phyloseq(day10DF,day15DF)
day10_15 <-  phyloseq::filter_taxa(day10_15, function(x) sum(x>0) > length(sample_data(day10_15)$PCRControl)*0.3, TRUE)


phylo_T <- t(abundances(day10_15))#OTU table where row=sample; column =ASV 
phylo_B <- bind_cols(data.frame(Sample.ID = rownames(phylo_T)), data.frame(phylo_T)) #### Adding Sample.ID as a vector


phylo_met <- meta(day10_15)
colnames(phylo_met)[which(names(phylo_met) == "X.SampleID")] <- "Sample.ID"

phylo_met$Sample.ID<-phylo_met$Featureid
phylo_met$Sample.ID==phylo_B$Sample.ID


system.time({comparison_test_day10_15<- ANCOM.main(OTUdat=phylo_B, #OTU table
                           Vardat=phylo_met, #metadata file
                           adjusted=FALSE, # yes to include covariates
                           repeated=FALSE, # No repeated measurements
                           main.var="Stage1",
                           adj.formula=NULL,
                           longitudinal = FALSE, # No longitudinal data
                           random.formula=NULL, #follows the nlme syntax; isn't taken into consideration unless you have longitudinal data!!!
                           repeat.var=NULL, # No repeated measurements
                           multcorr=2, #1 is conservative, 2 is moderate and 3 is no correction
                           sig=0.05, #significance level
prev.cut=1) # OTUs with proportion of zeroes greater than prev.cut are not included in the analysis, so here we will not exclude any OTUs
})

## read table
volcano_df_day10_15<-read.table("Ancom_day10_15.txt", header=TRUE, sep="")


Wclr_day10_15<-ggplot(volcano_df_day10_15) +
        geom_point(aes(x=estimates, y=W_stat,color=Class)) +
        scale_fill_gradient(limits=c(0, 50)) +
        ggtitle("Day 10- Day 15") +
        xlab("clr estimates") + 
        ylab("W")


########### ANCOM cloaca samples day 15 - day 20 ###########

day20DF<-punifilis_pure %>% subset_samples(Stage1 == "day20") 
day15_20<-merge_phyloseq(day15DF,day20DF)

day15_20 <-  phyloseq::filter_taxa(day15_20, function(x) sum(x>0) > length(sample_data(day15_20)$PCRControl)*0.3, TRUE)


phylo_T <- t(abundances(day15_20))#OTU table where row=sample; column =ASV 
phylo_B <- bind_cols(data.frame(Sample.ID = rownames(phylo_T)), data.frame(phylo_T)) #### Adding Sample.ID as a vector

phylo_met <- meta(day15_20)
colnames(phylo_met)[which(names(phylo_met) == "X.SampleID")] <- "Sample.ID"

phylo_met$Sample.ID<-phylo_met$Featureid
phylo_met$Sample.ID==phylo_B$Sample.ID


system.time({comparison_test_day15_20<- ANCOM.main(OTUdat=phylo_B, #OTU table
                           Vardat=phylo_met, #metadata file
                           adjusted=FALSE, # yes to include covariates
                           repeated=FALSE, # No repeated measurements
                           main.var="Stage1",
                           adj.formula=NULL,
                           longitudinal = FALSE, # No longitudinal data
                           random.formula=NULL, #follows the nlme syntax; isn't taken into consideration unless you have longitudinal data!!!
                           repeat.var=NULL, # No repeated measurements
                           multcorr=2, #1 is conservative, 2 is moderate and 3 is no correction
                           sig=0.05, #significance level
prev.cut=1) # OTUs with proportion of zeroes greater than prev.cut are not included in the analysis, so here we will not exclude any OTUs
})

#read table
volcano_df_day15_20<-read.table("Ancom_day15_20.txt", header=TRUE, sep="")

ggplot(volcano_df_day15_20) +
        geom_point(aes(x=estimates, y=W_stat,color=Class)) +
        scale_fill_gradient(limits=c(0, 50)) +
        ggtitle("Day 15- Day 20") +
        xlab("clr estimates") + 
        ylab("W")





########### ANCOM cloaca samples day 20 - day 25 ###########

day25DF<-punifilis_pure %>% subset_samples(Stage1 == "day25") 
day20_25<-merge_phyloseq(day20DF,day25DF)
day20_25 <-  phyloseq::filter_taxa(day20_25, function(x) sum(x>0) > length(sample_data(day20_25)$PCRControl)*0.3, TRUE)


phylo_T <- t(abundances(day20_25))#OTU table where row=sample; column =ASV 
phylo_B <- bind_cols(data.frame(Sample.ID = rownames(phylo_T)), data.frame(phylo_T)) #### Adding Sample.ID as a vector

phylo_met <- meta(day20_25)
colnames(phylo_met)[which(names(phylo_met) == "X.SampleID")] <- "Sample.ID"

phylo_met$Sample.ID<-phylo_met$Featureid
phylo_met$Sample.ID==phylo_B$Sample.ID


system.time({comparison_test_day20_25<- ANCOM.main(OTUdat=phylo_B, #OTU table
                           Vardat=phylo_met, #metadata file
                           adjusted=FALSE, # yes to include covariates
                           repeated=FALSE, # No repeated measurements
                           main.var="Stage1",
                           adj.formula=NULL,
                           longitudinal = FALSE, # No longitudinal data
                           random.formula=NULL, #follows the nlme syntax; isn't taken into consideration unless you have longitudinal data!!!
                           repeat.var=NULL, # No repeated measurements
                           multcorr=2, #1 is conservative, 2 is moderate and 3 is no correction
                           sig=0.05, #significance level
prev.cut=1) # OTUs with proportion of zeroes greater than prev.cut are not included in the analysis, so here we will not exclude any OTUs
})

## read table
volcano_df_day20_25<-read.table("Ancom_day20_25.txt", header=TRUE, sep="")

ggplot(volcano_df_day20_25) +
        geom_point(aes(x=estimates, y=W_stat,color=Class)) +
        scale_fill_gradient(limits=c(0, 50)) +
        ggtitle("Day 20- Day 25") +
        xlab("clr estimates") + 
        ylab("W")


