
#This code is organized as follows:
#section 1: Plant species composition
#section 2: H. incana growth and phenology
#section 3: H. incana fitness

##################################
########### section 1 ############
##################################

######## Packages & data #########

library(openxlsx)
library(ade4)
library(vegan)
library(FD)
library(devtools)
library(pairwiseAdonis)
library(Rtsne)
library(ggplot2)
library(ggdendro)
library(dendextend)
library(indicspecies)
library(plyr)
library(dplyr)
library(lares)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)
library(emmeans)
library(glmmTMB)
library(fitdistrplus)
library(DHARMa)
library(readxl)

#loading data-set
vegetation <- read.xlsx("percent_cover_summeddata_only.xlsx", colNames = TRUE, rowNames = FALSE, detectDates = TRUE)


####### data clean-up #########

#make distance an ordered factor
vegetation$distance <- factor(vegetation$distance, ordered = TRUE, 
                              levels = c("0", "2.5", "5", "10"))

#replace NA w/ 0
vegetation[is.na(vegetation)] <- 0

#remove site info columns
vegetation2 = vegetation[,6:25]


##### overall species info ######

#number of samples where species are present
spe.pres2 <- apply(vegetation2 > 0, 2, sum)
sort(spe.pres2)

#frequency of species present
spe.relf2 <- 100*spe.pres2/nrow(vegetation2)
round(sort(spe.relf2), 1)


###### diversity indices ######

SR <- rowSums(vegetation2 > 0) ## Species richness
H <- diversity(vegetation2) ## Shannon entropy
Div1 <- exp(H) ## Shannon's diversity (number of abundant species)
Div2 <- diversity(vegetation2, "inv") ## Simpson diversity (number of dominant species)
Eve <- H/log(SR) ## Pielou evenness
E10 <- Div1/SR ## Shannon evenness (Hill's ratio)
E20 <- Div2/SR ## Simpson evenness (Hill's ratio)
div <- data.frame(SR, H, Div1, Div2, E10, E20, Eve)## create a dataframe

#Comparing SR and evenness among the distances

#model decision making
hist(Eve) #histogram
ggdensity(Eve, 
          main = "Density plot of evenness",
          xlab = "evenness") #density plot
ggqqplot(Eve) #Q-Q plot
shapiro.test(Eve) #non-normal
library(goft)
gamma_test(Eve) #not gamma

#adding Eve to data frame
vegetation4 <- cbind.data.frame(vegetation, Eve = Eve)
#remove rows with NA
vegetation5 <- na.omit(vegetation4)

descdist(vegetation5$Eve) #looks like beta


Eve.aov = glmmTMB(Eve ~ distance + (1|nest), 
                  family = ordbeta(link = "logit"), data = vegetation5) #model
summary(Eve.aov)
Anova(Eve.aov)


simulationOutput <- simulateResiduals(fittedModel = Eve.aov) #model fit
plot(simulationOutput)

#emmeans: post-hoc pairwise comparison
comp <- emmeans(Eve.aov, pairwise ~ distance)
pairs(comp)


#model descision making
hist(SR) #histogram
ggdensity(SR, 
          main = "Density plot of SR",
          xlab = "evenness") #density plot
ggqqplot(SR) #Q-Q plot
shapiro.test(SR) #non-normal


#adding SR to data frame
vegetation5 <- cbind.data.frame(vegetation, SR = SR)
#remove rows with NA
vegetation6 <- na.omit(vegetation5)
descdist(vegetation6$SR) #looks like normal


SR.aov = glmmTMB(SR ~ distance + (1|nest), data = vegetation6)
SR2.aov = glmmTMB(SR ~ distance + (1|nest), family = poisson(link = "log"), data = vegetation6)
#even though the data is technically poisson, a  normal distribution looks better in the residuals
#and the AIC is lower for the normal

simulationOutput <- simulateResiduals(fittedModel = SR.aov)
plot(simulationOutput)

summary(SR.aov)
Anova(SR.aov)

emmeans(SR.aov, pairwise ~ distance)


###### beta diversity ######

## Bray-Curtis dissimilarity matrix on relative abundances
spe.rel <- decostand(vegetation2, "total")
head(spe.rel)
spe.dbrel <- vegdist(spe.rel, 'bray')
str(spe.dbrel)

VegBeta = betadisper(spe.dbrel, vegetation$distance)
VegBeta

## we can quickly test if there are different
anova(VegBeta)
TukeyHSD(VegBeta) #post-hoc


## now lets visualize our patterns
plot(VegBeta)
plot(VegBeta, main = "Beta Diversity Plot", cex = 0.5)

#boxplot
boxplot(VegBeta, xlab="Sample location")


####### PERMANOVA #######

adonis2(spe.dbrel ~ vegetation$distance + (1|vegetation$nest), permutations = 999)

library(pairwiseAdonis)
pairwise.adonis(spe.dbrel,vegetation$distance, 
                p.adjust.m = "fdr", perm = 999)


####### Species Indicator Analysis #######

library(indicspecies)

mig.ind.sp = multipatt(vegetation2, vegetation$distance,
                       func = "r.g", control = how(nperm=999))

summary(mig.ind.sp)


###### NMDS ######

#show species as arrows

colnames(vegetation2) <- gsub("\\.", " ", colnames(vegetation2))

SITEinfo <- vegetation[c(1:5)]

ant.nmds=metaMDS(vegetation2, distance="bray", try=100)
stressplot(ant.nmds)
plot(ant.nmds)

data.scores <- scores(ant.nmds)  
data.scores$distance <- as.data.frame(data.scores$distance)

data.scores$site <- rownames(data.scores)  
data.scores$grp<-SITEinfo$distance
species.scores <- as.data.frame(scores(ant.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores
coord_equal()
SN.envfit <- envfit(ant.nmds, SITEinfo, permutations = 999) # this fits environmental vectors
SN.spp.fit <- envfit(ant.nmds, vegetation2, permutations = 999) # this fits species vectors
#if you run 'SN.spp.fit' you can see which species are significant

site.scrs <- as.data.frame(scores(ant.nmds, display = "sites")) #save NMDS results into dataframe
site.scrs <- cbind(site.scrs, Distance = SITEinfo$distance) #add grouping variable "Management" to dataframe
spp.scrs <- as.data.frame(scores(SN.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
spp.scrs <- cbind(spp.scrs, pval = SN.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show only species significant at 0.05, can change the value if you want to show generally all spp
env.scores.SN <- as.data.frame(scores(SN.envfit, display = "vectors"))


library(Rmisc)


all_spp<-ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  stat_ellipse(aes(x=NMDS1,y=NMDS2,color=as.factor(SITEinfo$distance)),level = 0.95) +
  geom_point(aes(NMDS1, NMDS2,  colour = as.factor(SITEinfo$distance)), size = 4)+ #adds site points to plot, shape determined by Landuse, colour determined by Management
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))+ # add legend at right of plot
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.2, "cm")), colour = "blue", lwd=1) + 
  ggrepel::geom_text_repel(data = sig.spp.scrs,
                           aes(label = Species), cex = 5, direction = "both", segment.size = 0.25,
                           nudge_x=0.29,nudge_y=0.05)
all_spp+
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.8, "cm")), colour = "grey10", lwd=0.8) + #add vector arrows of significant species
  ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = SITEinfo$distance), cex = 5, direction = "both", segment.size = 0.6) #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap


all_spp2<-ggplot(site.scrs, aes(x=NMDS1, y=NMDS2))+ #sets up the plot
  stat_ellipse(aes(x=NMDS1,y=NMDS2,color=as.factor(SITEinfo$distance)),level = 0.95) +
  geom_point(aes(NMDS1, NMDS2, colour = as.factor(SITEinfo$distance))) +
  coord_fixed()+
  theme_classic()+ 
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "top", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10))+ # add legend at right of plot
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.2, "cm")), colour = "blue", lwd=1) + 
  ggrepel::geom_text_repel(data = sig.spp.scrs,
                           aes(label = paste0("italic(", gsub(" ", "~", Species), ")")), 
                           cex = 5, direction = "both", segment.size = 0.25,
                           nudge_x=0.29,nudge_y=0.05)+
  theme(axis.text = element_text(size = 18))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 18))


#going to summarize the points to see a general trend
nm1<-summarySE(site.scrs, measurevar="NMDS1", groupvars = "Distance")
nm2<-summarySE(site.scrs, measurevar="NMDS2", groupvars = "Distance")
nm.graph<-data.frame(nm1, nm2)

sum <- ggplot(nm.graph, aes(x=NMDS1, y=NMDS2))+
  geom_point(size = 10, aes(color=as.factor(Distance), shape=as.factor(Distance)))+
  geom_errorbar(aes(ymin=NMDS2-se, ymax=NMDS2+se, color=as.factor(Distance)), width=0.07,size=1)+
  geom_errorbarh(aes(xmin=NMDS1-se, xmax=NMDS1+se, color=as.factor(Distance)), height=0.07,size=1)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.2, "cm")), colour = "blue", lwd=1) + 
  ggrepel::geom_text_repel(data = sig.spp.scrs,
                           aes(label = paste0("italic(", gsub(" ", "~", Species), ")")), 
                           cex = 5, direction = "both", segment.size = 0.25,
                           nudge_x=0.29,nudge_y=0.05)+ 
  scale_shape_manual(values=c(15,17,20,25))+
  coord_fixed()+
  theme_classic()+
  theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
  theme(legend.position = "top")+
  theme(axis.text = element_text(size = 18))+
  theme(axis.title = element_text(size = 20))+
  theme(legend.text = element_text(size = 18))


sum <- ggplot(nm.graph, aes(x = NMDS1, y = NMDS2)) +
  geom_point(size = 10, aes(color = as.factor(Distance), shape = as.factor(Distance))) +
  geom_errorbar(aes(ymin = NMDS2 - se, ymax = NMDS2 + se, color = as.factor(Distance)), width = 0.07, size = 1) +
  geom_errorbarh(aes(xmin = NMDS1 - se, xmax = NMDS1 + se, color = as.factor(Distance)), height = 0.07, size = 1) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_segment(
    data = sig.spp.scrs,
    aes(x = 0, xend = NMDS1, y = 0, yend = NMDS2),
    arrow = arrow(length = unit(0.2, "cm")),
    colour = "blue", lwd = 1
  ) + 
  ggrepel::geom_text_repel(
    data = sig.spp.scrs,
    aes(label = paste0("italic(", Species, ")")),
    cex = 5, direction = "both", segment.size = 0.25,
    nudge_x = 0.20, nudge_y = 0.17,
    parse = TRUE
  ) +
  scale_shape_manual(values = c(15, 17, 20, 18)) +
  coord_fixed() +
  theme_classic() +
  labs(
    color = "Distance (m)", # Set the title for the color legend
    shape = "Distance (m)"  # Set the title for the shape legend
  ) +
  theme(
    panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"),
    legend.position = "right", # Move legend to the right
    legend.direction = "vertical", # Display legend items vertically
    legend.key = element_rect(fill = NA, color = NA), # Remove the black box around shapes
    axis.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )


#no convergence is a warning
NMDS.1= metaMDS(vegetation2, distance = "bray")
NMDS.1.1= NMDS.1$points[,1] #making NMDS into dataframe
NMDS.1.2 = NMDS.1$points[,2]

NMDS.new = data.frame(MDS1 = NMDS.1.1, MDS2 = NMDS.1.2)

#extract NMDS scores (x and y coordinates)
NMDS.new2 = as.data.frame(scores(NMDS.new))
NMDS.new2$distance<- vegetation$distance

NMDS.new2$nest<- as.character(vegetation$nest)

ggplot(NMDS.new2, aes(x=NMDS.1.1, y=NMDS.1.2, col=distance)) +
  geom_point() + xlab("NMDS1") + ylab("NMDS2") + stat_ellipse(lwd = 1.2) 
#Is each dot a different sample or plant?

ggplot(NMDS.new2, aes(x=NMDS.1.1, y=NMDS.1.2, col=nest)) +
  geom_point() + xlab("NMDS1") + ylab("NMDS2") + stat_ellipse(lwd = 1.2) 

species.scores <- as.data.frame(scores(NMDS.1, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores)  # create a column of species, from the rownames of species.scores

SN.spp.fit <- envfit(ant.nmds, vegetation2, permutations = 999) # this fits species vectors
#if you run 'SN.spp.fit' you can see which species are significant


######## Percent Cover of H. incana GLMM ########

hist(vegetation$Hirschfeldia.incana)
ggdensity(vegetation$Hirschfeldia.incana, 
          main = "Density plot of h.incana",
          xlab = "h.incana")
ggqqplot(vegetation$Hirschfeldia.incana)
shapiro.test(vegetation$Hirschfeldia.incana) #not normal
library(fitdistrplus)
descdist(vegetation$Hirschfeldia.incana) #looks like beta

vegetation$logincana = log(vegetation$Hirschfeldia.incana)
hist(vegetation$logincana)

#testing whether total cover of all other species is correlated with H. incana
# Add a new column 'total_cover' which is the sum of selected columns
vegetation$total_cover <- rowSums(vegetation[, c("Pectocarya.linearis", "Lasthenia.californica", 
                                                 "Bromus.rubens", "Erodium.cicutarium", 
                                                 "Oncosiphon.piluliferum", "Amsinckia.menziesii", 
                                                 "Salvia.columbariae", "Plagiobothrys.collinus", 
                                                 "Encelia.farinosa", "Eriogonum.fasciculatum", 
                                                 "Cryptantha.muricata", "Camissoniopsis.bistorta", 
                                                 "Opuntia.parryi", "Mimulus.aurantiacus", 
                                                 "Corethrogyne.filaginifolia", "Tropidocarpum.gracile", 
                                                 "Stephanomeria.spp.", "Eriogonum.gracile", 
                                                 "Emmenanthe.penduliflora")], na.rm = TRUE)

model_tot <- glmmTMB(Hirschfeldia.incana ~ distance + total_cover + (1|nest),
                     family = ordbeta(link = "logit"),
                     data = vegetation)

simulationOutput <- simulateResiduals(fittedModel = model_tot)
plot(simulationOutput)

summary(model_tot)

Anova(model_tot)

emmeans(model_tot, pairwise ~ distance)
#sig negative relationship with total cover


##################################
########### section 2 ############
##################################

######## Packages ########

library(openxlsx)
library(plyr)
library(dplyr)
library(lares)
library(ggplot2)
library(rstatix)
library(ggpubr)
library(tidyverse)
library(data.table)
library(lme4)
library(lmerTest)
library(emmeans)
library(brms)
library(jtools)
library(devtools)
library("glmmTMB")
library(ggeffects)
library(DHARMa)
library(goft)
library(fitdistrplus)
library(AICcmodavg)
library(DescTools)
library(car)
library(interactions)
library(jtools)

se <- function(x) (sd(x)/sqrt(length(x))) # add the se function
percent <- function(x) (mean(x)*100)
percent_without_NA <- function(x) (mean(x, na.rm = TRUE) * 100)

####### load data #######

H.incana.2022.raw <- read.xlsx("H.incana_phenology_2022_only.xlsx", colNames = TRUE, rowNames = FALSE, detectDates = TRUE)

####### data clean-up #######

#changing column names
colnames(H.incana.2022.raw)[c(6:7)] <- c("Distance.m","Height.cm") # fix the weird names
colnames(H.incana.2022.raw)[c(13)] <- c("rainfall")

#converting factos
H.incana.2022.raw$Nest <- as.character(H.incana.2022.raw$Nest)

#add H. incana abundance to the data set
options(dplyr.summarise.inform = FALSE)

H.incana.abun2022 <- H.incana.2022.raw %>%
  group_by(Nest, Date, Degree, Distance.m) %>%
  summarise(abundance=length(Height.cm))

H.incana.2022 <- merge(H.incana.2022.raw, H.incana.abun2022, by = c("Nest", "Date", "Degree", "Distance.m"))

#Convert relevant yes/no to 1/0
H.incana.flowering2 <- H.incana.2022$Flowering<-ifelse(H.incana.2022$Flowering=="yes",1,0)
H.incana.senesced2 <- H.incana.2022$Premature.senensce<-ifelse(H.incana.2022$Premature.senensce=="yes",1,0)


#duplicate distance and make it categorical
H.incana.2022$Distance.n = H.incana.2022$Distance.m
H.incana.2022$Distance.n <- factor(H.incana.2022$Distance.n, ordered = TRUE, 
                                   levels = c("0", "2.5", "5"))

###### height 2022 model descision making

#looking at our data distribution
hist(H.incana.2022.raw$Height.cm) #this does not look normal to me
#density plot
ggdensity(H.incana.2022.raw$Height.cm, 
          main = "Density plot of height",
          xlab = "height") #does not look bell shaped
#Q-Q plot
ggqqplot(H.incana.2022.raw$Height.cm) #looks logarithmic?

#shapiro test
shapiro.test(H.incana.2022.raw$Height.cm)
#indeed, the dist is non-normal

gamma_fit(H.incana.2022.raw$Height.cm)
gamma_test(H.incana.2022.raw$Height.cm)

####### height 2022 model #######

P.5 <- glmmTMB(Height.cm ~ Distance.n + Days.Elapsed + abundance + (1|Nest) +(1|rainfall), family=Gamma(link="log"), data = H.incana.2022)
simulationOutput <- simulateResiduals(fittedModel = P.5)
plot(simulationOutput)
#the problem lies in the outliers, there are ~80

#I will "winsorize" the response
H.incana.2022$wHeight.cm <- Winsorize(H.incana.2022$Height.cm, probs = c(0.05, 0.95))

hist(H.incana.2022$wHeight.cm)
descdist(H.incana.2022$wHeight.cm) #close to gamma but not quite it

#comparing AICs among models
P.1 <- glmmTMB(wHeight.cm ~ Days.Elapsed + (1|Nest) + (1|rainfall), family=Gamma(link="log"), data = H.incana.2022)
#AIC = 5001.146, BIC = 5027.768
P.2 <- glmmTMB(wHeight.cm ~ Days.Elapsed + abundance + (1|Nest) +(1|rainfall), family=Gamma(link="log"), data = H.incana.2022)
#AIC = 4978.91, BIC = 5010.856
P.3 <- glmmTMB(wHeight.cm ~ Distance.n + Days.Elapsed + (1|Nest) +(1|rainfall), family=Gamma(link="log"), data = H.incana.2022)
#AIC = 4941.622, BIC = 4978.893
P.4 <- glmmTMB(wHeight.cm ~ Distance.n + Days.Elapsed + abundance + (1|Nest) +(1|rainfall), family=Gamma(link="log"), data = H.incana.2022)
#AIC = 4922.241, BIC = 4964.837 <- selected model
P.5 <- glmmTMB(wHeight.cm ~ Distance.n + Days.Elapsed*abundance + (1|Nest) +(1|rainfall), family=Gamma(link="log"), data = H.incana.2022)
#AIC = 4922.095, BIC = 4970.015 
P.6 <- glmmTMB(wHeight.cm ~ Days.Elapsed*abundance + (1|Nest) +(1|rainfall), family=Gamma(link="log"), data = H.incana.2022)
#AIC = 4979.818, BIC = 5017.09
P.7 <- glmmTMB(wHeight.cm ~ Distance.n*Days.Elapsed + abundance + (1|Nest) +(1|rainfall), family=Gamma(link="log"), data = H.incana.2022)
#AIC = 4919.909, BIC = 4973.154 <- this model has a lower AIC but still has outlier issues
P.8 <- glmmTMB(wHeight.cm ~ Distance.n*Days.Elapsed*abundance + (1|Nest) +(1|rainfall), family=Gamma(link="log"), data = H.incana.2022)
#AIC = 4923.633, BIC = 5003.5

simulationOutput <- simulateResiduals(fittedModel = P.4)
plot(simulationOutput)

emmeans(P.4, pairwise ~ Distance.n)

summary(P.4)
Anova(P.4)

####### height 2022 plots #########

interact_plot(P.4, pred = Days.Elapsed, modx = Distance.n, plot.points = FALSE, 
                      interval = TRUE, legend.main = "Distance (m)") + 
  ggtitle("(a) Interaction Plot: Days \nElapsed and Distance (2022)") +
  labs(
    x = "Days Elapsed", 
    y = "Winsorized height (cm)") +
  theme(
    plot.title = element_text(size = 12, face = "bold", family = "Arial", color = "black"),
    axis.title.x = element_text(size = 10, family = "Arial", color = "black"),
    axis.title.y = element_text(size = 10, family = "Arial", color = "black"),
    axis.text = element_text(size = 10, family = "Arial", color = "black"),
    legend.title = element_text(size = 10, family = "Arial", color = "black"),
    legend.text = element_text(size = 10, family = "Arial", color = "black"),
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank()   # Removes minor grid lines
  )


#density X height

interact_plot(P.4, pred = abundance, modx = Distance.n, plot.points = FALSE, 
                       interval = TRUE, legend.main = "Distance (m)") + 
  ggtitle("(b) Interaction Plot: Density \nand Distance (2022)") +
  labs(
    x = expression(italic("H. incana") ~ "Density"), 
    y = "Winsorized height (cm)") +
  theme(
    plot.title = element_text(size = 12, face = "bold", family = "Arial", color = "black"),
    axis.title.x = element_text(size = 10, family = "Arial", color = "black"),
    axis.title.y = element_text(size = 10, family = "Arial", color = "black"),
    axis.text = element_text(size = 10, family = "Arial", color = "black"),
    legend.title = element_text(size = 10, family = "Arial", color = "black"),
    legend.text = element_text(size = 10, family = "Arial", color = "black"),
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank()   # Removes minor grid lines
  )


###### presence of flowers model 2022 #######

clean_flow2 <- H.incana.2022[complete.cases(H.incana.2022[, "Flowering"]), ]

#NOTE: THESE ARE NOT THE MODELS USED IN THE PAPER, KEEP SCROLLING TO SEE WORKFLOW TO DEAL WITH COLINERITY
#comparing AICs among models
Q.1 <- glmmTMB(Flowering ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 325.5566, BIC = 352.0319
Q.2 <- glmmTMB(Flowering ~ Days.Elapsed + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 299.6296, BIC = 320.8098
Q.3 <- glmmTMB(Flowering ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 312.5132, BIC = 333.6935
Q.4 <- glmmTMB(Flowering ~ Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.2106, BIC = 194.3909
Q.5 <- glmmTMB(Flowering ~ Distance.n + Days.Elapsed + abundance + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.9205, BIC = 216.281
Q.6 <- glmmTMB(Flowering ~ Days.Elapsed + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.2215, BIC = 199.6968
Q.7 <- glmmTMB(Flowering ~ Distance.n + abundance + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.9933, BIC = 211.0587
Q.8 <- glmmTMB(Flowering ~ Days.Elapsed*Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 148.1908, BIC = 179.9611 <- selected model
Q.9 <- glmmTMB(Flowering ~ Distance.n + Days.Elapsed*Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 147.364, BIC = 189.7245
Q.10 <- glmmTMB(Flowering ~ Distance.n*Days.Elapsed + Days.Elapsed*Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 144.8363, BIC = 197.7868 <- newly selected model

simulationOutput <- simulateResiduals(fittedModel = Q.8)
plot(simulationOutput)

# Extract the fixed effects formula
fixed_effects_formula <- reformulate(attr(terms(Q.8), "term.labels"), response = as.character(formula(Q.8)[[2]]))

# Fit the linear model with the same fixed effects
lm_model <- lm(fixed_effects_formula, data = clean_flow2)

# Calculate the VIF values
vif_values <- vif(lm_model)

# Print VIF values
print(vif_values) #height is at 16

#making a correlation matrix
# Creating a correlation matrix
cor_matrix <- cor(clean_flow2[c("Days.Elapsed", "Height.cm")])
#height and days elapsed have a slightly heigh positive correlation of 0.33

########I will now regress height and add it back to the model as residuals
# Fit a linear model
height_model2 <- lm(Height.cm ~ Days.Elapsed, data = clean_flow2)

# Extract the residuals
clean_flow2$Height_residuals <- residuals(height_model2)

Q.1 <- glmmTMB(Flowering ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 325.5566, BIC = 352.0319
Q.2 <- glmmTMB(Flowering ~ Days.Elapsed + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 299.6296, BIC = 320.8098
Q.3 <- glmmTMB(Flowering ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 312.5132, BIC = 333.6935
Q.4 <- glmmTMB(Flowering ~ Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 175.2178, BIC = 196.398 went up by 2
Q.5 <- glmmTMB(Flowering ~ Distance.n + Days.Elapsed + abundance + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.9205, BIC = 216.281 stayed the same
Q.6 <- glmmTMB(Flowering ~ Days.Elapsed + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.2215, BIC = 199.6968 stayed the same
Q.7 <- glmmTMB(Flowering ~ Distance.n + abundance + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 176.0059, BIC = 213.0713 went up by 2
Q.8 <- glmmTMB(Flowering ~ Days.Elapsed*Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 148.7697, BIC = 180.5401 <- selected modelwent up slightly
Q.9 <- glmmTMB(Flowering ~ Distance.n + Days.Elapsed*Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 147.9462, BIC = 190.3066 went up slightly
Q.10 <- glmmTMB(Flowering ~ Distance.n*Days.Elapsed + Days.Elapsed*Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 145.4443, BIC = 198.3948 <- newly selected model went up slightly


simulationOutput <- simulateResiduals(fittedModel = Q.8)
plot(simulationOutput)

#IMPORTANT NOTE: AT THIS POINT WORRIED ABOUT COLINERITY, HERE WAS MY WORK-FLOW TO SOLVE:

# Extract the fixed effects formula
fixed_effects_formula <- reformulate(attr(terms(Q.8), "term.labels"), response = as.character(formula(Q.8)[[2]]))

# Fit the linear model with the same fixed effects
lm_model <- lm(fixed_effects_formula, data = clean_flow2)

# Calculate the VIF values
vif_values <- vif(lm_model)

# Print VIF values
print(vif_values) #VIF value only went down slightly 16 -> 15

cor_matrix <- cor(clean_flow2[c("Days.Elapsed", "Height_residuals")]) #still have a neg correlation

####HERE IS WHERE I DECIDED TO USE 'CENTERING', THESE ARE THE MODELS IN THE PAPER
#centering
clean_flow2$Height_centered <- scale(clean_flow2$Height.cm, center = TRUE, scale = FALSE)
clean_flow2$Days_elapsed_centered <- scale(clean_flow2$Days.Elapsed, center = TRUE, scale = FALSE)


Q.1 <- glmmTMB(Flowering ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 325.5566, BIC = 352.0319
Q.2 <- glmmTMB(Flowering ~ Days_elapsed_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 299.6296, BIC = 320.8098
Q.3 <- glmmTMB(Flowering ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 312.5132, BIC = 333.6935
Q.4 <- glmmTMB(Flowering ~ Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.2106, BIC = 194.3909
Q.5 <- glmmTMB(Flowering ~ Distance.n + Days_elapsed_centered + abundance + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.9205, BIC = 216.281
Q.6 <- glmmTMB(Flowering ~ Days_elapsed_centered + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.2215, BIC = 199.6968
Q.7 <- glmmTMB(Flowering ~ Distance.n + abundance + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 173.9933, BIC = 211.0587
Q.8 <- glmmTMB(Flowering ~ Days_elapsed_centered*Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#AIC = 148.1908, BIC = 179.9611 <- selected model
Q.9 <- glmmTMB(Flowering ~ Distance.n + Days_elapsed_centered*Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#no convergence
Q.10 <- glmmTMB(Flowering ~ Distance.n*Days_elapsed_centered + Days_elapsed_centered*Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow2)
#no convergence


simulationOutput <- simulateResiduals(fittedModel = Q.8)
plot(simulationOutput)

summary(Q.8)
Anova(Q.8)

# Check VIF values for the centered model
fixed_effects_formula_centered <- reformulate(attr(terms(Q.8), "term.labels"), response = as.character(formula(Q.8)[[2]]))
lm_model_centered <- lm(fixed_effects_formula_centered, data = clean_flow2)
vif_values_centered <- vif(lm_model_centered, type = "predictor")
print(vif_values_centered) #vif is now 1


####### load data ########
H.incana.2023.raw <- read.xlsx("H.incana_phenology_2023_only.xlsx", colNames = TRUE, rowNames = FALSE, detectDates = TRUE)

####### data clean-up #######

#changing column names
colnames(H.incana.2023.raw)[c(5:6)] <- c("Distance.m","Height.cm") # fix the weird names
colnames(H.incana.2023.raw)[c(12)] <- c("rainfall")

#converting nest ID to character
H.incana.2023.raw$Nest <- as.character(H.incana.2023.raw$Nest)

#add H. incana abundance to the data set
options(dplyr.summarise.inform = FALSE)

H.incana.abun <- H.incana.2023.raw %>%
  group_by(Nest, Date, Degree, Distance.m) %>%
  summarise(abundance=length(Height.cm))

H.incana.2023 <- merge(H.incana.2023.raw, H.incana.abun, by = c("Nest", "Date", "Degree", "Distance.m"))

#Convert relevant yes/no to 1/0
H.incana.flowering <- H.incana.2023$Flowering<-ifelse(H.incana.2023$Flowering=="yes",1,0)
H.incana.seed <- H.incana.2023$Set.Seed<-ifelse(H.incana.2023$Set.Seed=="yes",1,0)
H.incana.senesced <- H.incana.2023$Senesced<-ifelse(H.incana.2023$Senesced=="yes",1,0)

#duplicate distance and make it an ordered factor
H.incana.2023$Distance.n = H.incana.2023$Distance.m
H.incana.2023$Distance.n <- factor(H.incana.2023$Distance.n, ordered = TRUE, 
                                   levels = c("0", "2.5", "5"))


####### height 2023 model decision making #######

#looking at our data distribution
hist(H.incana.2023.raw$Height.cm) #this does not look normal to me
#density plot
ggdensity(H.incana.2023.raw$Height.cm, 
          main = "Density plot of height",
          xlab = "height") #does not look bell shaped
#Q-Q plot
ggqqplot(H.incana.2023.raw$Height.cm) 

#shapiro test
shapiro.test(H.incana.2023.raw$Height.cm)
#indeed, the dist is non-normal

gamma_fit(H.incana.2023.raw$Height.cm)
gamma_test(H.incana.2023.raw$Height.cm)
descdist(H.incana.2023.raw$Height.cm)#Height has a student's T distribution

###### height 2023 model #######

clean_height <- H.incana.2023[complete.cases(H.incana.2023[, "Height.cm"]), ] #remove blank entries

#comparing AICs among models
H.1 <- glmmTMB(Height.cm ~ Days.elapsed + (1|Nest) +(1|rainfall),
               family=t_family(link="identity"), data = clean_height)
#AIC = 16774.63, BIC = 16807.8
H.2 <- glmmTMB(Height.cm ~ Days.elapsed + abundance + (1|Nest) +(1|rainfall),
               family=t_family(link="identity"), data = clean_height)
#AIC = 16756.74, BIC = 16795.44
H.3 <- glmmTMB(Height.cm ~ Distance.n + Days.elapsed + (1|Nest) +(1|rainfall),
               family=t_family(link="identity"), data = clean_height)
#AIC = 16765.88, BIC = 16810.1
H.4 <- glmmTMB(Height.cm ~ Distance.n*Days.elapsed + (1|Nest) +(1|rainfall),
               family=t_family(link="identity"), data = clean_height)
#AIC = 16753.89, BIC = 16809.17
H.5 <- glmmTMB(Height.cm ~ Distance.n + Days.elapsed + abundance + (1|Nest) +(1|rainfall),
               family=t_family(link="identity"), data = clean_height)
#AIC = 16750.11, BIC = 16799.86
H.6 <- glmmTMB(Height.cm ~ Distance.n*Days.elapsed + abundance + (1|Nest) +(1|rainfall),
               family=t_family(link="identity"), data = clean_height)
#AIC = 16742.17, BIC = 16802.98
H.7 <- glmmTMB(Height.cm ~ Distance.n + Days.elapsed*abundance + (1|Nest) +(1|rainfall),
               family=t_family(link="identity"), data = clean_height)
#AIC = 16710.14, BIC = 16765.42
H.8 <- glmmTMB(Height.cm ~Days.elapsed*abundance + (1|Nest) +(1|rainfall),
               family=t_family(link="identity"),  data = clean_height)
#AIC = 16712.61, BIC = 16756.83

simulationOutput <- simulateResiduals(fittedModel = H.8)
plot(simulationOutput)

summary(H.8)
Anova(H.8)


####### presence of flowers model 2023 ########

clean_flow <- H.incana.2023[complete.cases(H.incana.2023[, "Flowering"]), ]

#NOTE: THESE ARE NOT THE MODELS USED IN THE PAPER, KEEP SCROLLING TO SEE WORKFLOW TO DEAL WITH COLINERITY
#comparing AICs among models
F.1 <- glmmTMB(Flowering ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1353.241, BIC = 1380.603
F.2 <- glmmTMB(Flowering ~ Days.elapsed + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1215.39, BIC = 1237.28
F.3 <- glmmTMB(Flowering ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1341.288, BIC = 1363.178 <-
F.4 <- glmmTMB(Flowering ~ Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1114.302, BIC = 1136.192
F.5 <- glmmTMB(Flowering ~ Distance.n + Days.elapsed + abundance + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1116.606, BIC = 1160.386 <-
F.6 <- glmmTMB(Flowering ~ Days.elapsed + abundance + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1115.516, BIC = 1148.351 <-
F.7 <- glmmTMB(Flowering ~ Distance.n + Days.elapsed + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1115.626, BIC = 1153.933
F.8 <- glmmTMB(Flowering ~ Days.elapsed + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1114.552, BIC = 1141.915
F.9 <- glmmTMB(Flowering ~ Days.elapsed*Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 825.3307, BIC = 858.1657
F.10 <- glmmTMB(Flowering ~ abundance + Days.elapsed*Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 822.749, BIC = 861.0565 
F.11 <- glmmTMB(Flowering ~ Distance.n + Days.elapsed*Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 826.64, BIC = 870.42
F.12 <- glmmTMB(Flowering ~ Distance.n*Days.elapsed+ Days.elapsed*Height.cm + abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 814.8886, BIC = 875.0861 <- selected model

simulationOutput <- simulateResiduals(fittedModel = F.12)
plot(simulationOutput)

emmeans(F.12, pairwise ~ Distance.n) 

#IMPORTANT NOTE: AT THIS POINT WORRIED ABOUT COLINERITY, HERE WAS MY WORK-FLOW TO SOLVE:

# Extract the fixed effects formula
fixed_effects_formula <- reformulate(attr(terms(F.12), "term.labels"), response = as.character(formula(F.12)[[2]]))

# Fit the linear model with the same fixed effects
lm_model <- lm(fixed_effects_formula, data = clean_flow)

# Calculate the VIF values
vif_values <- vif(lm_model)

# Print VIF values
print(vif_values)

#making a correlation matrix
# Creating a correlation matrix
cor_matrix <- cor(clean_flow[c("Days.elapsed", "Height.cm", "abundance")])
#height and days elapsed have a strong positive correlation of about 0.84

########I will now regress height and add it back to the model as residuals
# Fit a linear model
height_model <- lm(Height.cm ~ Days.elapsed, data = clean_flow)

# Extract the residuals
clean_flow$Height_residuals <- residuals(height_model)

F.1 <- glmmTMB(Flowering ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1353.241, BIC = 1380.603
F.2 <- glmmTMB(Flowering ~ Days.elapsed + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1215.39, BIC = 1237.28
F.3 <- glmmTMB(Flowering ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1341.288, BIC = 1363.178 <-
F.4 <- glmmTMB(Flowering ~ Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1259.241 #went up and is now higher than F.3
F.5 <- glmmTMB(Flowering ~ Distance.n + Days.elapsed + abundance + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1116.606 stayed the same
F.6 <- glmmTMB(Flowering ~ Days.elapsed + abundance + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1115.516 STAYED THE SAME
F.7 <- glmmTMB(Flowering ~ Distance.n + Days.elapsed + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1115.626 STAYED THE SAME
F.8 <- glmmTMB(Flowering ~ Days.elapsed + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1114.552, STAYED THE SAME
F.9 <- glmmTMB(Flowering ~ Days.elapsed*Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1055.782 <- best model!!!  
F.10 <- glmmTMB(Flowering ~ abundance + Days.elapsed*Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1057.479 
F.11 <- glmmTMB(Flowering ~ Distance.n + Days.elapsed*Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1057.691
F.12 <- glmmTMB(Flowering ~ Distance.n*Days.elapsed+ Days.elapsed*Height_residuals + abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1059.504


simulationOutput <- simulateResiduals(fittedModel = F.9)
plot(simulationOutput)

# Extract the fixed effects formula
fixed_effects_formula <- reformulate(attr(terms(F.9), "term.labels"), response = as.character(formula(F.9)[[2]]))

# Fit the linear model with the same fixed effects
lm_model <- lm(fixed_effects_formula, data = clean_flow)

# Calculate the VIF values
vif_values2 <- vif(lm_model)

# Print VIF values
print(vif_values2)

cor_matrix2 <- cor(clean_flow[c("Days.elapsed", "Height_residuals")])

####HERE IS WHERE I DECIDED TO USE 'CENTERING', THESE ARE THE MODELS IN THE PAPER
#centering
clean_flow$Height_centered <- scale(clean_flow$Height.cm, center = TRUE, scale = FALSE)
clean_flow$Days_elapsed_centered <- scale(clean_flow$Days.elapsed, center = TRUE, scale = FALSE)


F.1 <- glmmTMB(Flowering ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1353.241, BIC = 1380.603
F.2 <- glmmTMB(Flowering ~ Days_elapsed_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1215.39, BIC = 1237.28
F.3 <- glmmTMB(Flowering ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1341.288, BIC = 1363.178 <-
F.4 <- glmmTMB(Flowering ~ Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1114.302, BIC = 1136.192
F.5 <- glmmTMB(Flowering ~ Distance.n + Days_elapsed_centered + abundance + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1116.606, BIC = 1160.386
F.6 <- glmmTMB(Flowering ~ Days_elapsed_centered + abundance + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1115.516, BIC = 1153.933
F.7 <- glmmTMB(Flowering ~ Distance.n + Days_elapsed_centered + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1115.626, BIC = 1153.933
F.8 <- glmmTMB(Flowering ~ Days_elapsed_centered + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 1114.552, BIC = 1141.915
F.9 <- glmmTMB(Flowering ~ Days_elapsed_centered*Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 825.3307, BIC =  858.1657
F.10 <- glmmTMB(Flowering ~ abundance + Days_elapsed_centered*Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 822.749, BIC = 861.0565
F.11 <- glmmTMB(Flowering ~ Distance.n + Days_elapsed_centered*Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 826.64, BIC = 870.42
F.12 <- glmmTMB(Flowering ~ Distance.n*Days_elapsed_centered + Days_elapsed_centered*Height_centered + abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_flow)
#AIC = 814.8886, BIC = 875.0861


simulationOutput <- simulateResiduals(fittedModel = F.12_centered)
plot(simulationOutput)

summary(F.12)
Anova(F.12)

emmeans(F.12, pairwise ~ Distance.n)

# Check VIF values for the centered model
fixed_effects_formula_centered <- reformulate(attr(terms(F.12), "term.labels"), response = as.character(formula(F.12)[[2]]))
lm_model_centered <- lm(fixed_effects_formula_centered, data = clean_flow)
vif_values_centered <- vif(lm_model_centered, type = "predictor")
print(vif_values_centered)

cor_matrix3 <- cor(clean_flow[c("Days_elapsed_centered", "Height_centered")])


####### presence of flowers 2023 plot ########

plot2 <- interact_plot(F.12, pred = Days_elapsed_centered, modx = Distance.n, plot.points = FALSE, 
                       interval = TRUE, legend.main = "Distance (m)") + 
  ggtitle("(c) Interaction Plot: Days \nElapsed and Distance (2023)") +
  labs(
    x = "Days Elapsed", 
    y = "% Flowering") +
  theme(
    plot.title = element_text(size = 12, face = "bold", family = "Arial", color = "black"),
    axis.title.x = element_text(size = 10, family = "Arial", color = "black"),
    axis.title.y = element_text(size = 10, family = "Arial", color = "black"),
    axis.text = element_text(size = 10, family = "Arial", color = "black"),
    legend.title = element_text(size = 10, family = "Arial", color = "black"),
    legend.text = element_text(size = 10, family = "Arial", color = "black"),
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank()   # Removes minor grid lines
  )


######## presence of seeds 2023 model ########

clean_seed <- H.incana.2023[complete.cases(H.incana.2023[, "Set.Seed"]), ]


#NOTE: THESE ARE NOT THE MODELS USED IN THE PAPER, KEEP SCROLLING TO SEE WORKFLOW TO DEAL WITH COLINERITY

S.1 <- glmmTMB(Set.Seed ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 1079.349, BIC = 1106.712
S.2 <- glmmTMB(Set.Seed ~ Days.elapsed + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 378.2842, BIC = 400.1742
S.3 <- glmmTMB(Set.Seed ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 1081.241, BIC = 1103.131
S.4 <- glmmTMB(Set.Seed ~ Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 509.8007, BIC = 531.6907
S.5 <- glmmTMB(Set.Seed ~ Distance.n + Days.elapsed + abundance + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 240.4627, BIC = 284.2427
S.6 <- glmmTMB(Set.Seed ~ Days.elapsed + abundance + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 245.3173, BIC = 278.1523
S.7 <- glmmTMB(Set.Seed ~ Distance.n + Days.elapsed + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 238.9697, BIC = 277.2772
S.8 <- glmmTMB(Set.Seed ~ Days.elapsed + Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#no convergence, need to add a model here with terms interacting
S.9 <- glmmTMB(Set.Seed ~ Days.elapsed*Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 230.9139, BIC = 263.7489
S.10 <- glmmTMB(Set.Seed ~ Distance.n + Days.elapsed*Height.cm + abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 227.7391, BIC = 276.9916
S.11 <- glmmTMB(Set.Seed ~ Distance.n + Days.elapsed*Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 225.8641, BIC = 269.6441 <- selected model
S.12 <- glmmTMB(Set.Seed ~ Distance.n*Days.elapsed + Days.elapsed*Height.cm + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 225.9503, #BIC = 280.6753


simulationOutput <- simulateResiduals(fittedModel = S.11)
plot(simulationOutput)

emmeans(S.11, pairwise ~ Distance.n)

#IMPORTANT NOTE: AT THIS POINT WORRIED ABOUT COLINERITY, HERE WAS MY WORK-FLOW TO SOLVE:

# Extract the fixed effects formula
fixed_effects_formula <- reformulate(attr(terms(S.11), "term.labels"), response = as.character(formula(S.11)[[2]]))

# Fit the linear model with the same fixed effects
lm_model <- lm(fixed_effects_formula, data = clean_seed)

# Calculate the VIF values
vif_values <- vif(lm_model)

# Print VIF values
print(vif_values)

cor_matrix <- cor(clean_seed[c("Days.elapsed", "Height.cm")])
#high positive correlation at 0.84

########I will now regress height and add it back to the model as residuals
# Fit a linear model
height_model <- lm(Height.cm ~ Days.elapsed, data = clean_seed)

# Extract the residuals
clean_seed$Height_residuals <- residuals(height_model)

S.1 <- glmmTMB(Set.Seed ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 1079.349, BIC = 1106.712
S.2 <- glmmTMB(Set.Seed ~ Days.elapsed + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 378.2842, BIC = 400.1742
S.3 <- glmmTMB(Set.Seed ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 1081.241, BIC = 1103.131
S.4 <- glmmTMB(Set.Seed ~ Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 993.3111, BIC = 1015.201
S.5 <- glmmTMB(Set.Seed ~ Distance.n + Days.elapsed + abundance + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 240.4627, BIC = 284.2427
S.6 <- glmmTMB(Set.Seed ~ Days.elapsed + abundance + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 245.3173, BIC = 278.1523
S.7 <- glmmTMB(Set.Seed ~ Distance.n + Days.elapsed + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 238.9697, BIC = 277.2772
S.8 <- glmmTMB(Set.Seed ~ Days.elapsed + Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#no convergence, need to add a model here with terms interacting
S.9 <- glmmTMB(Set.Seed ~ Days.elapsed*Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 230.6037, BIC = 263.4387
S.10 <- glmmTMB(Set.Seed ~ Distance.n + Days.elapsed*Height_residuals + abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 226.8691, BIC = 276.1216
S.11 <- glmmTMB(Set.Seed ~ Distance.n + Days.elapsed*Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 224.8811, BIC = 268.6611 <- selected model
S.12 <- glmmTMB(Set.Seed ~ Distance.n*Days.elapsed + Days.elapsed*Height_residuals + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 223.997, BIC = 278.722


simulationOutput <- simulateResiduals(fittedModel = S.11)
plot(simulationOutput) #residuals actually look better than before

# Extract the fixed effects formula
fixed_effects_formula <- reformulate(attr(terms(S.11), "term.labels"), response = as.character(formula(S.11)[[2]]))

# Fit the linear model with the same fixed effects
lm_model <- lm(fixed_effects_formula, data = clean_seed)

# Calculate the VIF values
vif_values <- vif(lm_model)

# Print VIF values
print(vif_values) #values lowered to 2.8

cor_matrix <- cor(clean_seed[c("Days.elapsed", "Height_residuals")])
#matrix is showing extremely large negative values, may want to try centering

#CENTERING: MODELS USED IN THE PAPER
clean_seed$Height_centered <- scale(clean_seed$Height.cm, center = TRUE, scale = FALSE)
clean_seed$Days_elapsed_centered <- scale(clean_seed$Days.elapsed, center = TRUE, scale = FALSE)

S.1 <- glmmTMB(Set.Seed ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 1079.349, BIC = 1106.712
S.2 <- glmmTMB(Set.Seed ~ Days_elapsed_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 378.2842, BIC = 400.1742
S.3 <- glmmTMB(Set.Seed ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 1081.241, BIC = 1103.131
S.4 <- glmmTMB(Set.Seed ~ Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 509.8007, BIC = 531.6907
S.5 <- glmmTMB(Set.Seed ~ Distance.n + Days_elapsed_centered + abundance + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 240.4627, BIC = 284.2427
S.6 <- glmmTMB(Set.Seed ~ Days_elapsed_centered + abundance + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 245.3173, BIC = 278.1523
S.7 <- glmmTMB(Set.Seed ~ Distance.n + Days_elapsed_centered + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 238.9697, BIC = 277.2772
S.8 <- glmmTMB(Set.Seed ~ Days_elapsed_centered + Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 245.027, BIC = 272.3895
S.9 <- glmmTMB(Set.Seed ~ Days_elapsed_centered*Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 230.9139, BIC = 263.7489
S.10 <- glmmTMB(Set.Seed ~ Distance.n + Days_elapsed_centered*Height_centered + abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#no convergence
S.11 <- glmmTMB(Set.Seed ~ Distance.n + Days_elapsed_centered*Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 225.8641, BIC = 269.6441 <- selected model
S.12 <- glmmTMB(Set.Seed ~ Distance.n*Days_elapsed_centered + Days_elapsed_centered*Height_centered + (1|Nest) + (1|rainfall), family = "binomial", data = clean_seed)
#AIC = 225.9503, BIC = 280.6753


simulationOutput <- simulateResiduals(fittedModel = S.11)
plot(simulationOutput) #redisuals look the same as the original model

emmeans(S.11, pairwise ~ Distance.n)

summary(S.11)
Anova(S.11)

# Check VIF values for the centered model
fixed_effects_formula_centered <- reformulate(attr(terms(S.11), "term.labels"), response = as.character(formula(S.11)[[2]]))
lm_model_centered <- lm(fixed_effects_formula_centered, data = clean_seed)
vif_values_centered <- vif(lm_model_centered, type = "predictor")
print(vif_values_centered)

###### Presence of seeds 2023 plot #######

interact_plot(S.11, pred = Days_elapsed_centered, modx = Distance.n, plot.points = FALSE, 
              interval = TRUE, legend.main = "Distance (m)") + 
  ggtitle("(d) Interaction Plot: Days \nElapsed and Distance (2023)") +
  labs(
    x = "Days Elapsed", 
    y = "% Siliques Present") +
  theme(
    plot.title = element_text(size = 12, face = "bold", family = "Arial", color = "black"),
    axis.title.x = element_text(size = 10, family = "Arial", color = "black"),
    axis.title.y = element_text(size = 10, family = "Arial", color = "black"),
    axis.text = element_text(size = 10, family = "Arial", color = "black"),
    legend.title = element_text(size = 10, family = "Arial", color = "black"),
    legend.text = element_text(size = 10, family = "Arial", color = "black"),
    panel.grid.major = element_blank(),  # Removes major grid lines
    panel.grid.minor = element_blank()   # Removes minor grid lines
  )
######### senescence model 2023 ##########

clean_sen <- H.incana.2023[complete.cases(H.incana.2023[, "Senesced"]), ]

#not including height in the model, because I don't think height is a predictor of senescence

SE.1 <- glmmTMB(Senesced ~ Distance.n + (1|Nest) + (1|rainfall), family = "binomial", data = clean_sen)
#AIC = 577.5908, BIC = 605.2325
SE.2 <- glmmTMB(Senesced ~ Days.elapsed + (1|Nest) + (1|rainfall), family = "binomial", data = clean_sen)
#AIC = 397.8201, BIC = 419.9334
SE.3 <- glmmTMB(Senesced ~ abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_sen)
#AIC = 574.0573, BIC = 596.1707
SE.4 <- glmmTMB(Senesced ~ Distance.n + Days.elapsed + abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_sen)
#AIC = 380.029, BIC = 418.7274
SE.5 <- glmmTMB(Senesced ~ Days.elapsed + abundance + (1|Nest) + (1|rainfall), family = "binomial", data = clean_sen)
#AIC = 380.1472, BIC = 407.7889 <- selected model
SE.6 <- glmmTMB(Senesced ~ Distance.n + Days.elapsed + (1|Nest) + (1|rainfall), family = "binomial", data = clean_sen)
#AIC = 400.6964, BIC = 433.8664
SE.7 <- glmmTMB(Senesced ~ Distance.n + Days.elapsed*abundance + (1|Nest) +(1|rainfall) ,family = "binomial", data = clean_sen)
#AIC = 378.5575, BIC = 422.7842
SE.8 <- glmmTMB(Senesced ~ Distance.n*Days.elapsed*abundance + (1|Nest) +(1|rainfall) ,family = "binomial", data = clean_sen)
#AIC = 380.8871, BIC = 458.2838
SE.9 <- glmmTMB(Senesced ~ Days.elapsed*abundance + (1|Nest) +(1|rainfall) ,family = "binomial", data = clean_sen)
#AIC = 379.8377, BIC = 413.0077

simulationOutput <- simulateResiduals(fittedModel = SE.5)
plot(simulationOutput)

summary(SE.5)
Anova(SE.5)




##################################
########### section 3 ############
##################################


######## loading data #########
seed22 <- read.xlsx("H.incana_seedweight_2022_only.xlsx", colNames = TRUE, rowNames = FALSE, detectDates = TRUE)

###### data clean-up #######
#change column names
colnames(seed22)[c(2)] <- c("distance")
colnames(seed22)[c(7)] <- c("weight")

#change nest to a character
seed22$nest <- as.character(seed22$nest)

#change distance to an ordered factor
seed22$distance.n = seed22$distance
seed22$distance.n <- factor(seed22$distance.n, ordered = TRUE, 
                            levels = c("0", "2.5", "5"))

# Summarize the average abundance based on the grouping variables
seed22 <- seed22 %>%
  group_by(nest, distance.n, plant.specimen) %>%
  summarize(avg_weight = mean(weight, na.rm = TRUE))

###### model decision making #########
#looking at our data distribution
hist(seed22$avg_weight) #looks pretty normal aside from a couple of gaps
#density plot
ggdensity(seed22$avg_weight) #looks bell shaped
#Q-Q plot
ggqqplot(seed22$avg_weight) #not too much dispersion
#shapiro test
shapiro.test(seed22$avg_weight) #sig different from normal (0.0002)
#gamma test
gamma_test(seed22$avg_weight) #significantly different from gamma (0.02)

descdist(seed22$avg_weight)


###### seed weight 2022 model #########

model1 <- glmmTMB(avg_weight ~ distance.n + (1|nest), family = Gamma(link = "log"), 
                  data = seed22)

simulationOutput <- simulateResiduals(fittedModel = model1)
plot(simulationOutput) #looks good

summary(model1)

Anova(model1, type = "II", test = "Chisq")

#trying a studen't t-distribution

model3 <- glmmTMB(avg_weight ~ distance.n + (1|nest),
                  family=t_family(link="identity"), data = seed22)
#this model has a much lower AIC than the gamma model

simulationOutput <- simulateResiduals(fittedModel = model3)
plot(simulationOutput) #looks good

summary(model3)

Anova(model3, type = "II", test = "Chisq")

######### loading data ########

seed23 <- read.xlsx("H.incana_seedweight_2023_only.xlsx", colNames = TRUE, rowNames = FALSE, detectDates = TRUE)

####### data clean-up #######

#change column names
colnames(seed23)[c(2)] <- c("distance")
colnames(seed23)[c(7)] <- c("weight")

#change nest to a character
seed23$nest <- as.character(seed23$nest)

#change distance to an ordered factor
seed23$distance.n = seed23$distance
seed23$distance.n <- factor(seed23$distance.n, ordered = TRUE, 
                            levels = c("0", "2.5", "5"))

# Using table() to view the frequency of all levels in the column
table(seed23$distance.n)

# Using sum() to get the count of 2.5 specifically
sum(seed23$distance.n == 2.5)

####### model descision making ##########

#looking at our data distribution
hist(seed23$weight) #looks pretty normal aside from a couple of gaps
#density plot
ggdensity(seed23$weight) #looks mostly bell shaped
#Q-Q plot
ggqqplot(seed23$weight) #not too much dispersion
#shapiro test
shapiro.test(seed23$weight) #not sig diff from normal (0.09)
#gamma test
gamma_test(seed23$weight) #not significantly different from gamma (0.2)

descdist(seed23$weight)#close to logistic

#model

model2 <- glmmTMB(weight ~ distance.n + (1|nest), family = Gamma(link = "log"), 
                  data = seed23)
#AIC = -1115.447

simulationOutput <- simulateResiduals(fittedModel = model2)
plot(simulationOutput) #looks good

summary(model2)

Anova(model2, type = "II", test = "Chisq")

