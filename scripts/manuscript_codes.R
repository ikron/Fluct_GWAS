#Association mapping with GAPIT 3, PCA and quantitative genetics
#for N. crassa growth rate data in Räsänen et al. 2024.

#Load packages
library(dplyr)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
library(gridExtra)
library(signs)
library(scales)
library(grid)
library(tidyr)
library(tidybayes)
library(coda)
library(MASS)
library(colorspace)
library(brms)#Interfacing Stan for Bayesian models

#Quantitative genetics:

#Data with all replicates is used in a multivariate model. Strains that are
#not measured in all treatments or didn't grow are marked with NA. There are
#same number of replicates for each strain (NAs added to data if there was not 
#equal amount of replicates measured).

fluct <- read.csv("natpop_fluct_mv.csv", header = T, dec = ".", sep = ",")

head(fluct)

str(fluct)

#Change genotype "Genot", temperature "Temp" and speed of fluctuations
#"Frequency" as factors.

fluct$Frequency <- factor(fluct$Frequency)
fluct$Temp <- factor(fluct$Temp)
fluct$Genot <- factor(fluct$Genot)

str(fluct)

#Plot raw data for phenotypes measured at flutuating temperatures

fluctonly <- filter(fluct, Temp %in% c("25-35", "32-42"))

labels <- as_labeller(c(`25-35` = "25-35 °C", `32-42` = "32-42 °C"))

rawdata <- ggplot(data=fluctonly, aes(x=Frequency, y=growthrate, group=Genot)) +
  geom_point(alpha = 0.1)+
  geom_line(alpha = 0.1)+
  facet_wrap(~ Temp, ncol=2, labeller = labels)+
  ylab("Growth rate (mm/h)")+
  scale_x_discrete(labels = c("Fast", "Slow"))+
  theme_classic()+
  theme(strip.text.x = element_text(size =14))+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))

print(rawdata)

#Multivariate data analysis with brms. In multivariate analysis you can set multiple
#factors to be explained. Note that the growth rate at each temperature are thought 
#as a separate trait.

#Select only columns that are needed

fluct_cut <- fluct[,c(2,3,13,15)]

gr25 <- filter(fluct_cut, Temp == "25")
gr30 <- filter(fluct_cut, Temp == "30")
gr32 <- filter(fluct_cut, Temp == "32")
gr35 <- filter(fluct_cut, Temp == "35")
gr37 <- filter(fluct_cut, Temp == "37")
gr42 <- filter(fluct_cut, Temp == "42")
gr2535F <- subset(fluct_cut, Temp == "25-35" & Frequency == "120")
gr2535S <- subset(fluct_cut, Temp == "25-35" & Frequency == "480")
gr3242F <- subset(fluct_cut, Temp == "32-42" & Frequency == "120")
gr3242S <- subset(fluct_cut, Temp == "32-42" & Frequency == "480")

#Combing growth rates to create the final multivariate data set
datamv <- data.frame(Genot = gr25$Genot, gr25 = gr25$growthrate, gr30 = gr30$growthrate, gr32 = gr32$growthrate, gr35 = gr35$growthrate,
                    gr37 = gr37$growthrate, gr42 = gr42$growthrate, gr2535F = gr2535F$growthrate, gr2535S = gr2535S$growthrate, 
                   gr3242F = gr3242F$growthrate, gr3242S = gr3242S$growthrate)

#Find maximum and minimum growth rates at each treatment

head(datamv)

Rawdata <- matrix(rep(0, 10*4), ncol = 4)
colnames(Rawdata) <- c("Temperature", "mean", "min", "max")
Rawdata[1,1] <- "25"
Rawdata[2,1] <- "30"
Rawdata[3,1] <- "32"
Rawdata[4,1] <- "35"
Rawdata[5,1] <- "37"
Rawdata[6,1] <- "42"
Rawdata[7,1] <- "25-35 Fast"
Rawdata[8,1] <- "25-35 Slow"
Rawdata[9,1] <- "32-42 Fast"
Rawdata[10,1] <- "32-42 Slow"

Rawdata[1,2] <- mean(datamv$gr25, na.rm = TRUE)
Rawdata[1,3] <- min(datamv$gr25, na.rm = TRUE)
Rawdata[1,4] <- max(datamv$gr25, na.rm = TRUE)

Rawdata[2,2] <- mean(datamv$gr30, na.rm = TRUE)
Rawdata[2,3] <- min(datamv$gr30, na.rm = TRUE)
Rawdata[2,4] <- max(datamv$gr30, na.rm = TRUE)

Rawdata[3,2] <- mean(datamv$gr32, na.rm = TRUE)
Rawdata[3,3] <- min(datamv$gr32, na.rm = TRUE)
Rawdata[3,4] <- max(datamv$gr32, na.rm = TRUE)

Rawdata[4,2] <- mean(datamv$gr35, na.rm = TRUE)
Rawdata[4,3] <- min(datamv$gr35, na.rm = TRUE)
Rawdata[4,4] <- max(datamv$gr35, na.rm = TRUE)

Rawdata[5,2] <- mean(datamv$gr37, na.rm = TRUE)
Rawdata[5,3] <- min(datamv$gr37, na.rm = TRUE)
Rawdata[5,4] <- max(datamv$gr37, na.rm = TRUE)

Rawdata[6,2] <- mean(datamv$gr42, na.rm = TRUE)
Rawdata[6,3] <- min(datamv$gr42, na.rm = TRUE)
Rawdata[6,4] <- max(datamv$gr42, na.rm = TRUE)

Rawdata[7,2] <- mean(datamv$gr2535F, na.rm = TRUE)
Rawdata[7,3] <- min(datamv$gr2535F, na.rm = TRUE)
Rawdata[7,4] <- max(datamv$gr2535F, na.rm = TRUE)

Rawdata[8,2] <- mean(datamv$gr2535S, na.rm = TRUE)
Rawdata[8,3] <- min(datamv$gr2535S, na.rm = TRUE)
Rawdata[8,4] <- max(datamv$gr2535S, na.rm = TRUE)

Rawdata[9,2] <- mean(datamv$gr3242F, na.rm = TRUE)
Rawdata[9,3] <- min(datamv$gr3242F, na.rm = TRUE)
Rawdata[9,4] <- max(datamv$gr3242F, na.rm = TRUE)

Rawdata[10,2] <- mean(datamv$gr3242S, na.rm = TRUE)
Rawdata[10,3] <- min(datamv$gr3242S, na.rm = TRUE)
Rawdata[10,4] <- max(datamv$gr3242S, na.rm = TRUE)

Rawdata

#Save data for results section 
#write.csv(Rawdata, file = "Growth rates.csv")

#Model fitting
model1 <- brm(mvbind(gr2535F, gr2535S, gr3242F, gr3242S, gr32, gr37, gr42, gr25, gr30, gr35) ~ 1 + (1|p|Genot), data = datamv, warmup = 1000, iter = 6000, thin = 2, chains = 2, cores = 2)

#Check model output
summary(model1)

#Plot model
plot(model1)

#Check priors
prior_summary(model1)

#save model1 as .RData for later use
save(model1, file= "multivariate_model.RData")
#Load .RData as needed
load("multivariate_model.RData")

##Extract posterior samples
mvpost <- posterior_samples(model1, pars = c("^b", "^sd", "cor", "sigma"))

##Change names to better ones
nimet <- gsub("b_", "", colnames(mvpost))
nimet <- gsub("_Intercept", "", nimet)
nimet <- gsub("_Genot", "", nimet)
nimet <- gsub("__", "_", nimet)
colnames(mvpost) <- nimet

#Check posterior data
head(mvpost)
summary(mvpost)
str(mvpost)

#Posterior calculations for the results section.
#Report posterior mean and 95% HPDI

#Difference within temperature range, between frequencies

#25-35 °C

mean(mvpost$gr2535S-mvpost$gr2535F)

quantile((mvpost$gr2535S-mvpost$gr2535F), probs = c(0.025, 0.975))

#32-42 °C

mean(mvpost$gr3242S-mvpost$gr3242F)

quantile((mvpost$gr3242S-mvpost$gr3242F), probs = c(0.025, 0.975))

#Difference within frequencies, between temperature ranges

#Fast fluctuation

mean(mvpost$gr2535F-mvpost$gr3242F)

quantile((mvpost$gr2535F-mvpost$gr3242F), probs = c(0.025, 0.975))

#Slow fluctuation

mean(mvpost$gr2535S-mvpost$gr3242S)

quantile((mvpost$gr2535S-mvpost$gr3242S), probs = c(0.025, 0.975))

#between optimal and upper extreme

mean(mvpost$gr35-mvpost$gr42)

quantile((mvpost$gr35-mvpost$gr42), probs = c(0.025, 0.975))


#Count heritabilities and save data frame

herit <- matrix(rep(0, 10*4), ncol = 4)

colnames(herit) <- c("Temperature", "mean", "lower", "upper")

herit <- data.frame(herit)

herit[1,1] <- "25"
herit[2,1] <- "30"
herit[3,1] <- "32"
herit[4,1] <- "35"
herit[5,1] <- "37"
herit[6,1] <- "42"
herit[7,1] <- "25-35 Fast"
herit[8,1] <- "25-35 Slow"
herit[9,1] <- "32-42 Fast"
herit[10,1] <- "32-42 Slow"
herit

#25 °C

heritability25 <- mvpost$sd_gr25^2/(mvpost$sd_gr25^2+mvpost$sigma_gr25^2)
herit25 <- quantile(heritability25, probs = c(0.025, 0.975))
herit25 <- as.data.frame(herit25)
herit[1, c(3,4)] <- herit25[c(1,2), 1]
herit[1, 2] <- mean(heritability25)

#30 °C

heritability30 <- mvpost$sd_gr30^2/(mvpost$sd_gr30^2+mvpost$sigma_gr30^2)
herit30 <- quantile(heritability30, probs = c(0.025, 0.975))
herit30 <- as.data.frame(herit30)
herit[2, c(3,4)] <- herit30[c(1,2), 1]
herit[2, 2] <- mean(heritability30)

#32 °C

heritability32 <- mvpost$sd_gr32^2/(mvpost$sd_gr32^2+mvpost$sigma_gr32^2)
herit32 <- quantile(heritability32, probs = c(0.025, 0.975))
herit32 <- as.data.frame(herit32)
herit[3, c(3,4)] <- herit32[c(1,2), 1]
herit[3, 2] <- mean(heritability32)

#35 °C

heritability35 <- mvpost$sd_gr35^2/(mvpost$sd_gr35^2+mvpost$sigma_gr35^2)
herit35 <- quantile(heritability35, probs = c(0.025, 0.975))
herit35 <- as.data.frame(herit35)
herit[4, c(3,4)] <- herit35[c(1,2), 1]
herit[4, 2] <- mean(heritability35)

#37 °C

heritability37 <- mvpost$sd_gr37^2/(mvpost$sd_gr37^2+mvpost$sigma_gr37^2)
herit37 <- quantile(heritability37, probs = c(0.025, 0.975))
herit37 <- as.data.frame(herit37)
herit[5, c(3,4)] <- herit37[c(1,2), 1]
herit[5, 2] <- mean(heritability37)

#42 °C

heritability42 <- mvpost$sd_gr42^2/(mvpost$sd_gr42^2+mvpost$sigma_gr42^2)
herit42 <- quantile(heritability42, probs = c(0.025, 0.975))
herit42 <- as.data.frame(herit42)
herit[6, c(3,4)] <- herit42[c(1,2), 1]
herit[6, 2] <- mean(heritability42)

#25-35 °C fast

heritability2535F <- mvpost$sd_gr2535F^2/(mvpost$sd_gr2535F^2+mvpost$sigma_gr2535F^2)
herit2535F <- quantile(heritability2535F, probs = c(0.025, 0.975))
herit2535F <- as.data.frame(herit2535F)
herit[7, c(3,4)] <- herit2535F[c(1,2), 1]
herit[7, 2] <- mean(heritability2535F)

#25-35 °C slow

heritability2535S <- mvpost$sd_gr2535S^2/(mvpost$sd_gr2535S^2+mvpost$sigma_gr2535S^2)
herit2535S <- quantile(heritability2535S, probs = c(0.025, 0.975))
herit2535S <- as.data.frame(herit2535S)
herit[8, c(3,4)] <- herit2535S[c(1,2), 1]
herit[8, 2] <- mean(heritability2535S)

#32-42 °C fast

heritability3242F <- mvpost$sd_gr3242F^2/(mvpost$sd_gr3242F^2+mvpost$sigma_gr3242F^2)
herit3242F <- quantile(heritability3242F, probs = c(0.025, 0.975))
herit3242F <- as.data.frame(herit3242F)
herit[9, c(3,4)] <- herit3242F[c(1,2), 1]
herit[9, 2] <- mean(heritability3242F)

#32-42 °C slow

heritability3242S <- mvpost$sd_gr3242S^2/(mvpost$sd_gr3242S^2+mvpost$sigma_gr3242S^2)
herit3242S <- quantile(heritability3242S, probs = c(0.025, 0.975))
herit3242S <- as.data.frame(herit3242S)
herit[10, c(3,4)] <- herit3242S[c(1,2), 1]
herit[10, 2] <- mean(heritability3242S)
herit

write.csv(herit, file = "heritabilities.csv")

#Make plot for results

herit$Temperature <- factor(herit$Temperature, levels = c("25", "30", "35", "25-35 Fast", "25-35 Slow", "32", "37", "42", "32-42 Fast", "32-42 Slow"))

heritplot <- ggplot(data=herit, aes(x=Temperature, y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange()+
  ylab("Heritability")+
  xlab("Temperature (°C)") +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6, size = 14),
        axis.title.x = element_text(size = 16))+
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))

print(heritplot)

#Count coefficients of variation

#Genetic variation

varmat <- matrix(rep(0, 10*4), ncol = 4)

colnames(varmat) <- c("Temperature", "mean", "lower","upper")

varmat <- data.frame(varmat)

varmat[1,1] <- "25"
varmat[2,1] <- "30"
varmat[3,1] <- "32"
varmat[4,1] <- "35"
varmat[5,1] <- "37"
varmat[6,1] <- "42"
varmat[7,1] <- "25-35 Fast"
varmat[8,1] <- "25-35 Slow"
varmat[9,1] <- "32-42 Fast"
varmat[10,1] <- "32-42 Slow"
varmat

#25 °C

Coefg25 <- 100*(mvpost$sd_gr25/mvpost$gr25)
CVg25 <- quantile(Coefg25, probs = c(0.025, 0.975))
CVg25 <- as.data.frame(CVg25)
CVg25
varmat[1, c(3,4)] <- CVg25[c(1,2), 1]
varmat[1, 2] <- mean(Coefg25)

#30 °C

Coefg30 <- 100*(mvpost$sd_gr30/mvpost$gr30)
CVg30 <- quantile(Coefg30, probs = c(0.025, 0.975))
CVg30 <- as.data.frame(CVg30)
CVg30
varmat[2, c(3,4)] <- CVg30[c(1,2), 1]
varmat[2, 2] <- mean(Coefg30)

#32 °C

Coefg32 <- 100*(mvpost$sd_gr32/mvpost$gr32)
CVg32 <- quantile(Coefg32, probs = c(0.025, 0.975))
CVg32 <- as.data.frame(CVg32)
CVg32
varmat[3, c(3,4)] <- CVg32[c(1,2), 1]
varmat[3, 2] <- mean(Coefg32)

#35 °C

Coefg35 <- 100*(mvpost$sd_gr35/mvpost$gr35)
CVg35 <- quantile(Coefg35, probs = c(0.025, 0.975))
CVg35 <- as.data.frame(CVg35)
CVg35
varmat[4, c(3,4)] <- CVg35[c(1,2), 1]
varmat[4, 2] <- mean(Coefg35)

#37 °C

Coefg37 <- 100*(mvpost$sd_gr37/mvpost$gr37)
CVg37 <- quantile(Coefg37, probs = c(0.025, 0.975))
CVg37 <- as.data.frame(CVg37)
CVg37
varmat[5, c(3,4)] <- CVg37[c(1,2), 1]
varmat[5, 2] <- mean(Coefg37)

#42 °C

Coefg42 <- 100*(mvpost$sd_gr42/mvpost$gr42)
CVg42 <- quantile(Coefg42, probs = c(0.025, 0.975))
CVg42 <- as.data.frame(CVg42)
CVg42
varmat[6, c(3,4)] <- CVg42[c(1,2), 1]
varmat[6, 2] <- mean(Coefg42)

#25-35 °C fast

Coefg2535F <- 100*(mvpost$sd_gr2535F/mvpost$gr2535F)
CVg2535F <- quantile(Coefg2535F, probs = c(0.025, 0.975))
CVg2535F <- as.data.frame(CVg2535F)
CVg2535F
varmat[7, c(3,4)] <- CVg2535F[c(1,2), 1]
varmat[7, 2] <- mean(Coefg2535F)

#25-35 °C slow

Coefg2535S <- 100*(mvpost$sd_gr2535S/mvpost$gr2535S)
CVg2535S <- quantile(Coefg2535S, probs = c(0.025, 0.975))
CVg2535S <- as.data.frame(CVg2535S)
CVg2535S
varmat[8, c(3,4)] <- CVg2535S[c(1,2), 1]
varmat[8, 2] <- mean(Coefg2535S)

#32-42 °C fast

Coefg3242F <- 100*(mvpost$sd_gr3242F/mvpost$gr3242F)
CVg3242F <- quantile(Coefg3242F, probs = c(0.025, 0.975))
CVg3242F <- as.data.frame(CVg3242F)
CVg3242F
varmat[9, c(3,4)] <- CVg3242F[c(1,2), 1]
varmat[9, 2] <- mean(Coefg3242F)

#32-42 °C slow

Coefg3242S <- 100*(mvpost$sd_gr3242S/mvpost$gr3242S)
CVg3242S <- quantile(Coefg3242S, probs = c(0.025, 0.975))
CVg3242S <- as.data.frame(CVg3242S)
CVg3242S
varmat[10, c(3,4)] <- CVg3242S[c(1,2), 1]
varmat[10, 2] <- mean(Coefg3242S)

varmat

write.csv(varmat, file = "CVg.csv")

#Environmental variation

varmatE <- matrix(rep(0, 10*4), ncol = 4)

colnames(varmatE) <- c("Temperature", "mean", "lower","upper")

varmatE <- data.frame(varmatE)

varmatE[1,1] <- "25"
varmatE[2,1] <- "30"
varmatE[3,1] <- "32"
varmatE[4,1] <- "35"
varmatE[5,1] <- "37"
varmatE[6,1] <- "42"
varmatE[7,1] <- "25-35 Fast"
varmatE[8,1] <- "25-35 Slow"
varmatE[9,1] <- "32-42 Fast"
varmatE[10,1] <- "32-42 Slow"
varmatE

#25 °C

Coefe25 <- 100*(mvpost$sigma_gr25/mvpost$gr25)
CVg25 <- quantile(Coefe25, probs = c(0.025, 0.975))
CVg25 <- as.data.frame(CVg25)
CVg25
varmatE[1, c(3,4)] <- CVg25[c(1,2), 1]
varmatE[1, 2] <- mean(Coefe25)

#30 °C

Coefe30 <- 100*(mvpost$sigma_gr30/mvpost$gr30)
CVg30 <- quantile(Coefe30, probs = c(0.025, 0.975))
CVg30 <- as.data.frame(CVg30)
CVg30
varmatE[2, c(3,4)] <- CVg30[c(1,2), 1]
varmatE[2, 2] <- mean(Coefe30)

#32 °C

Coefe32 <- 100*(mvpost$sigma_gr32/mvpost$gr32)
CVg32 <- quantile(Coefe32, probs = c(0.025, 0.975))
CVg32 <- as.data.frame(CVg32)
CVg32
varmatE[3, c(3,4)] <- CVg32[c(1,2), 1]
varmatE[3, 2] <- mean(Coefe32)

#35 °C

Coefe35 <- 100*(mvpost$sigma_gr35/mvpost$gr35)
CVg35 <- quantile(Coefe35, probs = c(0.025, 0.975))
CVg35 <- as.data.frame(CVg35)
CVg35
varmatE[4, c(3,4)] <- CVg35[c(1,2), 1]
varmatE[4, 2] <- mean(Coefe35)

#37 °C

Coefe37 <- 100*(mvpost$sigma_gr37/mvpost$gr37)
CVg37 <- quantile(Coefe37, probs = c(0.025, 0.975))
CVg37 <- as.data.frame(CVg37)
CVg37
varmatE[5, c(3,4)] <- CVg37[c(1,2), 1]
varmatE[5, 2] <- mean(Coefe37)

#42 °C

Coefe42 <- 100*(mvpost$sigma_gr42/mvpost$gr42)
CVg42 <- quantile(Coefe42, probs = c(0.025, 0.975))
CVg42 <- as.data.frame(CVg42)
CVg42
varmatE[6, c(3,4)] <- CVg42[c(1,2), 1]
varmatE[6, 2] <- mean(Coefe42)

#25-35 °C fast

Coefe2535F <- 100*(mvpost$sigma_gr2535F/mvpost$gr2535F)
CVg2535F <- quantile(Coefe2535F, probs = c(0.025, 0.975))
CVg2535F <- as.data.frame(CVg2535F)
CVg2535F
varmatE[7, c(3,4)] <- CVg2535F[c(1,2), 1]
varmatE[7, 2] <- mean(Coefe2535F)

#25-35 °C slow

Coefe2535S <- 100*(mvpost$sigma_gr2535S/mvpost$gr2535S)
CVg2535S <- quantile(Coefe2535S, probs = c(0.025, 0.975))
CVg2535S <- as.data.frame(CVg2535S)
CVg2535S
varmatE[8, c(3,4)] <- CVg2535S[c(1,2), 1]
varmatE[8, 2] <- mean(Coefe2535S)

#32-42 °C fast

Coefe3242F <- 100*(mvpost$sigma_gr3242F/mvpost$gr3242F)
CVg3242F <- quantile(Coefe3242F, probs = c(0.025, 0.975))
CVg3242F <- as.data.frame(CVg3242F)
CVg3242F
varmatE[9, c(3,4)] <- CVg3242F[c(1,2), 1]
varmatE[9, 2] <- mean(Coefe3242F)

#32-42 °C slow

Coefe3242S <- 100*(mvpost$sigma_gr3242S/mvpost$gr3242S)
CVg3242S <- quantile(Coefe3242S, probs = c(0.025, 0.975))
CVg3242S <- as.data.frame(CVg3242S)
CVg3242S
varmatE[10, c(3,4)] <- CVg3242S[c(1,2), 1]
varmatE[10, 2] <- mean(Coefe3242S)
varmatE

write.csv(varmatE, file = "CVe.csv")

#Plot and combine with heritability

varmat$Temperature <- factor(varmat$Temperature, levels = c("25", "30", "35", 
"25-35 Fast", "25-35 Slow", "32", "37", "42", "32-42 Fast", "32-42 Slow"))

varmatE$Temperature <- factor(varmatE$Temperature, levels = c("25", "30", "35",
"25-35 Fast", "25-35 Slow", "32", "37", "42", "32-42 Fast", "32-42 Slow"))

varmatall <- rbind(varmat, varmatE)
varmatall$type <- factor(c(rep("CVg", 10), rep("CVe", 10)))

varmatall

varcoefplot <-  ggplot(data = varmatall, aes(x = Temperature, y = mean, ymin = lower, ymax = upper)) +
  geom_pointrange() +
  ylab("Coefficient of variation") +
  xlab("Temperature (°C)") +
  facet_grid(type ~ .) +
  theme_classic()+
  theme(strip.text.y = element_text(size =16))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6, size = 14),
        axis.title.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))
  
print(varcoefplot)

fullplot <- plot_grid(heritplot, varcoefplot, labels= c("A", "B"), nrow=1, ncol=2)

save_plot("herit_CV_facet.png", fullplot, base_height = 6, base_width = 17)

#Genetic variances, covariances, correlations, and environmental variances for 
#growth rate estimated from the posteriors of multivariate model

#95 % HPDI

mvpost_quantile <- t(apply(mvpost, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975)))

mvpost_quantile <- round(mvpost_quantile, digits = 2)

mvpost_quantile <- as.data.frame(mvpost_quantile)

#posterior means

mvpost_mean <- apply(mvpost, MARGIN = 2, FUN=mean)

mvpost_mean <- round(mvpost_mean, digits = 2)

mvpost_mean <- as.data.frame(mvpost_mean)

#combined data for posterior means and quantiles

mvpost_all <- cbind(mvpost_mean, mvpost_quantile)

mvpost_all

#Plot posterior estimates and combine with plot from raw data 

#25-35 °C
fluct2535 <- data.frame(mvpost_all[c(1,2),])
colnames(fluct2535) <- c("mean", "lower", "upper")
fluct2535$freq <- c(120, 480)

const2535 <- data.frame(mvpost_all[c(8,9,10),])
colnames(const2535) <- c("mean", "lower", "upper")

const2535mean <- data.frame(xstart = c(rep(120,3)), xfin = c(rep(480,3)), ystart = const2535[,1], yfin = const2535[,1])
const2535lower <- data.frame(xstart = c(rep(120,3)), xfin = c(rep(480,3)), ystart = const2535[,2], yfin = const2535[,2])
const2535upper <- data.frame(xstart = c(rep(120,3)), xfin = c(rep(480,3)), ystart = const2535[,3], yfin = const2535[,3])

#32-42 °C
fluct3242 <- data.frame(mvpost_all[c(3,4),])
colnames(fluct3242) <- c("mean", "lower", "upper")
fluct3242$freq <- c(120, 480)

const3242 <- data.frame(mvpost_all[c(5,6,7),])
colnames(const3242) <- c("mean", "lower", "upper")

const3242mean <- data.frame(xstart = c(rep(120,3)), xfin = c(rep(480,3)), ystart = const3242[,1], yfin = const3242[,1])
const3242lower<- data.frame(xstart = c(rep(120,3)), xfin = c(rep(480,3)), ystart = const3242[,2], yfin = const3242[,2])
const3242upper <- data.frame(xstart = c(rep(120,3)), xfin = c(rep(480,3)), ystart = const3242[,3], yfin = const3242[,3])

#Combine data frames for fluctuations and add column for temperature treatment

fluctrange<- rbind(fluct2535, fluct3242)
fluctrange$Temp <- factor(c(rep("25-35", 2), rep("32-42", 2)))
fluctrange

const2535mean$Temp <- factor(c(rep("25-35", 3)))
const2535lower$Temp <- factor(c(rep("25-35", 3)))
const2535upper$Temp <- factor(c(rep("25-35", 3)))

const3242mean$Temp <- factor(c(rep("32-42", 3)))
const3242lower$Temp <- factor(c(rep("32-42", 3)))
const3242upper$Temp <- factor(c(rep("32-42", 3)))

#Plot estimates

labels <- as_labeller(c(`25-35` = "25-35 °C", `32-42` = "32-42 °C"))

facetrange <- ggplot() +
  geom_pointrange(data = fluctrange, aes(x = freq, y = mean, ymin = lower, ymax = upper), size=6, fatten =0.3) +
  geom_line(data = fluctrange, aes(x = freq, y = mean)) +
  geom_segment(data = const2535mean, aes(x = xstart, xend = xfin, y = ystart, yend = yfin), colour = c("#0072B2", "#009E73", "#D55E00")) + #Draw constant mean lines
  geom_segment(data = const2535lower, aes(x = xstart, xend = xfin, y = ystart, yend = yfin), colour = c("#0072B2", "#009E73", "#D55E00"), lty = "dotted") +
  geom_segment(data = const2535upper, aes(x = xstart, xend = xfin, y = ystart, yend = yfin), colour = c("#0072B2", "#009E73", "#D55E00"), lty = "dotted") +   
  geom_segment(data = const3242mean, aes(x = xstart, xend = xfin, y = ystart, yend = yfin), colour = c("#0072B2", "#009E73", "#D55E00")) + #Draw constant mean lines
  geom_segment(data = const3242lower, aes(x = xstart, xend = xfin, y = ystart, yend = yfin), colour = c("#0072B2", "#009E73", "#D55E00"), lty = "dotted") +
  geom_segment(data = const3242upper, aes(x = xstart, xend = xfin, y = ystart, yend = yfin), colour = c("#0072B2", "#009E73", "#D55E00"), lty = "dotted") +   
  ylab("Growth rate (mm/h)") +
  xlab("Frequency") +
  scale_x_continuous(breaks = c(120, 480), limits = c(60, 540), labels = c("Fast", "Slow"))  +
  scale_y_continuous(limits = c(0,4.5), breaks = seq(0,4.5,1))+
  theme_classic()+
  theme(strip.text.x = element_text(size =14))+
  theme(axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14))+
  facet_wrap(~ Temp, ncol=2, labeller = labels)
  
print(facetrange)

#Data frame for legend

legend_frame <- data.frame(x = rnorm(6), y = rnorm(6), type = c("Low", "Low", "Mean", "Mean", "High", "High"))

#Make plot to extract legend

legend_plot <- ggplot(legend_frame, aes(x = x, y = y, group = type, colour = type)) +
  geom_line() +
  scale_colour_manual(values = c("#D55E00", "#009E73", "#0072B2"), 
  name = "Constant (°C)", labels = c("High", "Mean", "Low"))+
  theme_classic()+
  theme(legend.position = "right")+
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=14))

legend_facet <- get_legend(legend_plot)

#combine with plot from raw data

finalfacet <- plot_grid(rawdata, facetrange, legend_facet, nrow=1, rel_widths = c(1, 1, 0.5), labels =c("A", "B"))

save_plot("Rawdata_modelestimates.png", finalfacet,  base_height = 4, base_width = 10)


#Data frame for quantitative genetics results

Rmatall <- matrix(rep(0, 10*10), ncol = 10)

colnames(Rmatall) <- c("25◦C", "30◦C", "32◦C", "35◦C", "37◦C", "42◦C", "25-35◦C Fast", "25-35◦C Slow", "32-42◦C Fast", "32-42◦C Slow")

rownames(Rmatall) <- c("25◦C", "30◦C", "32◦C", "35◦C", "37◦C", "42◦C", "25-35◦C Fast", "25-35◦C Slow", "32-42◦C Fast", "32-42◦C Slow")

Rmatall

#Pick genetic correlations

Rmatall[1,2] <- paste(mvpost_all[36,1], " (", mvpost_all[36,2], "-", mvpost_all[36,3], ")", sep = "")
Rmatall[1,3] <- paste(mvpost_all[26,1], " (", mvpost_all[26,2], "-", mvpost_all[26,3], ")", sep = "")
Rmatall[1,4] <- paste(mvpost_all[44,1], " (", mvpost_all[44,2], "-", mvpost_all[44,3], ")", sep = "")
Rmatall[1,5] <- paste(mvpost_all[27,1], " (", mvpost_all[27,2], "-", mvpost_all[27,3], ")", sep = "")
Rmatall[1,6] <- paste(mvpost_all[28,1], " (", mvpost_all[28,2], "-", mvpost_all[28,3], ")", sep = "")
Rmatall[1,7] <- paste(mvpost_all[22,1], " (", mvpost_all[22,2], "-", mvpost_all[22,3], ")", sep = "")
Rmatall[1,8] <- paste(mvpost_all[23,1], " (", mvpost_all[23,2], "-", mvpost_all[23,3], ")", sep = "")
Rmatall[1,9] <- paste(mvpost_all[24,1], " (", mvpost_all[24,2], "-", mvpost_all[24,3], ")", sep = "")
Rmatall[1,10] <- paste(mvpost_all[25,1], " (", mvpost_all[25,2], "-", mvpost_all[25,3], ")", sep = "")

Rmatall[2,3] <- paste(mvpost_all[33,1], " (", mvpost_all[33,2], "-", mvpost_all[33,3], ")", sep = "")
Rmatall[2,4] <- paste(mvpost_all[45,1], " (", mvpost_all[45,2], "-", mvpost_all[45,3], ")", sep = "")
Rmatall[2,5] <- paste(mvpost_all[34,1], " (", mvpost_all[34,2], "-", mvpost_all[34,3], ")", sep = "")
Rmatall[2,6] <- paste(mvpost_all[35,1], " (", mvpost_all[35,2], "-", mvpost_all[35,3], ")", sep = "")
Rmatall[2,7] <- paste(mvpost_all[29,1], " (", mvpost_all[29,2], "-", mvpost_all[29,3], ")", sep = "")
Rmatall[2,8] <- paste(mvpost_all[30,1], " (", mvpost_all[30,2], "-", mvpost_all[30,3], ")", sep = "")
Rmatall[2,9] <- paste(mvpost_all[31,1], " (", mvpost_all[31,2], "-", mvpost_all[31,3], ")", sep = "")
Rmatall[2,10] <- paste(mvpost_all[32,1], " (", mvpost_all[32,2], "-", mvpost_all[32,3], ")", sep = "")

Rmatall[3,4] <- paste(mvpost_all[41,1], " (", mvpost_all[41,2], "-", mvpost_all[41,3], ")", sep = "")
Rmatall[3,5] <- paste(mvpost_all[15,1], " (", mvpost_all[15,2], "-", mvpost_all[15,3], ")", sep = "")
Rmatall[3,6] <- paste(mvpost_all[20,1], " (", mvpost_all[20,2], "-", mvpost_all[20,3], ")", sep = "")
Rmatall[3,7] <- paste(mvpost_all[7,1], " (", mvpost_all[7,2], "-", mvpost_all[7,3], ")", sep = "")
Rmatall[3,8] <- paste(mvpost_all[8,1], " (", mvpost_all[8,2], "-", mvpost_all[8,3], ")", sep = "")
Rmatall[3,9] <- paste(mvpost_all[9,1], " (", mvpost_all[9,2], "-", mvpost_all[9,3], ")", sep = "")
Rmatall[3,10] <- paste(mvpost_all[10,1], " (", mvpost_all[10,2], "-", mvpost_all[10,3], ")", sep = "")

Rmatall[4,5] <- paste(mvpost_all[42,1], " (", mvpost_all[42,2], "-", mvpost_all[42,3], ")", sep = "")
Rmatall[4,6] <- paste(mvpost_all[43,1], " (", mvpost_all[43,2], "-", mvpost_all[43,3], ")", sep = "")
Rmatall[4,7] <- paste(mvpost_all[37,1], " (", mvpost_all[37,2], "-", mvpost_all[37,3], ")", sep = "")
Rmatall[4,8] <- paste(mvpost_all[38,1], " (", mvpost_all[38,2], "-", mvpost_all[38,3], ")", sep = "")
Rmatall[4,9] <- paste(mvpost_all[39,1], " (", mvpost_all[39,2], "-", mvpost_all[39,3], ")", sep = "")
Rmatall[4,10] <- paste(mvpost_all[40,1], " (", mvpost_all[40,2], "-", mvpost_all[40,3], ")", sep = "")

Rmatall[5,6] <- paste(mvpost_all[21,1], " (", mvpost_all[21,2], "-", mvpost_all[21,3], ")", sep = "")
Rmatall[5,7] <- paste(mvpost_all[11,1], " (", mvpost_all[11,2], "-", mvpost_all[11,3], ")", sep = "")
Rmatall[5,8] <- paste(mvpost_all[12,1], " (", mvpost_all[12,2], "-", mvpost_all[12,3], ")", sep = "")
Rmatall[5,9] <- paste(mvpost_all[13,1], " (", mvpost_all[13,2], "-", mvpost_all[13,3], ")", sep = "")
Rmatall[5,10] <- paste(mvpost_all[14,1], " (", mvpost_all[14,2], "-", mvpost_all[14,3], ")", sep = "")

Rmatall[6,7] <- paste(mvpost_all[16,1], " (", mvpost_all[16,2], "-", mvpost_all[16,3], ")", sep = "")
Rmatall[6,8] <- paste(mvpost_all[17,1], " (", mvpost_all[17,2], "-", mvpost_all[17,3], ")", sep = "")
Rmatall[6,9] <- paste(mvpost_all[18,1], " (", mvpost_all[18,2], "-", mvpost_all[18,3], ")", sep = "")
Rmatall[6,10] <- paste(mvpost_all[19,1], " (", mvpost_all[19,2], "-", mvpost_all[19,3], ")", sep = "")

Rmatall[7,8] <- paste(mvpost_all[1,1], " (", mvpost_all[1,2], "-", mvpost_all[1,3], ")", sep = "")
Rmatall[7,9] <- paste(mvpost_all[2,1], " (", mvpost_all[2,2], "-", mvpost_all[2,3], ")", sep = "")
Rmatall[7,10] <- paste(mvpost_all[4,1], " (", mvpost_all[4,2], "-", mvpost_all[4,3], ")", sep = "")

Rmatall[8,9] <- paste(mvpost_all[3,1], " (", mvpost_all[3,2], "-", mvpost_all[3,3], ")", sep = "")
Rmatall[8,10] <- paste(mvpost_all[5,1], " (", mvpost_all[5,2], "-", mvpost_all[5,3], ")", sep = "")

Rmatall[9,10] <- paste(mvpost_all[6,1], " (", mvpost_all[6,2], "-", mvpost_all[6,3], ")", sep = "")

Rmatall

#Count genetic covariances

#25 °C

cov25_30 <- (mvpost$sd_gr25)*(mvpost$sd_gr30)*(mvpost$cor_gr25_gr30)
covq25_30 <- quantile(cov25_30, probs = c(0.025, 0.975))
covq25_30 <- round(covq25_30, digits = 2) 
covq25_30 <- as.data.frame(covq25_30)
covm25_30 <- mean(cov25_30)
covm25_30 <- round(covm25_30, digits = 2)
Rmatall[2,1] <- paste(covm25_30, " (", covq25_30[1,1], "-", covq25_30[2,1], ")", sep = "")

cov25_32 <- (mvpost$sd_gr25)*(mvpost$sd_gr32)*(mvpost$cor_gr32_gr25)
covq25_32 <- quantile(cov25_32, probs = c(0.025, 0.975))
covq25_32 <- round(covq25_32, digits = 2) 
covq25_32 <- as.data.frame(covq25_32)
covm25_32 <- mean(cov25_32)
covm25_32 <- round(covm25_32, digits = 2)
Rmatall[3,1] <- paste(covm25_32, " (", covq25_32[1,1], "-", covq25_32[2,1], ")", sep = "")

cov25_35 <- (mvpost$sd_gr25)*(mvpost$sd_gr35)*(mvpost$cor_gr25_gr35)
covq25_35 <- quantile(cov25_35, probs = c(0.025, 0.975))
covq25_35 <- round(covq25_35, digits = 2) 
covq25_35 <- as.data.frame(covq25_35)
covm25_35 <- mean(cov25_35)
covm25_35 <- round(covm25_35, digits = 2)
Rmatall[4,1] <- paste(covm25_35, " (", covq25_35[1,1], "-", covq25_35[2,1], ")", sep = "")

cov25_37 <- (mvpost$sd_gr25)*(mvpost$sd_gr37)*(mvpost$cor_gr37_gr25)
covq25_37 <- quantile(cov25_37, probs = c(0.025, 0.975))
covq25_37 <- round(covq25_37, digits = 2) 
covq25_37 <- as.data.frame(covq25_37)
covm25_37 <- mean(cov25_37)
covm25_37 <- round(covm25_37, digits = 2)
Rmatall[5,1] <- paste(covm25_37, " (", covq25_37[1,1], "-", covq25_37[2,1], ")", sep = "")

cov25_42 <- (mvpost$sd_gr25)*(mvpost$sd_gr42)*(mvpost$cor_gr42_gr25)
covq25_42 <- quantile(cov25_42, probs = c(0.025, 0.975))
covq25_42 <- round(covq25_42, digits = 2) 
covq25_42 <- as.data.frame(covq25_42)
covm25_42 <- mean(cov25_42)
covm25_42 <- round(covm25_42, digits = 2)
Rmatall[6,1] <- paste(covm25_42, " (", covq25_42[1,1], "-", covq25_42[2,1], ")", sep = "")

cov25_2535F <- (mvpost$sd_gr25)*(mvpost$sd_gr2535F)*(mvpost$cor_gr2535F_gr25)
covq25_2535F <- quantile(cov25_2535F, probs = c(0.025, 0.975))
covq25_2535F <- round(covq25_2535F, digits = 2) 
covq25_2535F <- as.data.frame(covq25_2535F)
covm25_2535F <- mean(cov25_2535F)
covm25_2535F <- round(covm25_2535F, digits = 2)
Rmatall[7,1] <- paste(covm25_2535F, " (", covq25_2535F[1,1], "-", covq25_2535F[2,1], ")", sep = "")

cov25_2535S <- (mvpost$sd_gr25)*(mvpost$sd_gr2535S)*(mvpost$cor_gr2535S_gr25)
covq25_2535S <- quantile(cov25_2535S, probs = c(0.025, 0.975))
covq25_2535S <- round(covq25_2535S, digits = 2) 
covq25_2535S <- as.data.frame(covq25_2535S)
covm25_2535S <- mean(cov25_2535S)
covm25_2535S <- round(covm25_2535S, digits = 2)
Rmatall[8,1] <- paste(covm25_2535S, " (", covq25_2535S[1,1], "-", covq25_2535S[2,1], ")", sep = "")

cov25_3242F <- (mvpost$sd_gr25)*(mvpost$sd_gr3242F)*(mvpost$cor_gr3242F_gr25)
covq25_3242F <- quantile(cov25_3242F, probs = c(0.025, 0.975))
covq25_3242F <- round(covq25_3242F, digits = 2) 
covq25_3242F <- as.data.frame(covq25_3242F)
covm25_3242F <- mean(cov25_3242F)
covm25_3242F <- round(covm25_3242F, digits = 2)
Rmatall[9,1] <- paste(covm25_3242F, " (", covq25_3242F[1,1], "-", covq25_3242F[2,1], ")", sep = "")

cov25_3242S <- (mvpost$sd_gr25)*(mvpost$sd_gr3242S)*(mvpost$cor_gr3242S_gr25)
covq25_3242S <- quantile(cov25_3242S, probs = c(0.025, 0.975))
covq25_3242S <- round(covq25_3242S, digits = 2) 
covq25_3242S <- as.data.frame(covq25_3242S)
covm25_3242S <- mean(cov25_3242S)
covm25_3242S <- round(covm25_3242S, digits = 2)
Rmatall[10,1] <- paste(covm25_3242S, " (", covq25_3242S[1,1], "-", covq25_3242S[2,1], ")", sep = "")

#30 °C

cov30_32 <- (mvpost$sd_gr30)*(mvpost$sd_gr32)*(mvpost$cor_gr32_gr30)
covq30_32 <- quantile(cov30_32, probs = c(0.025, 0.975))
covq30_32 <- round(covq30_32, digits = 2) 
covq30_32 <- as.data.frame(covq30_32)
covm30_32 <- mean(cov30_32)
covm30_32 <- round(covm30_32, digits = 2)
Rmatall[3,2] <- paste(covm30_32, " (", covq30_32[1,1], "-", covq30_32[2,1], ")", sep = "")

cov30_35 <- (mvpost$sd_gr30)*(mvpost$sd_gr35)*(mvpost$cor_gr30_gr35)
covq30_35 <- quantile(cov30_35, probs = c(0.025, 0.975))
covq30_35 <- round(covq30_35, digits = 2) 
covq30_35 <- as.data.frame(covq30_35)
covm30_35 <- mean(cov30_35)
covm30_35 <- round(covm30_35, digits = 2)
Rmatall[4,2] <- paste(covm30_35, " (", covq30_35[1,1], "-", covq30_35[2,1], ")", sep = "")

cov30_37 <- (mvpost$sd_gr30)*(mvpost$sd_gr37)*(mvpost$cor_gr37_gr30)
covq30_37 <- quantile(cov30_37, probs = c(0.025, 0.975))
covq30_37 <- round(covq30_37, digits = 2) 
covq30_37 <- as.data.frame(covq30_37)
covm30_37 <- mean(cov30_37)
covm30_37 <- round(covm30_37, digits = 2)
Rmatall[5,2] <- paste(covm30_37, " (", covq30_37[1,1], "-", covq30_37[2,1], ")", sep = "")

cov30_42 <- (mvpost$sd_gr30)*(mvpost$sd_gr42)*(mvpost$cor_gr42_gr30)
covq30_42 <- quantile(cov30_42, probs = c(0.025, 0.975))
covq30_42 <- round(covq30_42, digits = 2) 
covq30_42 <- as.data.frame(covq30_42)
covm30_42 <- mean(cov30_42)
covm30_42 <- round(covm30_42, digits = 2)
Rmatall[6,2] <- paste(covm30_42, " (", covq30_42[1,1], "-", covq30_42[2,1], ")", sep = "")

cov30_2535F <- (mvpost$sd_gr30)*(mvpost$sd_gr2535F)*(mvpost$cor_gr2535F_gr30)
covq30_2535F <- quantile(cov30_2535F, probs = c(0.025, 0.975))
covq30_2535F <- round(covq30_2535F, digits = 2) 
covq30_2535F <- as.data.frame(covq30_2535F)
covm30_2535F <- mean(cov30_2535F)
covm30_2535F <- round(covm30_2535F, digits = 2)
Rmatall[7,2] <- paste(covm30_2535F, " (", covq30_2535F[1,1], "-", covq30_2535F[2,1], ")", sep = "")

cov30_2535S <- (mvpost$sd_gr30)*(mvpost$sd_gr2535S)*(mvpost$cor_gr2535S_gr30)
covq30_2535S <- quantile(cov30_2535S, probs = c(0.025, 0.975))
covq30_2535S <- round(covq30_2535S, digits = 2) 
covq30_2535S <- as.data.frame(covq30_2535S)
covm30_2535S <- mean(cov30_2535S)
covm30_2535S <- round(covm30_2535S, digits = 2)
Rmatall[8,2] <- paste(covm30_2535S, " (", covq30_2535S[1,1], "-", covq30_2535S[2,1], ")", sep = "")

cov30_3242F <- (mvpost$sd_gr30)*(mvpost$sd_gr3242F)*(mvpost$cor_gr3242F_gr30)
covq30_3242F <- quantile(cov30_3242F, probs = c(0.025, 0.975))
covq30_3242F <- round(covq30_3242F, digits = 2) 
covq30_3242F <- as.data.frame(covq30_3242F)
covm30_3242F <- mean(cov30_3242F)
covm30_3242F <- round(covm30_3242F, digits = 2)
Rmatall[9,2] <- paste(covm30_3242F, " (", covq30_3242F[1,1], "-", covq30_3242F[2,1], ")", sep = "")

cov30_3242S <- (mvpost$sd_gr30)*(mvpost$sd_gr3242S)*(mvpost$cor_gr3242S_gr30)
covq30_3242S <- quantile(cov30_3242S, probs = c(0.025, 0.975))
covq30_3242S <- round(covq30_3242S, digits = 2) 
covq30_3242S <- as.data.frame(covq30_3242S)
covm30_3242S <- mean(cov30_3242S)
covm30_3242S <- round(covm30_3242S, digits = 2)
Rmatall[10,2] <- paste(covm30_3242S, " (", covq30_3242S[1,1], "-", covq30_3242S[2,1], ")", sep = "")

#32 °C

cov32_35 <- (mvpost$sd_gr32)*(mvpost$sd_gr35)*(mvpost$cor_gr32_gr35)
covq32_35 <- quantile(cov32_35, probs = c(0.025, 0.975))
covq32_35 <- round(covq32_35, digits = 2) 
covq32_35 <- as.data.frame(covq32_35)
covm32_35 <- mean(cov32_35)
covm32_35 <- round(covm32_35, digits = 2)
Rmatall[4,3] <- paste(covm32_35, " (", covq32_35[1,1], "-", covq32_35[2,1], ")", sep = "")

cov32_37 <- (mvpost$sd_gr32)*(mvpost$sd_gr37)*(mvpost$cor_gr32_gr37)
covq32_37 <- quantile(cov32_37, probs = c(0.025, 0.975))
covq32_37 <- round(covq32_37, digits = 2) 
covq32_37 <- as.data.frame(covq32_37)
covm32_37 <- mean(cov32_37)
covm32_37 <- round(covm32_37, digits = 2)
Rmatall[5,3] <- paste(covm32_37, " (", covq32_37[1,1], "-", covq32_37[2,1], ")", sep = "")

cov32_42 <- (mvpost$sd_gr32)*(mvpost$sd_gr42)*(mvpost$cor_gr32_gr42)
covq32_42 <- quantile(cov32_42, probs = c(0.025, 0.975))
covq32_42 <- round(covq32_42, digits = 2) 
covq32_42 <- as.data.frame(covq32_42)
covm32_42 <- mean(cov32_42)
covm32_42 <- round(covm32_42, digits = 2)
Rmatall[6,3] <- paste(covm32_42, " (", covq32_42[1,1], "-", covq32_42[2,1], ")", sep = "")

cov32_2535F <- (mvpost$sd_gr32)*(mvpost$sd_gr2535F)*(mvpost$cor_gr2535F_gr32)
covq32_2535F <- quantile(cov32_2535F, probs = c(0.025, 0.975))
covq32_2535F <- round(covq32_2535F, digits = 2) 
covq32_2535F <- as.data.frame(covq32_2535F)
covm32_2535F <- mean(cov32_2535F)
covm32_2535F <- round(covm32_2535F, digits = 2)
Rmatall[7,3] <- paste(covm32_2535F, " (", covq32_2535F[1,1], "-", covq32_2535F[2,1], ")", sep = "")

cov32_2535S <- (mvpost$sd_gr32)*(mvpost$sd_gr2535S)*(mvpost$cor_gr2535S_gr32)
covq32_2535S <- quantile(cov32_2535S, probs = c(0.025, 0.975))
covq32_2535S <- round(covq32_2535S, digits = 2) 
covq32_2535S <- as.data.frame(covq32_2535S)
covm32_2535S <- mean(cov32_2535S)
covm32_2535S <- round(covm32_2535S, digits = 2)
Rmatall[8,3] <- paste(covm32_2535S, " (", covq32_2535S[1,1], "-", covq32_2535S[2,1], ")", sep = "")

cov32_3242F <- (mvpost$sd_gr32)*(mvpost$sd_gr3242F)*(mvpost$cor_gr3242F_gr32)
covq32_3242F <- quantile(cov32_3242F, probs = c(0.025, 0.975))
covq32_3242F <- round(covq32_3242F, digits = 2) 
covq32_3242F <- as.data.frame(covq32_3242F)
covm32_3242F <- mean(cov32_3242F)
covm32_3242F <- round(covm32_3242F, digits = 2)
Rmatall[9,3] <- paste(covm32_3242F, " (", covq32_3242F[1,1], "-", covq32_3242F[2,1], ")", sep = "")

cov32_3242S <- (mvpost$sd_gr32)*(mvpost$sd_gr3242S)*(mvpost$cor_gr3242S_gr32)
covq32_3242S <- quantile(cov32_3242S, probs = c(0.025, 0.975))
covq32_3242S <- round(covq32_3242S, digits = 2) 
covq32_3242S <- as.data.frame(covq32_3242S)
covm32_3242S <- mean(cov32_3242S)
covm32_3242S <- round(covm32_3242S, digits = 2)
Rmatall[10,3] <- paste(covm32_3242S, " (", covq32_3242S[1,1], "-", covq32_3242S[2,1], ")", sep = "")

#35 °C

cov35_37 <- (mvpost$sd_gr35)*(mvpost$sd_gr37)*(mvpost$cor_gr37_gr35)
covq35_37 <- quantile(cov35_37, probs = c(0.025, 0.975))
covq35_37 <- round(covq35_37, digits = 2) 
covq35_37 <- as.data.frame(covq35_37)
covm35_37 <- mean(cov35_37)
covm35_37 <- round(covm35_37, digits = 2)
Rmatall[5,4] <- paste(covm35_37, " (", covq35_37[1,1], "-", covq35_37[2,1], ")", sep = "")

cov35_42 <- (mvpost$sd_gr35)*(mvpost$sd_gr42)*(mvpost$cor_gr42_gr35)
covq35_42 <- quantile(cov35_42, probs = c(0.025, 0.975))
covq35_42 <- round(covq35_42, digits = 2) 
covq35_42 <- as.data.frame(covq35_42)
covm35_42 <- mean(cov35_42)
covm35_42 <- round(covm35_42, digits = 2)
Rmatall[6,4] <- paste(covm35_42, " (", covq35_42[1,1], "-", covq35_42[2,1], ")", sep = "")

cov35_2535F <- (mvpost$sd_gr35)*(mvpost$sd_gr2535F)*(mvpost$cor_gr2535F_gr35)
covq35_2535F <- quantile(cov35_2535F, probs = c(0.025, 0.975))
covq35_2535F <- round(covq35_2535F, digits = 2) 
covq35_2535F <- as.data.frame(covq35_2535F)
covm35_2535F <- mean(cov35_2535F)
covm35_2535F <- round(covm35_2535F, digits = 2)
Rmatall[7,4] <- paste(covm35_2535F, " (", covq35_2535F[1,1], "-", covq35_2535F[2,1], ")", sep = "")

cov35_2535S <- (mvpost$sd_gr35)*(mvpost$sd_gr2535S)*(mvpost$cor_gr2535S_gr35)
covq35_2535S <- quantile(cov35_2535S, probs = c(0.025, 0.975))
covq35_2535S <- round(covq35_2535S, digits = 2) 
covq35_2535S <- as.data.frame(covq35_2535S)
covm35_2535S <- mean(cov35_2535S)
covm35_2535S <- round(covm35_2535S, digits = 2)
Rmatall[8,4] <- paste(covm35_2535S, " (", covq35_2535S[1,1], "-", covq35_2535S[2,1], ")", sep = "")

cov35_3242F <- (mvpost$sd_gr35)*(mvpost$sd_gr3242F)*(mvpost$cor_gr3242F_gr35)
covq35_3242F <- quantile(cov35_3242F, probs = c(0.025, 0.975))
covq35_3242F <- round(covq35_3242F, digits = 2) 
covq35_3242F <- as.data.frame(covq35_3242F)
covm35_3242F <- mean(cov35_3242F)
covm35_3242F <- round(covm35_3242F, digits = 2)
Rmatall[9,4] <- paste(covm35_3242F, " (", covq35_3242F[1,1], "-", covq35_3242F[2,1], ")", sep = "")

cov35_3242S <- (mvpost$sd_gr35)*(mvpost$sd_gr3242S)*(mvpost$cor_gr3242S_gr35)
covq35_3242S <- quantile(cov35_3242S, probs = c(0.025, 0.975))
covq35_3242S <- round(covq35_3242S, digits = 2) 
covq35_3242S <- as.data.frame(covq35_3242S)
covm35_3242S <- mean(cov35_3242S)
covm35_3242S <- round(covm35_3242S, digits = 2)
Rmatall[10,4] <- paste(covm35_3242S, " (", covq35_3242S[1,1], "-", covq35_3242S[2,1], ")", sep = "")

#37 °C

cov37_42 <- (mvpost$sd_gr37)*(mvpost$sd_gr42)*(mvpost$cor_gr37_gr42)
covq37_42 <- quantile(cov37_42, probs = c(0.025, 0.975))
covq37_42 <- round(covq37_42, digits = 2) 
covq37_42 <- as.data.frame(covq37_42)
covm37_42 <- mean(cov37_42)
covm37_42 <- round(covm37_42, digits = 2)
Rmatall[6,5] <- paste(covm37_42, " (", covq37_42[1,1], "-", covq37_42[2,1], ")", sep = "")

cov37_2535F <- (mvpost$sd_gr37)*(mvpost$sd_gr2535F)*(mvpost$cor_gr2535F_gr37)
covq37_2535F <- quantile(cov37_2535F, probs = c(0.025, 0.975))
covq37_2535F <- round(covq37_2535F, digits = 2) 
covq37_2535F <- as.data.frame(covq37_2535F)
covm37_2535F <- mean(cov37_2535F)
covm37_2535F <- round(covm37_2535F, digits = 2)
Rmatall[7,5] <- paste(covm37_2535F, " (", covq37_2535F[1,1], "-", covq37_2535F[2,1], ")", sep = "")

cov37_2535S <- (mvpost$sd_gr37)*(mvpost$sd_gr2535S)*(mvpost$cor_gr2535S_gr37)
covq37_2535S <- quantile(cov37_2535S, probs = c(0.025, 0.975))
covq37_2535S <- round(covq37_2535S, digits = 2) 
covq37_2535S <- as.data.frame(covq37_2535S)
covm37_2535S <- mean(cov37_2535S)
covm37_2535S <- round(covm37_2535S, digits = 2)
Rmatall[8,5] <- paste(covm37_2535S, " (", covq37_2535S[1,1], "-", covq37_2535S[2,1], ")", sep = "")

cov37_3242F <- (mvpost$sd_gr37)*(mvpost$sd_gr3242F)*(mvpost$cor_gr3242F_gr37)
covq37_3242F <- quantile(cov37_3242F, probs = c(0.025, 0.975))
covq37_3242F <- round(covq37_3242F, digits = 2) 
covq37_3242F <- as.data.frame(covq37_3242F)
covm37_3242F <- mean(cov37_3242F)
covm37_3242F <- round(covm37_3242F, digits = 2)
Rmatall[9,5] <- paste(covm37_3242F, " (", covq37_3242F[1,1], "-", covq37_3242F[2,1], ")", sep = "")

cov37_3242S <- (mvpost$sd_gr37)*(mvpost$sd_gr3242S)*(mvpost$cor_gr3242S_gr37)
covq37_3242S <- quantile(cov37_3242S, probs = c(0.025, 0.975))
covq37_3242S <- round(covq37_3242S, digits = 2) 
covq37_3242S <- as.data.frame(covq37_3242S)
covm37_3242S <- mean(cov37_3242S)
covm37_3242S <- round(covm37_3242S, digits = 2)
Rmatall[10,5] <- paste(covm37_3242S, " (", covq37_3242S[1,1], "-", covq37_3242S[2,1], ")", sep = "")

#42 °C

cov42_2535F <- (mvpost$sd_gr42)*(mvpost$sd_gr2535F)*(mvpost$cor_gr2535F_gr42)
covq42_2535F <- quantile(cov42_2535F, probs = c(0.025, 0.975))
covq42_2535F <- round(covq42_2535F, digits = 2) 
covq42_2535F <- as.data.frame(covq42_2535F)
covm42_2535F <- mean(cov42_2535F)
covm42_2535F <- round(covm42_2535F, digits = 2)
Rmatall[7,6] <- paste(covm42_2535F, " (", covq42_2535F[1,1], "-", covq42_2535F[2,1], ")", sep = "")

cov42_2535S <- (mvpost$sd_gr42)*(mvpost$sd_gr2535S)*(mvpost$cor_gr2535S_gr42)
covq42_2535S <- quantile(cov42_2535S, probs = c(0.025, 0.975))
covq42_2535S <- round(covq42_2535S, digits = 2) 
covq42_2535S <- as.data.frame(covq42_2535S)
covm42_2535S <- mean(cov42_2535S)
covm42_2535S <- round(covm42_2535S, digits = 2)
Rmatall[8,6] <- paste(covm42_2535S, " (", covq42_2535S[1,1], "-", covq42_2535S[2,1], ")", sep = "")

cov42_3242F <- (mvpost$sd_gr42)*(mvpost$sd_gr3242F)*(mvpost$cor_gr3242F_gr42)
covq42_3242F <- quantile(cov42_3242F, probs = c(0.025, 0.975))
covq42_3242F <- round(covq42_3242F, digits = 2) 
covq42_3242F <- as.data.frame(covq42_3242F)
covm42_3242F <- mean(cov42_3242F)
covm42_3242F <- round(covm42_3242F, digits = 2)
Rmatall[9,6] <- paste(covm42_3242F, " (", covq42_3242F[1,1], "-", covq42_3242F[2,1], ")", sep = "")

cov42_3242S <- (mvpost$sd_gr42)*(mvpost$sd_gr3242S)*(mvpost$cor_gr3242S_gr42)
covq42_3242S <- quantile(cov42_3242S, probs = c(0.025, 0.975))
covq42_3242S <- round(covq42_3242S, digits = 2) 
covq42_3242S <- as.data.frame(covq42_3242S)
covm42_3242S <- mean(cov42_3242S)
covm42_3242S <- round(covm42_3242S, digits = 2)
Rmatall[10,6] <- paste(covm42_3242S, " (", covq42_3242S[1,1], "-", covq42_3242S[2,1], ")", sep = "")

#25-35 °C fast

cov2535F_2535S <- (mvpost$sd_gr2535F)*(mvpost$sd_gr2535S)*(mvpost$cor_gr2535F_gr2535S)
covq2535F_2535S <- quantile(cov2535F_2535S, probs = c(0.025, 0.975))
covq2535F_2535S <- round(covq2535F_2535S, digits = 2) 
covq2535F_2535S <- as.data.frame(covq2535F_2535S)
covm2535F_2535S <- mean(cov2535F_2535S)
covm2535F_2535S <- round(covm2535F_2535S, digits = 2)
Rmatall[8,7] <- paste(covm2535F_2535S, " (", covq2535F_2535S[1,1], "-", covq2535F_2535S[2,1], ")", sep = "")

cov2535F_3242F <- (mvpost$sd_gr2535F)*(mvpost$sd_gr3242F)*(mvpost$cor_gr2535F_gr3242F)
covq2535F_3242F <- quantile(cov2535F_3242F, probs = c(0.025, 0.975))
covq2535F_3242F <- round(covq2535F_3242F, digits = 2) 
covq2535F_3242F <- as.data.frame(covq2535F_3242F)
covm2535F_3242F <- mean(cov2535F_3242F)
covm2535F_3242F <- round(covm2535F_3242F, digits = 2)
Rmatall[9,7] <- paste(covm2535F_3242F, " (", covq2535F_3242F[1,1], "-", covq2535F_3242F[2,1], ")", sep = "")

cov2535F_3242S <- (mvpost$sd_gr2535F)*(mvpost$sd_gr3242S)*(mvpost$cor_gr2535F_gr3242S)
covq2535F_3242S <- quantile(cov2535F_3242S, probs = c(0.025, 0.975))
covq2535F_3242S <- round(covq2535F_3242S, digits = 2) 
covq2535F_3242S <- as.data.frame(covq2535F_3242S)
covm2535F_3242S <- mean(cov2535F_3242S)
covm2535F_3242S <- round(covm2535F_3242S, digits = 2)
Rmatall[10,7] <- paste(covm2535F_3242S, " (", covq2535F_3242S[1,1], "-", covq2535F_3242S[2,1], ")", sep = "")

#25-35 °C slow

cov2535S_3242F <- (mvpost$sd_gr2535S)*(mvpost$sd_gr3242F)*(mvpost$cor_gr2535S_gr3242F)
covq2535S_3242F <- quantile(cov2535S_3242F, probs = c(0.025, 0.975))
covq2535S_3242F <- round(covq2535S_3242F, digits = 2) 
covq2535S_3242F <- as.data.frame(covq2535S_3242F)
covm2535S_3242F <- mean(cov2535S_3242F)
covm2535S_3242F <- round(covm2535S_3242F, digits = 2)
Rmatall[9,8] <- paste(covm2535S_3242F, " (", covq2535S_3242F[1,1], "-", covq2535S_3242F[2,1], ")", sep = "")

cov2535S_3242S <- (mvpost$sd_gr2535S)*(mvpost$sd_gr3242S)*(mvpost$cor_gr2535S_gr3242S)
covq2535S_3242S <- quantile(cov2535S_3242S, probs = c(0.025, 0.975))
covq2535S_3242S <- round(covq2535S_3242S, digits = 2) 
covq2535S_3242S <- as.data.frame(covq2535S_3242S)
covm2535S_3242S <- mean(cov2535S_3242S)
covm2535S_3242S <- round(covm2535S_3242S, digits = 2)
Rmatall[10,8] <- paste(covm2535S_3242S, " (", covq2535S_3242S[1,1], "-", covq2535S_3242S[2,1], ")", sep = "")

#32-42 °C fast

cov3242F_3242S <- (mvpost$sd_gr3242F)*(mvpost$sd_gr3242S)*(mvpost$cor_gr3242F_gr3242S)
covq3242F_3242S <- quantile(cov3242F_3242S, probs = c(0.025, 0.975))
covq3242F_3242S <- round(covq3242F_3242S, digits = 2) 
covq3242F_3242S <- as.data.frame(covq3242F_3242S)
covm3242F_3242S <- mean(cov3242F_3242S)
covm3242F_3242S <- round(covm3242F_3242S, digits = 2)
Rmatall[10,9] <- paste(covm3242F_3242S, " (", covq3242F_3242S[1,1], "-", covq3242F_3242S[2,1], ")", sep = "")

Rmatall

#Pick genetic variances on the diagonal

Var_all <- mvpost_all[c(11:20, 111:120), 1:3]^2

Var_all <- round(Var_all, digits = 2)

Rmatall[1,1] <- paste(Var_all[8,1], " (", Var_all[8,2], "-", Var_all[8,3], ")", sep = "")
Rmatall[2,2] <- paste(Var_all[9,1], " (", Var_all[9,2], "-", Var_all[9,3], ")", sep = "")
Rmatall[3,3] <- paste(Var_all[5,1], " (", Var_all[5,2], "-", Var_all[5,3], ")", sep = "")
Rmatall[4,4] <- paste(Var_all[10,1], " (", Var_all[10,2], "-", Var_all[10,3], ")", sep = "")
Rmatall[5,5] <- paste(Var_all[6,1], " (", Var_all[6,2], "-", Var_all[6,3], ")", sep = "")
Rmatall[6,6] <- paste(Var_all[7,1], " (", Var_all[7,2], "-", Var_all[7,3], ")", sep = "")
Rmatall[7,7] <- paste(Var_all[1,1], " (", Var_all[1,2], "-", Var_all[1,3], ")", sep = "")
Rmatall[8,8] <- paste(Var_all[2,1], " (", Var_all[2,2], "-", Var_all[2,3], ")", sep = "")
Rmatall[9,9] <- paste(Var_all[3,1], " (", Var_all[3,2], "-", Var_all[3,3], ")", sep = "")
Rmatall[10,10] <- paste(Var_all[4,1], " (", Var_all[4,2], "-", Var_all[4,3], ")", sep = "")

Rmatall

#Add column and pick environmental variances

Rmatall <- cbind(Rmatall, VarE=NA)

Rmatall[1,11] <- paste(Var_all[18,1], " (", Var_all[18,2], "-", Var_all[18,3], ")", sep = "")
Rmatall[2,11] <- paste(Var_all[19,1], " (", Var_all[19,2], "-", Var_all[19,3], ")", sep = "")
Rmatall[3,11] <- paste(Var_all[15,1], " (", Var_all[15,2], "-", Var_all[15,3], ")", sep = "")
Rmatall[4,11] <- paste(Var_all[20,1], " (", Var_all[20,2], "-", Var_all[20,3], ")", sep = "")
Rmatall[5,11] <- paste(Var_all[16,1], " (", Var_all[16,2], "-", Var_all[16,3], ")", sep = "")
Rmatall[6,11] <- paste(Var_all[17,1], " (", Var_all[17,2], "-", Var_all[17,3], ")", sep = "")
Rmatall[7,11] <- paste(Var_all[11,1], " (", Var_all[11,2], "-", Var_all[11,3], ")", sep = "")
Rmatall[8,11] <- paste(Var_all[12,1], " (", Var_all[12,2], "-", Var_all[12,3], ")", sep = "")
Rmatall[9,11] <- paste(Var_all[13,1], " (", Var_all[13,2], "-", Var_all[13,3], ")", sep = "")
Rmatall[10,11] <- paste(Var_all[14,1], " (", Var_all[14,2], "-", Var_all[14,3], ")", sep = "")

Rmatall

write.csv(Rmatall, file= "Correlation_matrix.csv")


#Principal component analysis to mean growth rate data

fluct_PCA <- read.csv("natpop_fluct_means.csv", header = T, dec = ".", sep = ",")

#Combine temperature and frequency into one treatment factor

fluct_PCA$Treatment <- paste(fluct_PCA$Temp, fluct_PCA$Frequency)

#Remove previous factors

fluct_PCA <- fluct_PCA[, -c(2,3)]

#Invert data frame to multivariate means per genotype

datawide <- spread(fluct_PCA, Treatment, meangr)

#Remove missing values

datawide_check <- complete.cases(datawide)

datawide <- datawide[datawide_check, ]

#Use princomp function for running PCA

PCA <- princomp(datawide[,-1])

summary(PCA)

loadings(PCA)

#Predicted values to be used in GWAS

PCA_predicted <- predict(PCA)

#Only the first three components

PCA_less <- PCA_predicted[,c(1:3)]

PCA_less <- data.frame(PCA_less)

#Put genotype factor back into data frame in first column

PCA_less <- cbind(PCA_less, Genot=datawide$Genot)

PCA_final <- PCA_less[, c(4,1,2,3)]

#Choose components individually for phenotype data

PCA1 <- PCA_final[,c(1,2)]
PCA2 <- PCA_final[,c(1,3)]
PCA3 <- PCA_final[,c(1,4)]

#save these as individual files for GWAS

write.csv(PCA1, file = "asmap_PCA1.csv", quote = FALSE, row.names = FALSE)
write.csv(PCA2, file = "asmap_PCA2.csv", quote = FALSE, row.names = FALSE)
write.csv(PCA3, file = "asmap_PCA3.csv", quote = FALSE, row.names = FALSE)

#Plot proportions of variance explained by the different principal components and
#the loadings of PC1, PC2 and PC3 for each temperature treatment as posterior means 
#and 95 % HPDI. 

#load posteriors of multivariate model
load("multivariate_model.RData")

#Extract posterior samples for intercepts of each temperature and each genotype
koe <- posterior_samples(model1, pars = c("^b", "^r_Genot"))

#Get temperature intercepts
int.gr2535F <- koe[,1]
int.gr2535S <- koe[,2]
int.gr3242F <- koe[,3]
int.gr3242S <- koe[,4]
int.gr32 <- koe[,5]
int.gr37 <- koe[,6]
int.gr42 <- koe[,7]
int.gr25 <- koe[,8]
int.gr30 <- koe[,9]
int.gr35 <- koe[,10]
#There are 435 genotype levels
genot.gr2535F <- koe[,11:445]
genot.gr2535S <- koe[,446:880]
genot.gr3242F <- koe[,881:1315]
genot.gr3242S <- koe[,1316:1750]
genot.gr32 <- koe[,(10+1+435*4):(435*5+10)]
genot.gr37 <- koe[,(10+1+435*5):(435*6+10)]
genot.gr42 <- koe[,(10+1+435*6):(435*7+10)]
genot.gr25 <- koe[,(10+1+435*7):(435*8+10)]
genot.gr30 <- koe[,(10+1+435*8):(435*9+10)]
genot.gr35 <- koe[,(10+1+435*9):(435*10+10)]

#Combine temperature intercepts and genotypic effects
genot.gr2535F <- genot.gr2535F + int.gr2535F
genot.gr2535S <- genot.gr2535S + int.gr2535S
genot.gr3242F <- genot.gr3242F + int.gr3242F
genot.gr3242S <- genot.gr3242S + int.gr3242S
genot.gr32 <- genot.gr32 + int.gr32
genot.gr37 <- genot.gr37 + int.gr37
genot.gr42 <- genot.gr42 + int.gr42
genot.gr25 <- genot.gr25 + int.gr25
genot.gr30 <- genot.gr30 + int.gr30
genot.gr35 <- genot.gr35 + int.gr35

#Make a list of phenotypic data from the posteriors
foo <- list(0)
PCAdata <- rep(foo, 5000)

#Make datasets from the posterior
for(i in 1:5000) {
  PCAdata[[i]] <- cbind(t(genot.gr25[i,]), t(genot.gr2535F[i,]), t(genot.gr2535S[i,]), t(genot.gr30[i,]), t(genot.gr32[i,]), t(genot.gr3242F[i,]), t(genot.gr3242S[i,]), t(genot.gr35[i,]), t(genot.gr37[i,]), t(genot.gr42[i,]))
  colnames(PCAdata[[i]]) <- c("25 0", "25-35 120", "25-35 480", "30 0", "32 0", "32-42 120", "32-42 480", "35 0", "37 0", "42 0")
}

#Save the PCA data
save(PCAdata, file = "posterior_samples_for_PCA.RData")
#Load as necessary
#load("posterior_samples_for_PCA.RData")

#Make a dataframe for posterior distribution of PCA results
#Proportion of variance
post.PCA <- data.frame(matrix(rep(0,10*5000), ncol = 10))
colnames(post.PCA) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

#Loadings data needs to be a list of the 10 components
#Then the 10 phenotypes are the columns and 5000 rows the posterior samples
loadf <- data.frame(matrix(rep(0,10*5000), ncol = 10))
colnames(loadf) <- c("25 0", "25-35 120", "25-35 480", "30 0", "32 0", "32-42 120", "32-42 480", "35 0", "37 0", "42 0")
post.loadings <- rep(list(loadf),10)

#Do the PCAs
for(i in 1:5000) {
  PCApost <- princomp(PCAdata[[i]]) #Do the PCA
  post.PCA[i,] <- (PCApost$sdev^2)/sum(PCApost$sdev^2) #Calculate proportion of variance from sd's
  
  ### #Then store the loadings ####
  for(j in 1:10) {
    post.loadings[[j]][i,] <- PCApost$loadings[,j]
  }
}

#Calculate posterior means and HPDIs for the proportions of variance explained 

resmat.propvar.mean <- data.frame(proportion = rep(0,10), lower = rep(0, 10), upper = rep(0,10), Comp = paste(rep("PC", 10), 1:10, sep = ""))
resmat.propvar.mean[,1] <- apply(post.PCA, 2, mean) #Mean for each 
for(i in 1:10) { resmat.propvar.mean[i,2:3] <- HPDinterval(post.PCA[,i]) } #Get HPD intervals

#Calculate posterior means and HPDIs for loadings

resmat.loadings.mean <- data.frame(loading = rep(0,10*10), lower = rep(0,10*10), upper = rep(0,10*10), Temperature = factor(rep(c("25", "30", "35", "25-35 Fast", "25-35 Slow", "32", "37", "42", "32-42 Fast", "32-42 Slow"), 10), levels = c("25", "30", "35", "25-35 Fast", "25-35 Slow", "32", "37", "42", "32-42 Fast", "32-42 Slow")), pc = rep(c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), each = 10) )
for(i in 1:10) {
  resmat.loadings.mean[((i-1)*10+1):(i*10),1] <- apply(post.loadings[[i]], 2, mean)[c(1,4,8,2,3,5,9,10,6,7)] #Results in correct order
  resmat.loadings.mean[((i-1)*10+1):(i*10),2:3] <- HPDinterval(as.mcmc(post.loadings[[i]]))[c(1,4,8,2,3,5,9,10,6,7),]
}

#Save the posterior data
resmat.propvar.mean <- resmat.propvar.mean[1:6,]
resmat.loadings.mean <- resmat.loadings.mean[1:30,]

save(resmat.propvar.mean, resmat.loadings.mean, file = "HPDI_for_PCA.RData")

#Load as necessary
#load("HPDI_for_PCA.RData")

#Draw figures based on PCA made on posterior distribution

prop.mean <- ggplot(resmat.propvar.mean, aes(x=Comp, y=proportion, ymin = lower, ymax = upper)) + 
  geom_bar(stat="identity", colour = "black", fill = "grey") +
  geom_errorbar(width = 0.2) +
  ylab("Proportion of variance explained") +
  xlab(" ") +
  theme_classic()+
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))+
  scale_y_continuous(limits=c(0.0,1.0), breaks=seq(0.0,1.0,0.1), expand = c(0,0))

loadings.mean <- ggplot(data = resmat.loadings.mean, aes(x = Temperature, y = loading, ymin = lower, ymax = upper)) +
  geom_pointrange(size = 0.2) +
  ylab("Loadings") +
  xlab("Temperature (°C)") +
  facet_grid(pc ~ .) +
  theme_classic()+
  geom_hline(yintercept = 0, lty = "dashed") +
  theme(strip.text.y = element_text(size =16))+
  theme(axis.text.x = element_text(angle = 30, vjust = 0.6, size = 14),
        axis.title.x = element_text(size = 16)) +
  theme(axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16))

fullpc <- plot_grid(prop.mean, loadings.mean, axis = 'b', align = "h", labels = c("A", "B"), nrow=1, ncol=2)

save_plot("PC_prop_loadings.png", fullpc, base_height = 6, base_width = 17)

#Running GWAS with GAPIT3 on CSC supercomputer https://www.puhti.csc.fi/public/

#Example code for running BLINK for growth rates at 25-35 °C fast fluctuations with 4 PCs

#Load genotype data

myG <- read.table(file = "/example_path/Neuro_hapmap_natpop_all_filtered.txt", header = FALSE, stringsAsFactors = FALSE)
geno.names <- myG[1,-c(1:11)]

#Load phenotypic data

genomeans <- read.csv(file = "/example_path/natpop_fluct_means.csv", header = T, sep = ",")

gr2535F <- dplyr::filter(genomeans, Temp == "25-35", Frequency == "120")

myY <- gr2535F[gr2535F$Genot %in% geno.names,]
myY <- myY[,-c(2,3)] #Drop the temp column, not needed
colnames(myY) <- c("Taxa", "gr2535F")

#Run GAPIT

myGAPIT2535F <- GAPIT(Y = myY, PCA.total = 4, G = myG, model = "Blink")

#Data for all associations was filtered with 0.01 treshold and manhattan plots 
#were drawn with custom scripts

#Load data files

assodata <- "C:/example_path/"

gwasdata.gr25 <- read.csv(paste(assodata, "25/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr25.csv", sep = ""))

gwasdata.gr30 <- read.csv(paste(assodata, "30/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr30.csv", sep = ""))

gwasdata.gr32 <- read.csv(paste(assodata, "32/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr32.csv", sep = ""))

gwasdata.gr35 <- read.csv(paste(assodata, "35/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr35.csv", sep = ""))

gwasdata.gr37 <- read.csv(paste(assodata, "37/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr37.csv", sep = ""))

gwasdata.gr42 <- read.csv(paste(assodata, "42/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr42.csv", sep = ""))

gwasdata.gr2535F <- read.csv(paste(assodata, "2535F/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr2535F.csv", sep = ""))

gwasdata.gr2535S <- read.csv(paste(assodata, "2535S/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr2535S.csv", sep = ""))

gwasdata.gr3242F <- read.csv(paste(assodata, "3242F/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr3242F.csv", sep = ""))

gwasdata.gr3242S <- read.csv(paste(assodata, "3242S/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.gr3242S.csv", sep = ""))

gwasdata.PCA1 <- read.csv(paste(assodata, "PCA1/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.PCA1.csv", sep = ""))

gwasdata.PCA2 <- read.csv(paste(assodata, "PCA2/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.PCA2.csv", sep = ""))

gwasdata.PCA3 <- read.csv(paste(assodata, "PCA3/new genotypes/filtered/GAPIT.Association.GWAS_Results.BLINK.PCA3.csv", sep = ""))

#This function makes a manhattan plot using ggplot2 commands

manhattan.plot <- function(gwasdata, signift, suggest, mylabel = NULL, hl = NULL) {
  ##Making cumulative base pair position
  ##This assumes that chromosome is the second column
  nCHR <- length(unique(gwasdata$Chr))
  gwasdata$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in 1:nCHR) {
    nbp[i] <- max(gwasdata[gwasdata$Chr == i,]$Pos)
    gwasdata[gwasdata$Chr == i,"BPcum"] <- gwasdata[gwasdata$Chr == i,"Pos"] + s
    s <- s + nbp[i]
  }
  
  axis.set <- group_by(gwasdata, Chr)
  axis.set <- summarize(axis.set, center = (max(BPcum) + min(BPcum))/2)
  ylim <- abs(floor(log10(min(gwasdata$P.value)))) + 2
  ylabel <- expression(paste(-"log"[10], "(p)", sep = ""))
  hl.bpcum <- NULL
  if(is.null(hl) == F) {
    #Lines to highlight candidates
    hl.ind <- match(hl, gwasdata[,1])
    hl.bpcum <- gwasdata[hl.ind, "BPcum"] }
  
  ggplot(gwasdata, aes(x = BPcum, y = -log10(P.value), 
                       color = as.factor(Chr))) +
    geom_vline(xintercept = hl.bpcum, color = "grey40", linetype = "dotted") +
    geom_point(alpha = 0.75, size = 2.5) +
    geom_hline(yintercept = -log10(signift), color = "grey40") +
    scale_x_continuous(label = axis.set$Chr, breaks = axis.set$center, expand = c(0.01,0)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 32.5), breaks = c(seq(0,32,4))) +
    scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
    labs(x = NULL, y = ylabel, subtitle = mylabel) +
    theme_classic() +
    theme( 
      legend.position = "none"
    )
  }

#This function calculates the bonferroni significance threshold

bonf.threshold <- function(pvalues, threshold) {
  bonft <- threshold / length(pvalues)
  return(bonft)
}

#Significance threshold corresponding p-value of 6.78×10−9
#This is the same for associations at all temperatures, because of the number of SNPs is the same

bonft <- bonf.threshold(asso.gr35.BLINK[,4], 0.01)

#Filtering significant SNPs for results section

signif.SNPs.gr25 <- filter(gwasdata.gr25, P.value < bonft)

signif.SNPs.gr32 <- filter(gwasdata.gr32, P.value < bonft)

signif.SNPs.gr30 <- filter(gwasdata.gr30, P.value < bonft)

signif.SNPs.gr35 <- filter(gwasdata.gr35, P.value < bonft)

signif.SNPs.gr37 <- filter(gwasdata.gr37, P.value < bonft)

signif.SNPs.gr42 <- filter(gwasdata.gr42, P.value < bonft)

signif.SNPs.gr2535F <- filter(gwasdata.gr2535F, P.value < bonft)

signif.SNPs.gr2535S <- filter(gwasdata.gr2535S, P.value < bonft)

signif.SNPs.gr3242F <- filter(gwasdata.gr3242F, P.value < bonft)

signif.SNPs.gr3242S <- filter(gwasdata.gr3242S, P.value < bonft)

signif.SNPs.PCA1 <- filter(gwasdata.PCA1, P.value < bonft)

signif.SNPs.PCA2 <- filter(gwasdata.PCA2, P.value < bonft)

signif.SNPs.PCA3 <- filter(gwasdata.PCA3, P.value < bonft)

#dotted vertical lines to highlight the candidate SNPs that were significant in
#at least one of the temperature treatments within the range

#25-35 °C

hl.gr25 <- signif.SNPs.gr25[,1]
hl.gr30 <- signif.SNPs.gr30[,1]
hl.gr35 <- signif.SNPs.gr35[,1]
hl.gr2535F <- signif.SNPs.gr2535F[,1]
hl.gr2535S <- signif.SNPs.gr2535S[,1]

hl.2535 <- c(as.character(hl.gr25), as.character(hl.gr30), as.character(hl.gr35), as.character(hl.gr2535F), as.character(hl.gr2535S))

#32-42 °C

hl.gr32 <- signif.SNPs.gr32[,1]
hl.gr37 <- signif.SNPs.gr37[,1]
hl.gr42 <- signif.SNPs.gr42[,1]
hl.gr3242F <- signif.SNPs.gr3242F[,1]
hl.gr3242S <- signif.SNPs.gr3242S[,1]

hl.3242 <- c(as.character(hl.gr32), as.character(hl.gr37), as.character(hl.gr42), as.character(hl.gr3242F), as.character(hl.gr3242S))

#All PCs

hl.PCA1 <- signif.SNPs.PCA1[,1]
hl.PCA2 <- signif.SNPs.PCA2[,1]
hl.PCA3 <- signif.SNPs.PCA3[,1]

hl.PCA.all <- c(as.character(hl.PCA1), as.character(hl.PCA2), as.character(hl.PCA3) )

#Labels to plots

label25 <- expression(paste("Growth rate at 25 ", degree, "C"))
label30 <- expression(paste("Growth rate at 30 ", degree, "C"))
label32 <- expression(paste("Growth rate at 32 ", degree, "C"))
label35 <- expression(paste("Growth rate at 35 ", degree, "C"))
label37 <- expression(paste("Growth rate at 37 ", degree, "C"))
label42 <- expression(paste("Growth rate at 42 ", degree, "C"))
label2535F <- expression(paste("Growth rate at 25-35 ", degree, "C ", "Fast"))
label2535S <- expression(paste("Growth rate at 25-35 ", degree, "C ", "Slow"))
label3242F <- expression(paste("Growth rate at 32-42 ", degree, "C ", "Fast"))
label3242S <- expression(paste("Growth rate at 32-42 ", degree, "C ", "Slow"))
labelPCA1 <- expression(paste("Growth rate for PC1"))
labelPCA2 <- expression(paste("Growth rate for PC2"))
labelPCA3 <- expression(paste("Growth rate for PC3"))

#Individual manhattan plots

gr25.plot <- manhattan.plot(gwasdata.gr25, signift = bonft, mylabel = label25, hl = hl.2535)

gr30.plot <- manhattan.plot(gwasdata.gr30, signift = bonft, mylabel = label30, hl = hl.2535)

gr32.plot <- manhattan.plot(gwasdata.gr32, signift = bonft, mylabel = label32, hl = hl.3242)

gr35.plot <- manhattan.plot(gwasdata.gr35, signift = bonft, mylabel = label35, hl = hl.2535)

gr37.plot <- manhattan.plot(gwasdata.gr37, signift = bonft, mylabel = label37, hl = hl.3242)

gr42.plot <- manhattan.plot(gwasdata.gr42, signift = bonft, mylabel = label42, hl = hl.3242)

gr2535F.plot <- manhattan.plot(gwasdata.gr2535F, signift = bonft, mylabel = label2535F, hl = hl.2535)

gr2535S.plot <- manhattan.plot(gwasdata.gr2535S, signift = bonft, mylabel = label2535S, hl = hl.2535)

gr3242F.plot <- manhattan.plot(gwasdata.gr3242F, signift = bonft, mylabel = label3242F, hl = hl.3242)

gr3242S.plot <- manhattan.plot(gwasdata.gr3242S, signift = bonft, mylabel = label3242S, hl = hl.3242)

PCA1.plot <- manhattan.plot(gwasdata.PCA1, signift = bonft, mylabel = labelPCA1, hl = hl.PCA.all)

PCA2.plot <- manhattan.plot(gwasdata.PCA2, signift = bonft, mylabel = labelPCA2, hl = hl.PCA.all)

PCA3.plot <- manhattan.plot(gwasdata.PCA3, signift = bonft, mylabel = labelPCA3, hl = hl.PCA.all)

#Combine plots for ranges and PCs

#25-35 °C

manplot.alltraits.hl.2535 <- plot_grid(gr25.plot, gr30.plot, gr35.plot, gr2535F.plot, gr2535S.plot, align = "v", ncol = 1)

save_plot("./Fluct_manhattan_2535.png", manplot.alltraits.hl.2535, base_height = 16, base_width = 11)

#32-42 °C

manplot.alltraits.hl.3242 <- plot_grid(gr32.plot, gr37.plot, gr42.plot, gr3242F.plot, gr3242S.plot, align = "v", ncol = 1)

save_plot("./Fluct_manhattan_3242.png", manplot.alltraits.hl.3242, base_height = 16, base_width = 11)

#PCs
manplot.PCA <- plot_grid(PCA1.plot, PCA2.plot, PCA3.plot, align = "v", ncol = 1)

save_plot("./PCA_manhattan.png", manplot.PCA, base_height = 10, base_width = 11)
