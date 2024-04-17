##############################################################
# title: "Alpha diversity in R - qiime2 output"
# author: "ANSC595"
# date: "March 16, 2021"
##############################################################

#The first step is very important. You need to set your working 
#directory. Just like in unix we have to `cd` to where our data is, the 
#same thing is true for R.
##############################################

#So where `~/Desktop/ANSC595/moving-pictures` is in the code below, you 
#need to enter the path to where you saved the tutorial or your data.
getwd()
setwd("/Users/nooshanmir/Documents/Purdue Courses/Spring 2024 (2nd Semester)/ASNC 516/Microbiome Final Project/Microbiome_Data_From_Annabel/FinalFiles")

list.files()



# core-metrics-results/evenness_vector.qza (alpha diversity)
# core-metrics-results/faith_pd_vector.qza (alpha diversity)
# core-metrics-results/observed_features_vector.qza (alpha diversity)
# core-metrics-results/shannon_vector.qza (alpha diversity)
# sample-metadata.tsv (metadata)


# Data manipulation
## Load Packages

library(tidyverse)
library(qiime2R)
library(ggpubr)

##Load Data
# In the code, the text before = is what the file will be called in R. 
# Make this short but unique as this is how you will tell R to use this 
# file in later commands.

# header: tells R that the first row is column names, not data
# row.names: tells R that the first column is row names, not data
# sep: tells R that the data are tab-delimited. 
# If you had a comma-delimited file, you would us sep=","

# Load data

meta<-read_q2metadata("Metadata.ACNAC.tsv")
str(meta)



evenness = read_qza("core-metrics-results/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

observed_features = read_qza("core-metrics-results/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

shannon = read_qza("core-metrics-results/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged

faith_pd = read_qza("core-metrics-results/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID") # this moves the sample names to a new column that matches the metadata and allows them to be merged\
faith_pd=faith_pd[,-1]
colnames(faith_pd)[colnames(faith_pd) == "V1"] <- "SampleID"
colnames(faith_pd)[colnames(faith_pd) == "V2"] <- "faith_pd"

## Clean up the data
# You can look at your data by clicking on it in the upper-right 
# quadrant "Environment"

# You always need to check the data types in your tables to make 
# sure they are what you want. We will now change some data types 
# in the meta now

str(meta)
#observed_features$observed_features_num <- lapply(observed_features$observed_features, as.numeric)
#observed_features$observed_features <- as.numeric(observed_features$observed_features)
str(observed_features)





###Alpha Diversity tables
# These tables will be merged for convenience and added to the 
# metadata table as the original tutorial was organized.

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
meta = merge(meta, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(meta) <- meta$SampleID
#meta = meta[,-7]

str(meta)


#Alpha-diversity
# Alpha-diversity is within sample diversity. It is how many 
# different species (OTUs) are in each sample (richness) and how 
# evenly they are distributed (evenness), which together are diversity. 
# Each sample has one value for each metric.


##Explore alpha metrics
# Now we will start to look at our data. We will first start with 
# alpha-diversity and richness. 
#
# You want the data to be roughly normal so that you can run ANOVA 
# or t-tests. If it is not normally distributed, you will need to 
# consider if you should normalize the data or usenon-parametric 
# tests such as Kruskal-Wallis.

# Here, we see that none of the data are normally distributed, 
# with the exception of "Faith" and "Observed Features".


#Plots
hist(meta$shannon_entropy, main="Shannon diversity", xlab="", breaks=10)
hist(meta$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(meta$pielou_eve, main="Evenness", xlab="", breaks=10)
hist(as.numeric(meta$observed_features), main="Observed Features", xlab="", breaks=10)

#Plots the qq-plot for residuals
ggqqplot(meta$shannon_entropy, title = "Shannon")
ggqqplot(meta$faith_pd, title = "Faith PD")
ggqqplot(meta$pielou_e, title = "Evenness")
ggqqplot(meta$observed_features, title = "Observed Features")







install.packages("ggpubr")
library("ggpubr")

# To test for normalcy statistically, we can run the Shapiro-Wilk 
# test of normality.

shapiro.test(meta$shannon)
shapiro.test(meta$faith_pd)
shapiro.test(meta$pielou_e)
shapiro.test(meta$observed_features)

# The null hypothesis of these tests is that “sample distribution 
# is normal”. If the test is significant(less tha 0.05), the distribution is non-normal. 

# We see that, as expected from the graphs, shannon and observed feature 
# are normally distributed.


#Overall, for alpha-diversity:

# ANOVA, t-test, or general linear models with the normal distribution 
# are used when the data is roughly normal. Transforming the data to 
# achieve a normal distribution could also be completed.
#
# Kruskal-Wallis, Wilcoxon rank sum test, or general linear models 
# with another distribution are used when the data is not normal or if 
# the n is low, like less than 30.

# Our main variables of interest are

# body site: colon, cecum
# Diet: AC , NAC
# Phenotype: CKD , NL
#group:

## Categorical variables
# Now that we know which tests can be used, let's run them. 

## Normally distributed metrics

# Since it's the closest to normalcy, we will use **Evenness** as an 
#example. First, we will test body site, which is a categorical variable 
# with more than 2 levels. Thus, we run ANOVA. If age were only two 
# levels, we could run a t-test

# Does body site impact the Evenness of the microbiota?

#Run the ANOVA and save it as an object

#Shannon
aov.Shannon.Diet = aov(shannon_entropy ~ Diet, data=meta)
aov.Shannon.Phenotype = aov(shannon_entropy ~ Phenotype, data=meta)
aov.Shannon.Group = aov(shannon_entropy ~ Group, data=meta)
#Call for the summary of that ANOVA, which will include P-values
summary(aov.Shannon.Diet)
summary(aov.Shannon.Phenotype)
summary(aov.Shannon.Group)
TukeyHSD(aov.Shannon.Group)
#Shannon for both diet and phenotype are not significant.

#observed_features
aov.observed_features.Diet = aov(observed_features ~ Diet, data=meta)
aov.observed_features.Phenotype = aov(observed_features ~ Phenotype, data=meta)
aov.observed_features.Group= aov(observed_features ~ Group, data=meta)

summary(aov.observed_features.Diet)
summary(aov.observed_features.Phenotype)
summary(aov.observed_features.Group)

#observed_features for both diet and phenotype are not significant.
#To do all the pairwise comparisons between groups and correct for multiple comparisons, we run Tukey's honest significance test of our ANOVA.

TukeyHSD(aov.observed_features.Group)


#Plot

shannon_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(Diet) %>%   # the grouping variable
  summarise(mean_shannon = mean(shannon_entropy),  # calculates the mean of each group
            sd_shannon = sd(shannon_entropy), # calculates the standard deviation of each group
            n_shannon = n(),  # calculates the sample size per group
            se_shannon = sd(shannon_entropy)/sqrt(n())) # calculates the standard error of each group

# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

shannon_se <- ggplot(shannon_summary, aes(Diet, mean_shannon, fill = Diet)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="shannon  ± s.e.", x = "") 

ggsave("output/Dietshannon_se.png", shannon_se, height = 2.5, width = 3)






observed_features_summary <- meta %>% # the names of the new data frame and the data frame to be summarised
  group_by(Diet) %>%   # the grouping variable
  summarise(mean_observed_features = mean(observed_features),  # calculates the mean of each group
            sd_observed_features = sd(observed_features), # calculates the standard deviation of each group
            n_observed_features = n(),  # calculates the sample size per group
            se_observed_features = sd(observed_features)/sqrt(n())) # calculates the standard error of each group

# We can now make a bar plot of means vs body site, with standard 
# deviations or standard errors as the error bar. The following code 
# uses the standard deviations.

observed_features_se <- ggplot(observed_features_summary, aes(Diet, mean_observed_features, fill = Diet)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_observed_features - se_observed_features, ymax = mean_observed_features + se_observed_features), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="observed_features  ± s.e.", x = "") 

ggsave("output/Dietobserved_features_se.png", observed_features_se, height = 2.5, width = 3)
# We clearly see that the evenness between hands and gut are different. 
# When we plot the data, we see that evenness decreases in the gut 
# compared to palms.


#NONparametric (non-normal) results
kruskal.test(faith_pd ~ Diet, data=meta)
kruskal.test(faith_pd ~ Phenotype, data=meta)
kruskal.test(faith_pd ~ Group, data=meta)

kruskal.test(pielou_evenness ~ Diet, data=meta)
kruskal.test(pielou_evenness ~ Phenotype, data=meta)
kruskal.test(pielou_evenness ~ Group, data=meta)


We can test pairwise within the age groups with Wilcoxon Rank Sum 
# Tests. This test has a slightly different syntax than our other tests

pairwise.wilcox.test(meta$faith_pd, meta$Group, p.adjust.method="BH")



#Plot
boxplot(faith_pd ~ Diet, data=meta, ylab="Faith phylogenetic diversity")

# or with ggplot2

faith_pd_boxplot <- ggplot(meta, aes(Diet, faith_pd)) + 
  geom_boxplot(aes(color = Diet)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith Phylogenetic Diversity", x = "") 
ggsave("output/Dietpd.png", faith_pd_boxplot, height = 3, width = 3)




boxplot(pielou_evenness ~ Diet, data=meta, ylab="pielou_evenness_diversity")

# or with ggplot2

pielou_evenness_boxplot <- ggplot(meta, aes(Diet, pielou_evenness)) + 
  geom_boxplot(aes(color = Diet)) + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="pielou_evenness", x = "") 
ggsave("output/Dietpielou.png", pielou_evenness_boxplot, height = 3, width = 3)
