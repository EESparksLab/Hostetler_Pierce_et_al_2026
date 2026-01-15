#Script associated with Hostetler and Pierce et al., 2026####
#Last updated on 1/6/2026; 12/16/2025

#As of 1/9/26:
## work on updating Lodging and RSS figure (density plot) + add lines to other plot separating low, med, high RSS
##work on emergence 2022 data and stats 
#As of 1/14/26:
##just finished figure 3, need to update text in manuscript, update root angle data in supplemental table, 
library(ggplot2)
citation("ggplot2")
packageVersion("ggplot2")
library(dplyr)
citation("dplyr")
packageVersion("dplyr")
library(agricolae) #HSD Test
citation("agricolae")
packageVersion("agricolae")
library(rcompanion) #transformTukey Test
citation("rcompanion")
packageVersion("rcompanion")
library(ggtext) #coloring x axis
citation("ggtext")
packageVersion("ggtext")
library(lme4) #H2
citation("lme4")
packageVersion("lme4")
library(PerformanceAnalytics) #Correlations Matrix
citation("PerformanceAnalytics")
packageVersion("PerformanceAnalytics")
library(corrplot) #Correlations
citation("corrplot")
packageVersion("corrplot")
library(factoextra) #get_pca_var
citation("factoextra")
packageVersion("factoextra")
library(ggfortify) #autplot
citation("ggfortify")
packageVersion("ggfortify")
library(cowplot) #autplot
citation("cowplot")
packageVersion("cowplot")
library(tidyr) #pivotwider
citation("tidyr")
packageVersion("tidyr")
library(stringr) 
citation("stringr")
packageVersion("stringr")
library("janitor")
citation("janitor")
packageVersion("janitor")
library(tidymodels)
citation("tidymodels")
packageVersion("tidymodels")
library(ranger)
citation("ranger")
packageVersion("ranger")
library(kknn)
citation("kknn")
packageVersion("kknn")
library(randomForest)
citation("randomForest")
packageVersion("randomForest")

cat("\014")
rm(list=ls()) 
ls() 
setwd(dir = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Data/")
getwd()

#The root system stiffness varies among maize inbred genotypes####
data = read.csv("SMURFDatabase.csv", header = TRUE, na.strings = "NA")
df = subset(data, Year == "2021")
head(df)
##SMURF Classifications
X = unique(df$Genotype)
df1 = df
colnames(df1)
df1 = df1[,c(1,15)]
head(df1)
df1 = df1 %>%
  group_by(Genotype) %>%
  summarise(across(1, list(mean = ~mean(. , na.rm = TRUE), 
                           sd = ~sd(. , na.rm = TRUE))) 
  )
df1=df1[,c(1,2)]
df1 = as.data.frame(df1)
head(df1)
data12 = df1
colnames(data12)[2] = "SMURF_Slope_Mean"
colnames(data12)
data12$SMURF_Slope_Mean = scale(data12$SMURF_Slope_Mean, center=TRUE, scale=TRUE)
data12$SMURF_Slope_Mean = as.numeric(data12$SMURF_Slope_Mean)
mean = mean(data12$SMURF_Slope_Mean)
sd = sd(data12$SMURF_Slope_Mean)
high = mean+sd
low = mean-sd
data12$SMURF_cat = data12$SMURF_Slope_Mean
data12$SMURF_cat = as.numeric(data12$SMURF_cat)
for (i in 1:length(data12$SMURF_cat)){
  if (data12$SMURF_Slope_Mean[i] <= low) {
    data12$SMURF_cat[i] = "low"
  } else if (data12$SMURF_Slope_Mean[i] < high) {
    data12$SMURF_cat[i] = "average"
  } else if (data12$SMURF_Slope_Mean[i] >= high) {
    data12$SMURF_cat[i] = "high"
  }
}
head(data12)
data12 = data12[,c(1,3)]
#write.csv(data12, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Data/SMURFClassifications.csv", quote = F, row.names = F)
df = merge(df, data12, by = "Genotype")
data = merge(data, data12, by = "Genotype")
head(df)

ordered_genotypes = with(df, reorder(Genotype, line_raw_slope_N.m, function(x) mean(x, na.rm = TRUE)))
genotype_order = levels(ordered_genotypes)
genotype_order = as.list(genotype_order)
#saveRDS(genotype_order, "genotype_order.rds")
#ordered_genotypes = readRDS("genotype_order.rds")
df$Genotype_ordered = ordered_genotypes
head(df)
df$SMURF_cat = factor(df$SMURF_cat, levels = c("low","average", "high"))
Figure1=ggplot(df, aes(x = Genotype_ordered, y = line_raw_slope_N.m, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("Root System Stiffness (N/m)")+
  xlab("Maize Inbred Line")+
  scale_y_continuous(limits=c(0,8000), breaks=seq(0,8000,1000))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/Figure1.pdf", width=7, height=5)
Figure1
#dev.off()
head(df)

#Experimental Design Summary
summary_df = df %>%
  group_by(Genotype, Plot) %>%
  summarise(
    Replicates_per_Plot = n_distinct(Replicate),
    .groups = "drop"
  ) %>%
  group_by(Genotype) %>%
  summarise(
    Number_of_Plots = n_distinct(Plot),
    Avg_Replicates = mean(Replicates_per_Plot),
    Min_Replicates = min(Replicates_per_Plot),
    Max_Replicates = max(Replicates_per_Plot),
    .groups = "drop"
  )
summary_df = as.data.frame(summary_df)
head(summary_df)

#Stats
attach(df)
head(df)
lm_x = lm(line_raw_slope_N.m ~ Genotype)
anova(lm_x)
lm_x_aov=aov(lm_x) 
X = HSD.test(lm_x_aov, trt = c("Genotype"), console = TRUE)
Means = X$means
#write.csv(Means, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS3a.csv", row.names = TRUE)
par(mfrow=c(2,2))
plot(lm_x)
par(mfrow=c(2,1))
plot(df$line_raw_slope_N.m)
hist(df$line_raw_slope_N.m)
resid = residuals(object = lm_x_aov)
shapiro.test(x=resid)
#transformation data to meet assumptions
par(mfrow=c(2,2))
df$tukey = transformTukey(
  df$line_raw_slope_N.m,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 1,
  returnLambda = FALSE
)
hist(df$tukey)
lm_x2 = lm(df$tukey ~ Genotype)
Z = anova(lm_x2)
Z = as.data.frame(Z)
Z
#write.csv(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS2.csv", row.names = TRUE)
lm_x2_aov2=aov(lm_x2) 
Y = HSD.test(lm_x2_aov2, trt = c("Genotype"), console = TRUE)
Groups = Y$groups
#write.csv(Groups, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS3b.csv", row.names = TRUE)
df$tukey = NULL
detach(df)

rm(list = setdiff(ls(), c("data", "df", "ordered_genotypes")))
data1 = data %>%
  group_by(Genotype) %>%
  filter(n_distinct(Year) == 2) %>%
  ungroup()
data1 = as.data.frame(data1)
data1$Year = as.character(data1$Year)
genotype_levels = levels(ordered_genotypes)
data1$Genotype_ordered = factor(data1$Genotype, levels = genotype_levels)
data1$SMURF_cat = factor(data1$SMURF_cat, levels = c("low","average", "high"))
FigureS1=ggplot(data1, aes(x = Genotype_ordered, y = line_raw_slope_N.m, fill = Year)) +
  geom_boxplot()+
  scale_fill_manual(values = c("seagreen","khaki"))+
  labs(fill = "Year")+
  ylab("Root System Stiffness (N/m)")+
  xlab("Maize Inbred Line")+
  scale_y_continuous(limits=c(0,8000), breaks=seq(0,8000,1000))+
  theme_bw() +
  theme(
    axis.text.x = element_markdown(angle = 60, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS1.pdf", width=11, height=8.5)
FigureS1
#dev.off()
#Experimental Design Summary
head(data1)
df1 = subset(data1, Year == "2020")
summary_df = df1 %>%
  group_by(Genotype, Plot) %>%
  summarise(
    Replicates_per_Plot = n_distinct(Replicate),
    .groups = "drop"
  ) %>%
  group_by(Genotype) %>%
  summarise(
    Number_of_Plots = n_distinct(Plot),
    Avg_Replicates = mean(Replicates_per_Plot),
    Min_Replicates = min(Replicates_per_Plot),
    Max_Replicates = max(Replicates_per_Plot),
    .groups = "drop"
  )
summary_df = as.data.frame(summary_df)
head(summary_df)

#Stats
attach(data1)
lm_x = lm(line_raw_slope_N.m ~ Genotype*Year)
anova(lm_x)
lm_x_aov=aov(lm_x) 
X = HSD.test(lm_x_aov, trt = c("Genotype"), console = TRUE)
Means = X$means
#write.csv(Means, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS5a.csv", row.names = TRUE)
X = HSD.test(lm_x_aov, trt = c("Year"), console = TRUE)
Means = X$means
#write.csv(Means, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS6a.csv", row.names = TRUE)
par(mfrow=c(2,2))
plot(lm_x)
par(mfrow=c(2,1))
plot(data1$line_raw_slope_N.m)
hist(data1$line_raw_slope_N.m)
resid = residuals(object = lm_x_aov)
shapiro.test(x=resid)
#transformation data1 to meet assumptions
par(mfrow=c(2,2))
data1$tukey = transformTukey(
  data1$line_raw_slope_N.m,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 1,
  returnLambda = FALSE
)
hist(data1$tukey)
lm_x2 = lm(data1$tukey ~ Genotype*Year)
Z = anova(lm_x2)
Z = as.data.frame(Z)
#write.csv(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS4.csv", row.names = TRUE)
lm_x2_aov2=aov(lm_x2) 
Y = HSD.test(lm_x2_aov2, trt = c("Genotype"), console = TRUE)
Groups = Y$groups
#write.csv(Groups, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS5b.csv", row.names = TRUE)
Y = HSD.test(lm_x2_aov2, trt = c("Year"), console = TRUE)
Groups = Y$groups
#write.csv(Groups, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS6b.csv", row.names = TRUE)
resid = residuals(object = lm_x2_aov2)
shapiro.test(x=resid)
data1$tukey = NULL
detach(data1)

head(data1)
geno_averages = data1 %>%
  group_by(Genotype, Year) %>%
  summarise(mean_slope = mean(line_raw_slope_N.m, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Year, values_from = mean_slope)

fit = lm(`2021` ~ `2020`, data = geno_averages)
r2_value = summary(fit)$r.squared
r_value = cor(geno_averages$`2020`, geno_averages$`2021`, use = "complete.obs")
FigureS2 = ggplot(geno_averages, aes(x = `2020`, y = `2021`)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black", fullrange = TRUE) +
  geom_text(aes(label = Genotype), color="gray50",
            vjust = -1, show.legend = FALSE )+
  labs(x = "Average Root System Stiffness (N/m) in 2020", y = "Average Root System Stiffness (N/m) in 2021") +
  scale_y_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 1000)) +
  scale_x_continuous(limits = c(0, 5000), breaks = seq(0, 5000, 1000)) +
  annotate("text",
    x = 100, y = 5000,
    label = paste0("r = ", round(r_value, 2), "\nR² = ", round(r2_value, 2)),
    size = 4, fontface = "bold"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 10)
  )

#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS2.pdf", width=7, height=5)
FigureS2
#dev.off()

head(geno_averages)
geno_ranked = geno_averages %>%
  mutate(
    rank_2020 = min_rank(`2020`),
    rank_2021 = min_rank(`2021`)
  )

geno_ranked
geno_ranked = geno_ranked[,c(1,4,5)]
geno_ranked = geno_ranked %>%
  pivot_longer(
    cols = c(rank_2020, rank_2021),
    names_to = c("Rank", "Year"),
    names_pattern = "(.*)_(\\d+)$"
  )
geno_ranked = geno_ranked[,c(1,3,4)]
colnames(geno_ranked)[3] = "Rank"
geno_ranked
geno_ranked$Year = as.character(geno_ranked$Year)
geno_ranked$Genotype_ordered = factor(geno_ranked$Genotype, levels = genotype_levels)
FigureS3 = ggplot(geno_ranked, aes(x = Genotype_ordered, y = Rank, color = Year)) +
  geom_point()+
  scale_color_manual(values = c("seagreen","khaki3"))+
  labs(x = "Genotype", y = "Rank") +
  scale_y_continuous(limits = c(1, 13), breaks = seq(1, 13, 1)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 10)
  )
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS3.pdf", width=7, height=5)
FigureS3
#dev.off()
rm(list = setdiff(ls(), c("ordered_genotypes")))

data = read.csv("SMURFDatabase.csv", header = TRUE, na.strings = "NA")
data = subset(data, Year == "2021")
data_pl = read.csv("LodgingData_55&Friends.csv", header = TRUE, na.strings = "NA")
data = merge(data, data_pl, by = "Genotype")
df2 = read.csv("SMURFClassifications.csv", header = TRUE, na.strings = "NA")
data = merge(data, df2, by = "Genotype")
head(data)
data$PL_cat = factor(data$PL_cat, levels = c("low","mid", "high"))
FigureS4B = ggplot(data, aes(x = PL_cat, y = line_raw_slope_N.m, fill = PL_cat)) +
  geom_boxplot()+
  scale_fill_manual(values = c("low" = "lightblue1", "mid" = "lightblue3", "high" = "lightblue4"), 
                    labels = c("low" = "Low", "mid" = "Moderate", "high" = "Severe")) +
  ylab("Root System Stiffness (N/m)")+
  xlab("Percent Lodging Category")+
  scale_x_discrete(labels = c("low" = "Low", "mid" = "Moderate", "high" = "Severe")) +
  labs(fill = "Percent Lodging Category") +
  scale_y_continuous(limits=c(0,8000), breaks=seq(0,8000,1000))+
  theme_bw()+
  theme(legend.position = "none",
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
FigureS4B
#Stats
attach(data)
head(data)
lm_x = lm(line_raw_slope_N.m ~ PL_cat)
anova(lm_x)
lm_x_aov=aov(lm_x) 
X = HSD.test(lm_x_aov, trt = c("PL_cat"), console = TRUE)
Means = X$means
Means
#write.csv(Means, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS8a.csv", row.names = TRUE)
par(mfrow=c(2,2))
plot(lm_x)
par(mfrow=c(2,1))
plot(data$line_raw_slope_N.m)
hist(data$line_raw_slope_N.m)
resid = residuals(object = lm_x_aov)
shapiro.test(x=resid)
#transformation data to meet assumptions
par(mfrow=c(2,2))
data$tukey = transformTukey(
  data$line_raw_slope_N.m,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 1,
  returnLambda = FALSE
)
hist(data$tukey)
par(mfrow=c(1,1))
lm_x2 = lm(data$tukey ~ PL_cat)
Z = anova(lm_x2)
Z = as.data.frame(Z)
#write.csv(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS7.csv", row.names = TRUE)
lm_x2_aov2=aov(lm_x2) 
resid = residuals(object = lm_x2_aov2)
shapiro.test(x=resid)
Y = HSD.test(lm_x2_aov2, trt = c("PL_cat"), console = TRUE)
Groups = Y$groups
#write.csv(Groups, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS8b.csv", row.names = TRUE)
data$tukey = NULL
detach(data)

genotype_levels = levels(ordered_genotypes)
data$Genotype_ordered = factor(data$Genotype, levels = genotype_levels)
head(data)
unique(data$Year)
FigureS4A = ggplot(data, aes(x = Genotype_ordered, y = line_raw_slope_N.m, fill = PL_cat)) +
  geom_boxplot()+
  scale_fill_manual(values = c("low" = "lightblue1", "mid" = "lightblue3", "high" = "lightblue4"), 
                    labels = c("low" = "Low", "mid" = "Moderate", "high" = "Severe")) +
  ylab("Root System Stiffness (N/m)")+
  xlab("Percent Lodging Category")+
  scale_x_discrete(labels = c("low" = "Low", "mid" = "Moderate", "high" = "Severe")) +
  labs(fill = "Percent Lodging Category") +
  scale_y_continuous(limits=c(0,8000), breaks=seq(0,8000,1000))+
  theme_bw()+
  theme(legend.position = "none",
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
FigureS4A
head(data)
data_pct = data %>%
  count(SMURF_cat, PL_cat) %>%
  group_by(SMURF_cat) %>%
  mutate(
    prop = n / sum(n),
    label = scales::percent(prop, accuracy = 1),
    label = ifelse(prop < 0.05, "", label))
data_pct$PL_cat = factor(data_pct$PL_cat, levels = c("high","mid", "low"))
FigureS4C = ggplot(data_pct, aes(x = SMURF_cat, y = n, fill = PL_cat)) +
  geom_col(position = "fill") +
  geom_text(aes(label = label),position = position_fill(vjust = 0.5),size = 3,color = "black") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("low" = "lightblue1","mid" = "lightblue3","high" = "lightblue4"),
    labels = c("low" = "Low","mid" = "Moderate","high" = "Severe")) +
  scale_x_discrete(labels = c("low" = "Low", "average" = "Moderate","high" = "High")) +
  labs(
    y = "",
    x = "Root System Stiffness Category",
    fill = "Root Lodging Category") +
  theme_bw()+
  theme(legend.position = "none",
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.text = element_text(size=10),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
FigureS4C
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS4.pdf", width=8, height=6)
ggdraw() +
  draw_plot(FigureS4A, x = 0, y = 0, width = 1, height = 0.50) +
  draw_plot(FigureS4B, x = 0, y = 0.50, width = 0.50, height = 0.50) +
  draw_plot(FigureS4C, x = 0.50, y = 0.50, width = 0.50, height = 0.50) +
  draw_plot_label(label = c("A", "B", "C"), 
                  size = 10,
                  x = c(0,0.5,0), 
                  y = c(1,1,0.55))

#dev.off()
rm(list=ls()) 

#The root system architecture varies among maize inbreds####
par(mfrow=c(1,1))
data = read.csv("RhizovisionData2017_55&Friends.csv", header = TRUE, na.strings = "NA")
head(data)
summary_table = data %>%
  group_by(Genotype) %>%
  summarise(
    n_Plots = n_distinct(Plot),
    Plants_per_Plot = paste(unique(
      sapply(split(Plant, Plot), function(x) length(unique(x)))
    ), collapse = ", ")
  )
colnames(data)
data = data[,c(1,6:26)]
df1 = data %>%
  group_by(Genotype) %>%
  summarise(across(1:21, list(mean = ~mean(. , na.rm = TRUE)))
  )
df1 = as.data.frame(df1)
head(df1)
df2 = read.csv("SMURFClassifications.csv", header = TRUE, na.strings = "NA")
data_pl = read.csv("LodgingData_55&Friends.csv", header = TRUE, na.strings = "NA")
df1 = merge(df1, data_pl, by = "Genotype")
df1 = merge(df1, df2, by = "Genotype")

pca = prcomp(df1[,2:22], scale. = TRUE)
pca2 = as.data.frame(pca$rotation)
pca3 = as.data.frame(pca$x)
colnames(df1)
data2 = cbind(df1, pca3)
head(data2)
#write.table(pca2, file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS11.csv", sep = "," , row.names = TRUE, col.names = TRUE, append = FALSE)
#write.table(pca3, file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS12.csv", sep = "," , row.names = TRUE, col.names = TRUE, append = FALSE)
var = get_pca_var(pca)
data2$SMURF_cat = as.character(data2$SMURF_cat)
data2$SMURF_cat = factor(data2$SMURF_cat, levels = c("low","average", "high"))
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS6.pdf", width=7, height=5)
autoplot(pca, data = data2, colour="SMURF_cat", 
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.colour = 'black', loadings.label.size = 3, loadings.label.vjust=-.1,  loadings.label.hjust=-.05,
         scale=TRUE, frame=TRUE)+
  labs(fill = "RSS Category", color = "RSS Category") +
  scale_fill_manual(values = c("low" = "pink4", "average" = "steelblue2", "high" = "burlywood"), 
                    labels = c("low" = "Low", "average" = "Moderate", "high" = "High")) +
  scale_color_manual(values = c("low" = "pink4", "average" = "steelblue2", "high" = "burlywood"), 
                     labels = c("low" = "Low", "average" = "Moderate", "high" = "High")) +
  theme_bw()
#dev.off()

#Bar Chart Genotype Distributions
data4 = data
head(data4)
df2 = read.csv("SMURFClassifications.csv", header = TRUE, na.strings = "NA")
data_pl = read.csv("LodgingData_55&Friends.csv", header = TRUE, na.strings = "NA")
data4 = merge(data4, data_pl, by = "Genotype")
data4 = merge(data4, df2, by = "Genotype")
genotype_order = readRDS("genotype_order.rds")
data4$Genotype = factor(data4$Genotype, levels = genotype_order)
head(data4)

#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS5.pdf", width=11, height=8.5)
i = 0
for (trait in 2:22){
  i = i+1
  alp = LETTERS[(i)]
  data.na = data4[,c(1,trait,25)]
  col = colnames(data4[trait])
  print(ggplot(data = data.na, aes(x=Genotype, y=data.na[,2], fill = SMURF_cat))+ 
          geom_boxplot() +
          ylab(col)+
          xlab("Maize Inbred Line")+
          ggtitle(alp)+
          scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), limits = c("low", "average", "high"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
          labs(fill = "RSS Category")+
          theme_bw()+
          expand_limits(y = 0)+
          theme(
            axis.text.x = element_text(size=8, angle = 60, hjust=1),
            axis.text.y = element_text(size=8),
            axis.text = element_text(size=8),
            axis.title = element_text(face="bold", size=10),
            axis.title.x = element_text(face="bold", size=10),
            axis.title.y = element_text(face="bold", size=10))
  )
}
#dev.off()
head(data4)
#Set up data frame to fill in
SW = matrix(NA,nrow=2,ncol=2)
rownames(SW) = c("W","pvalue")
SW[1,1] = "W"
SW[2,1] = "pvalue"
Z = "" #placeholder
colnames(data4)
df2 = data4[,c(1:22)]
head(df2)
colnames(df2)
for (i in 2:22){
  df3 = df2[,c(1,i)]
  a = colnames(df3)[2]
  colnames(df3)[2] = "trait"
  aov = lm(trait ~ Genotype, data = df3)
  hsd = aov(aov)              
  resid = residuals(object = hsd)
  shap = shapiro.test(x=resid)
  SW[1,2] = shap$statistic 
  SW[2,2] = shap$p.value 
  colnames(SW)[2] = a 
  SW = as.data.frame(SW)
  if (shap$p.value > 0.05){  
    print(shapiro.test(x=resid)$p.value)
    aov = lm(trait ~ Genotype, data = df3)
    write.table(a, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS9.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    Z = ""
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS9.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS9.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)  
    hsd_test = HSD.test(hsd, trt = c("Genotype"), console = TRUE)  
    means1 = as.data.frame(hsd_test$means)
    colnames(means1)[1] = a 
    write.table(means1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS10a.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS10a.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = a 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS10b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS10b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
  #for data that needs normalization
  if (shap$p.value < 0.05){
    print(shapiro.test(x=resid)$p.value)
    norm_col_name = paste0("norm_", a)
    par(mfrow=c(3,1))
    df3$norm_col_name = transformTukey(df3$trait,
                                       start = -10,
                                       end = 10,
                                       int = 0.025,
                                       plotit = TRUE, #can be false 
                                       verbose = FALSE,
                                       quiet = FALSE,
                                       statistic = 1,
                                       returnLambda = FALSE  )
    aov = lm(norm_col_name ~ Genotype, data = df3)
    col = norm_col_name
    write.table(col, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS9.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS9.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS9.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)              
    hsd_test = HSD.test(hsd, trt = c("Genotype"), console = TRUE)  
    means1 = as.data.frame(hsd_test$means)
    colnames(means1)[1] = col 
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = col 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS10b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS10b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    aov = lm(trait ~ Genotype, data = df3)
    hsd = aov(aov)  
    hsd_test = HSD.test(hsd, trt = c("Genotype"), console = TRUE)  
    means1 = as.data.frame(hsd_test$means)
    colnames(means1)[1] = a 
    write.table(means1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS10a.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS10a.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
}


head(data4)
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS7.pdf", width=8, height=6)
i = 0
for (trait in 2:22){
  i = i+1
  alp = LETTERS[(i)]
  data.na = data4[,c(1,trait,25)]
  col = colnames(data4[trait])
  data.na$SMURF_cat = factor(data.na$SMURF_cat, levels = c("low","average", "high"))
  print(ggplot(data = data.na, aes(x=SMURF_cat, y=data.na[,2], fill = SMURF_cat))+ 
          geom_boxplot() +
          ylab(col)+
          xlab(NULL)+
          scale_x_discrete(labels = c("low" = "Low", 
                                      "average" = "Moderate", 
                                      "high" = "High"))+
          labs(fill = "RSS Category")+
          ggtitle(alp)+
          scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
          theme_bw()+
          expand_limits(y = 0)+
          theme(#legend.position = "none",
            axis.text.x = element_text(size=10),
            axis.text.y = element_text(size=10),
            axis.text = element_text(size=10),
            axis.title = element_text(face="bold", size=12),
            axis.title.x = element_text(face="bold", size=12),
            axis.title.y = element_text(face="bold", size=12))
  )
}
#dev.off()

#Set up data frame to fill in 2
SW = matrix(NA,nrow=2,ncol=2)
rownames(SW) = c("W","pvalue")
SW[1,1] = "W"
SW[2,1] = "pvalue"
Z = "" #placeholder
colnames(data4)
df2 = data4[,c(2:22,25)]
head(df2)
colnames(df2)
for (i in 1:21){
  df3 = df2[,c(22,i)]
  a = colnames(df3)[2]
  colnames(df3)[2] = "trait"
  aov = lm(trait ~ SMURF_cat, data = df3)
  hsd = aov(aov)              
  resid = residuals(object = hsd)
  shap = shapiro.test(x=resid)
  SW[1,2] = shap$statistic 
  SW[2,2] = shap$p.value 
  colnames(SW)[2] = a 
  SW = as.data.frame(SW)
  if (shap$p.value > 0.05){  
    print(shapiro.test(x=resid)$p.value)
    aov = lm(trait ~ SMURF_cat, data = df3)
    write.table(a, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS13.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    Z = ""
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS13.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS13.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)  
    hsd_test = HSD.test(hsd, trt = c("SMURF_cat"), console = TRUE)  
    means1 = as.data.frame(hsd_test$means)
    colnames(means1)[1] = a 
    write.table(means1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS14a.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS14a.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = a 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS14b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS14b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
  #for data that needs normalization
  if (shap$p.value < 0.05){
    print(shapiro.test(x=resid)$p.value)
    norm_col_name <- paste0("norm_", a)
    par(mfrow=c(3,1))
    df3$norm_col_name = transformTukey(df3$trait,
                                       start = -10,
                                       end = 10,
                                       int = 0.025,
                                       plotit = TRUE, #can be false 
                                       verbose = FALSE,
                                       quiet = FALSE,
                                       statistic = 1,
                                       returnLambda = FALSE  )
    aov = lm(norm_col_name ~ SMURF_cat, data = df3)
    col = norm_col_name
    write.table(col, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS13.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS13.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS13.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)              
    hsd_test = HSD.test(hsd, trt = c("SMURF_cat"), console = TRUE)  
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = col 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS14b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS14b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    aov = lm(trait ~ SMURF_cat, data = df3)
    hsd = aov(aov)  
    hsd_test = HSD.test(hsd, trt = c("SMURF_cat"), console = TRUE)  
    means1 = as.data.frame(hsd_test$means)
    colnames(means1)[1] = a 
    write.table(means1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS14a.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS14a.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
}

colnames(data4)
data4$SMURF_cat = factor(data4$SMURF_cat, levels = c("low","average", "high"))
Figure2A = ggplot(data = data4, aes(x=SMURF_cat, y=Solidity, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Solidity")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2B = ggplot(data = data4, aes(x=SMURF_cat, y=TotalRL_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Total Root Length (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2C = ggplot(data = data4, aes(x=SMURF_cat, y=LowerRootArea_cm2, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Lower Root Area (cm2)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2D = ggplot(data = data4, aes(x=SMURF_cat, y=Depth_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Depth (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2E = ggplot(data = data4, aes(x=SMURF_cat, y=MaxWidth_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Max Width (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2F = ggplot(data = data4, aes(x=SMURF_cat, y=NetworkArea_cm2, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Network Area (cm2)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2G = ggplot(data = data4, aes(x=SMURF_cat, y=ConvexArea_cm2, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Convex Area (cm2)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2H = ggplot(data = data4, aes(x=SMURF_cat, y=MediumAngleFrequency, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Med. Angle Freq.")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2I = ggplot(data = data4, aes(x=SMURF_cat, y=SteepAngleFrequency, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Steep Angle Freq.")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2J = ggplot(data = data4, aes(x=SMURF_cat, y=AverageDiameter_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Average Diameter (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2K = ggplot(data = data4, aes(x=SMURF_cat, y=MedianDiameter_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Median Diameter (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2L = ggplot(data = data4, aes(x=SMURF_cat, y=MaxDiameter_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Max Diameter (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2M = ggplot(data = data4, aes(x=SMURF_cat, y=Perimeter_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Perimeter (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2N = ggplot(data = data4, aes(x=SMURF_cat, y=SurfaceArea_cm2, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Surface Area (cm2)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2O = ggplot(data = data4, aes(x=SMURF_cat, y=Volume_cm3, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Volume (cm3)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Mod.", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Mod.", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure2P = autoplot(pca, data = data2, colour="SMURF_cat", 
         loadings = FALSE, loadings.colour = 'black',
         loadings.label = FALSE, loadings.label.colour = 'black', loadings.label.size = 3, loadings.label.vjust=-.1,  loadings.label.hjust=-.05,
         scale=TRUE, frame=TRUE)+
  labs(fill = "RSS Category", color = "RSS Category") +
  scale_fill_manual(values = c("low" = "pink4", "average" = "steelblue2", "high" = "burlywood"), 
                    labels = c("low" = "Low", "average" = "Mod.", "high" = "High")) +
  scale_color_manual(values = c("low" = "pink4", "average" = "steelblue2", "high" = "burlywood"), 
                     labels = c("low" = "Low", "average" = "Mod.", "high" = "High")) +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/Figure2.pdf", width=8, height=6)
ggdraw() +
  draw_plot(Figure2P, x = 0, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(Figure2A, x = 0.25, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(Figure2B, x = 0.50, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(Figure2C, x = 0.75, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(Figure2D, x = 0, y = 0.50, width = 0.25, height = 0.25) +
  draw_plot(Figure2E, x = 0.25, y = 0.50, width = 0.25, height = 0.25) +
  draw_plot(Figure2F, x = 0.50, y = 0.50, width = 0.25, height = 0.25) +
  draw_plot(Figure2G, x = 0.75, y = 0.50, width = 0.25, height = 0.25) +
  draw_plot(Figure2H, x = 0, y = 0.25, width = 0.25, height = 0.25) +
  draw_plot(Figure2I, x = 0.25, y = 0.25, width = 0.25, height = 0.25) +
  draw_plot(Figure2J, x = 0.50, y = 0.25, width = 0.25, height = 0.25) +
  draw_plot(Figure2K, x = 0.75, y = 0.25, width = 0.25, height = 0.25) +
  draw_plot(Figure2L, x = 0, y = 0, width = 0.25, height = 0.25) +
  draw_plot(Figure2M, x = 0.25, y = 0, width = 0.25, height = 0.25) +
  draw_plot(Figure2N, x = 0.50, y = 0, width = 0.25, height = 0.25) +
  draw_plot(Figure2O, x = 0.75, y = 0, width = 0.25, height = 0.25) +
  draw_plot_label(label = c("A", "B", "C", "D"), 
                  size = 10,
                  x = c(0,0.25,0,0), 
                  y = c(1,1,0.75,0.50))
#dev.off()
rm(list=ls()) 


#Brace root architecture varies among maize inbreds
par(mfrow=c(1,1))
rm(list=ls()) 
data = read.csv("RootPhotoData_2019.csv", header = TRUE, na.strings = "NA")
head(data)
colnames(data)
data$root_angle_deg2 = 90 - data$root_angle_deg
colnames(data)[15] = "root_angle_deg"
head(data)
data = data[,c(1:9,15,11:14)]
data = data[,c(1,5:14)]
colnames(data)
df1 = data %>%
  group_by(Genotype) %>%
  summarise(across(1:10, list(mean = ~mean(. , na.rm = TRUE)))
  )
df1 = as.data.frame(df1)
head(df1)
df2 = read.csv("SMURFClassifications.csv", header = TRUE, na.strings = "NA")
data_pl = read.csv("LodgingData_55&Friends.csv", header = TRUE, na.strings = "NA")
df1 = merge(df1, data_pl, by = "Genotype")
df1 = merge(df1, df2, by = "Genotype")
pca = prcomp(df1[,c(2,7:11)], scale. = TRUE)
pca2 = as.data.frame(pca$rotation)
pca3 = as.data.frame(pca$x)
colnames(df1)
data2 = cbind(df1, pca3)
head(data2)
#write.table(pca2, file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS17.csv", sep = "," , row.names = TRUE, col.names = TRUE, append = FALSE)
#write.table(pca3, file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS18.csv", sep = "," , row.names = TRUE, col.names = TRUE, append = FALSE)
var = get_pca_var(pca)
data2$SMURF_cat = as.character(data2$SMURF_cat)
data2$SMURF_cat = factor(data2$SMURF_cat, levels = c("low","average", "high"))
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS9.pdf", width=7, height=5)
autoplot(pca, data = data2, colour="SMURF_cat", 
         loadings = TRUE, loadings.colour = 'black',
         loadings.label = TRUE, loadings.label.colour = 'black', loadings.label.size = 3, loadings.label.vjust=-.1,  loadings.label.hjust=-.05,
         scale=TRUE, frame=TRUE)+
  labs(fill = "RSS Category", color = "RSS Category") +
  scale_fill_manual(values = c("low" = "pink4", "average" = "steelblue2", "high" = "burlywood"), 
                    labels = c("low" = "Low", "average" = "Moderate", "high" = "High")) +
  scale_color_manual(values = c("low" = "pink4", "average" = "steelblue2", "high" = "burlywood"), 
                     labels = c("low" = "Low", "average" = "Moderate", "high" = "High")) +
  theme_bw()
#dev.off()

head(data)
data4 = data[,c(1:2,7:11)]
head(data4)

df2 = read.csv("SMURFClassifications.csv", header = TRUE, na.strings = "NA")
data_pl = read.csv("LodgingData_55&Friends.csv", header = TRUE, na.strings = "NA")
data4 = merge(data4, data_pl, by = "Genotype")
data4 = merge(data4, df2, by = "Genotype")
genotype_order = readRDS("genotype_order.rds")
data4$Genotype = factor(data4$Genotype, levels = genotype_order)
head(data4)
head(data4)
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS8-2.pdf", width=6, height=4 )
i = 0
for (trait in 2:7){
  i = i+1
  alp = LETTERS[(i)]
  data.na = data4[,c(1,trait,10)]
  col = colnames(data4[trait])
  data.na$SMURF_cat = factor(data.na$SMURF_cat, levels = c("low","average", "high"))
  print(ggplot(data = data.na, aes(x=SMURF_cat, y=data.na[,2], fill = SMURF_cat))+ 
          geom_boxplot() +
          ylab(col)+
          xlab(NULL)+
          scale_x_discrete(labels = c("low" = "Low", 
                                      "average" = "Moderate", 
                                      "high" = "High"))+
          labs(fill = "RSS Category") +
          ggtitle(alp)+
          scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
          theme_bw()+
          expand_limits(y = 0)+
          theme(#legend.position = "none",
                axis.text.x = element_text(size=8),
                axis.text.y = element_text(size=8),
                axis.text = element_text(size=8),
                axis.title = element_text(face="bold", size=10),
                axis.title.x = element_text(face="bold", size=10),
                axis.title.y = element_text(face="bold", size=10))
  )
}
#dev.off()
head(data4)
#Set up data frame to fill in 2
SW = matrix(NA,nrow=2,ncol=2)
rownames(SW) = c("W","pvalue")
SW[1,1] = "W"
SW[2,1] = "pvalue"
Z = "" #placeholder
colnames(data4)
df2 = data4[,c(2:7,10)]
head(df2)
colnames(df2)
for (i in 1:6){
  df3 = df2[,c(7,i)]
  a = colnames(df3)[2]
  colnames(df3)[2] = "trait"
  aov = lm(trait ~ SMURF_cat, data = df3)
  hsd = aov(aov)              
  resid = residuals(object = hsd)
  shap = shapiro.test(x=resid)
  SW[1,2] = shap$statistic 
  SW[2,2] = shap$p.value 
  colnames(SW)[2] = a 
  SW = as.data.frame(SW)
  if (shap$p.value > 0.05){  
    print(shapiro.test(x=resid)$p.value)
    aov = lm(trait ~ SMURF_cat, data = df3)
    write.table(a, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS15.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    Z = ""
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS15.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS15.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)  
    hsd_test = HSD.test(hsd, trt = c("SMURF_cat"), console = TRUE)  
    means1 = as.data.frame(hsd_test$means)
    colnames(means1)[1] = a 
    write.table(means1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS16a.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS16a.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = a 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS16b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS16b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
  #for data that needs normalization
  if (shap$p.value < 0.05){
    print(shapiro.test(x=resid)$p.value)
    norm_col_name = paste0("norm_", a)
    par(mfrow=c(3,1))
    df3$norm_col_name = transformTukey(df3$trait,
                                       start = -10,
                                       end = 10,
                                       int = 0.025,
                                       plotit = TRUE, #can be false 
                                       verbose = FALSE,
                                       quiet = FALSE,
                                       statistic = 1,
                                       returnLambda = FALSE  )
    aov = lm(norm_col_name ~ SMURF_cat, data = df3)
    col = norm_col_name
    write.table(col, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS15.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS15.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS15.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)              
    hsd_test = HSD.test(hsd, trt = c("SMURF_cat"), console = TRUE)  
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = col 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS16b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS16b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    aov = lm(trait ~ SMURF_cat, data = df3)
    anova = anova(aov)
    anova = as.data.frame(anova)
    hsd = aov(aov)  
    hsd_test = HSD.test(hsd, trt = c("SMURF_cat"), console = TRUE)  
    means1 = as.data.frame(hsd_test$means)
    colnames(means1)[1] = a 
    write.table(means1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS16a.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS16a.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
}
colnames(data4)
data4$SMURF_cat = factor(data4$SMURF_cat, levels = c("low","average", "high"))
Figure3B = ggplot(data = data4, aes(x=SMURF_cat, y=root_angle_deg, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Root Angle (deg)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure3C = ggplot(data = data4, aes(x=SMURF_cat, y=stalk_width_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Stalk Width (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure3D = ggplot(data = data4, aes(x=SMURF_cat, y=single_root_width_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Single Root Width (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure3E = ggplot(data = data4, aes(x=SMURF_cat, y=spread_width_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Spread Width (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure3F = ggplot(data = data4, aes(x=SMURF_cat, y=root_heightonstalk_cm, fill = SMURF_cat))+ 
  geom_boxplot() +
  xlab(NULL)+
  ylab("Root Height on Stalk (cm)")+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  theme_bw()+
  expand_limits(y = 0)+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure3A = autoplot(pca, data = data2, colour="SMURF_cat", 
                   loadings = FALSE, loadings.colour = 'black',
                   loadings.label = FALSE, loadings.label.colour = 'black', loadings.label.size = 3, loadings.label.vjust=-.1,  loadings.label.hjust=-.05,
                   scale=TRUE, frame=TRUE)+
  labs(fill = NULL, color = NULL) +
  scale_fill_manual(values = c("low" = "pink4", "average" = "steelblue2", "high" = "burlywood"), 
                    labels = c("low" = "Low", "average" = "Moderate", "high" = "High")) +
  scale_color_manual(values = c("low" = "pink4", "average" = "steelblue2", "high" = "burlywood"), 
                     labels = c("low" = "Low", "average" = "Moderate", "high" = "High")) +
  theme_bw()+
  theme(legend.position = "none",
        legend.text = element_text(size = 8))
Figure3A
#Emergence Data#
cat("\014")
data = read.csv("Emergence_Data_55Friends_2021.csv", header = TRUE, na.strings = "NA")
head(data)
unique(data$Genotype)
PlotInfo = data %>%
  group_by(Genotype) %>%
  summarise(Num_Plots = n_distinct(Plot)) %>%
  arrange(desc(Num_Plots))
PlotInfo = as.data.frame(PlotInfo)
PlotInfo
dap_variability = data %>%
  group_by(Genotype) %>%
  summarise(
    Mean_DAP = mean(DAP, na.rm = TRUE),
    SD_DAP = sd(DAP, na.rm = TRUE),
    Min_DAP = min(DAP, na.rm = TRUE),
    Max_DAP = max(DAP, na.rm = TRUE),
    Range_DAP = Max_DAP - Min_DAP
  )
dap_variability = as.data.frame(dap_variability)
dap_variability

data1 = read.csv("Emergence_Data_55Friends_2022.csv", header = TRUE, na.strings = "NA")
head(data1)
PlotInfo = data1 %>%
  group_by(Genotype) %>%
  summarise(Num_Plots = n_distinct(X2022.Field.Location)) %>%
  arrange(desc(Num_Plots))
PlotInfo = as.data.frame(PlotInfo)
PlotInfo
head(data1)
data1$PlantingDate = "5/23/22"
data1 = data1 %>%
  mutate(
    PlantingDate = as.Date(PlantingDate, format = "%m/%d/%y"),
    Date.of.emergence = as.Date(na_if(Date.of.emergence, ""), format = "%m/%d/%y"),
    DAP = as.numeric(Date.of.emergence - PlantingDate)
  )
str(data)
dap_variability = data %>%
  group_by(Genotype) %>%
  summarise(
    Mean_DAP = mean(DAP, na.rm = TRUE),
    SD_DAP = sd(DAP, na.rm = TRUE),
    Min_DAP = min(DAP, na.rm = TRUE),
    Max_DAP = max(DAP, na.rm = TRUE),
    Range_DAP = Max_DAP - Min_DAP
  )
dap_variability = as.data.frame(dap_variability)
dap_variability

colnames(data1)[3] = "Plot"
data1$Year = "2022"
data$Year = "2021"
data1 = data1[,c(1,3:9)]
data = data[,c(1:7,9)]
head(data1)
head(data)
data = rbind(data, data1)

df2 = read.csv("SMURFClassifications.csv", header = TRUE, na.strings = "NA")
data = merge(data, df2, by = "Genotype")
genotype_order = readRDS("genotype_order.rds")
data$Genotype = factor(data$Genotype, levels = genotype_order)
head(data)

#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS10.pdf", width=11, height=8.5)
ggplot(data, aes(x = Genotype, y = DAP, color = SMURF_cat, fill = Year)) +
  geom_boxplot()+
  ylab("Brace Root Emergence (dap)")+
  xlab("Maize Inbred")+
  scale_y_continuous(limits=c(0,65), breaks=seq(0,65,7))+
  scale_color_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  scale_fill_manual(values = c("gray20","gray80"))+
  labs(fill = "Year", color = "RSS Category")+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))
#dev.off()
attach(data)
head(data)
lm_x = lm(DAP ~ Genotype*Year)
anova(lm_x)
lm_x_aov=aov(lm_x) 
X = HSD.test(lm_x_aov, trt = c("Genotype","Year"), console = TRUE)
Means = X$means
#write.csv(Means, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS20b.csv", row.names = TRUE)
par(mfrow=c(2,2))
plot(lm_x)
par(mfrow=c(2,1))
plot(data$DAP)
hist(data$DAP)
resid = residuals(object = lm_x_aov)
shapiro.test(x=resid)
#transformation data to meet assumptions
par(mfrow=c(2,2))
data$tukey = transformTukey(
  data$DAP,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 1,
  returnLambda = FALSE
)
hist(data$tukey)
lm_x2 = lm(data$tukey ~ Genotype*Year)
Z = anova(lm_x2)
Z = as.data.frame(Z)
Z
#write.csv(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS19.csv", row.names = TRUE)
lm_x2_aov2=aov(lm_x2) 
Y = HSD.test(lm_x2_aov2, trt = c("Genotype","Year"), console = TRUE)
Groups = Y$groups
#write.csv(Groups, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS20a.csv", row.names = TRUE)
data$tukey = NULL
detach(data)

head(data)
data$SMURF_cat = factor(data$SMURF_cat, levels = c("low","average", "high"))
Figure3G = ggplot(data, aes(x = SMURF_cat, y = DAP, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("Brace Root Emergence (dap)")+
  xlab(NULL)+
  expand_limits(y = 0)+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
Figure3G
attach(data)
head(data)
lm_x = lm(DAP ~ SMURF_cat*Year)
anova(lm_x)
lm_x_aov=aov(lm_x) 
X = HSD.test(lm_x_aov, trt = c("SMURF_cat"), console = TRUE)
Means = X$means
#write.csv(Means, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS22a.csv", row.names = TRUE)
X = HSD.test(lm_x_aov, trt = c("Year"), console = TRUE)
Means = X$means
#write.csv(Means, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS22a-2.csv", row.names = TRUE)
par(mfrow=c(2,2))
plot(lm_x)
par(mfrow=c(2,1))
plot(data$DAP)
hist(data$DAP)
resid = residuals(object = lm_x_aov)
shapiro.test(x=resid)
#transformation data to meet assumptions
par(mfrow=c(2,2))
data$tukey = transformTukey(
  data$DAP,
  start = -10,
  end = 10,
  int = 0.025,
  plotit = TRUE,
  verbose = FALSE,
  quiet = FALSE,
  statistic = 1,
  returnLambda = FALSE
)
hist(data$tukey)
lm_x2 = lm(data$tukey ~ SMURF_cat*Year)
Z = anova(lm_x2)
Z = as.data.frame(Z)
Z
#write.csv(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS21.csv", row.names = TRUE)
lm_x2_aov2=aov(lm_x2) 
Y = HSD.test(lm_x2_aov2, trt = c("SMURF_cat"), console = TRUE)
Groups = Y$groups
#write.csv(Groups, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS22b.csv", row.names = TRUE)
data$tukey = NULL
detach(data)
head(data)
data$Leaf.count = as.numeric(data$Leaf.count)
FigureS11B = ggplot(data, aes(x = SMURF_cat, y = Leaf.count, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("Leaf Number (count)")+
  xlab(NULL)+
  expand_limits(y = 0)+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
FigureS11C = ggplot(data, aes(x = DAP, y = Leaf.count, color = SMURF_cat)) +
  geom_point(size=1)+
  ylab("Leaf Number (count)")+
  xlab("Brace Root Emergence (dap)")+
  expand_limits(y = 0, x=0)+
  scale_color_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(color = "RSS Category")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=8),
        axis.title.x = element_text(face="bold", size=8),
        axis.title.y = element_text(face="bold", size=8))
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS11.pdf", width=8, height=6)
ggdraw() +
  draw_plot(Figure3G, x = 0, y = 0.66, width = 0.5, height = 0.33) +
  draw_plot(FigureS11B, x = 0.50, y = 0.66, width = 0.5, height = 0.33) +
  draw_plot(FigureS11C, x = 0, y = 0.33, width = 0.66, height = 0.33) +
  draw_plot_label(label = c("A", "B", "C"), 
                  size = 10,
                  x = c(0,0.50,0), 
                  y = c(1,1,0.70))
#dev.off()

#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/Figure3-3.pdf", width=7, height=5)
ggdraw() +
  draw_plot(Figure3A, x = 0, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(Figure3B, x = 0.25, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(Figure3C, x = 0.50, y = 0.75, width = 0.25, height = 0.25) +
  draw_plot(Figure3D, x = 0, y = 0.50, width = 0.25, height = 0.25) +
  draw_plot(Figure3E, x = 0.25, y = 0.50, width = 0.25, height = 0.25) +
  draw_plot(Figure3F, x = 0.50, y = 0.50, width = 0.25, height = 0.25) +
  draw_plot(Figure3G, x = 0.75, y = 0.50, width = 0.25, height = 0.25) +
  draw_plot_label(label = c("A", "B", "C", "D","E","F","G"), 
                  size = 10,
                  x = c(0,0.25,0.50,0,0.25,0.50,0.75), 
                  y = c(1,1,1,0.75,0.75,0.75,0.75))
#dev.off()
head(data)
dap_variability = data %>%
  group_by(Genotype, Year) %>%
  summarise(
    Mean_DAP = mean(DAP, na.rm = TRUE),
    SD_DAP = sd(DAP, na.rm = TRUE),
    Min_DAP = min(DAP, na.rm = TRUE),
    Max_DAP = max(DAP, na.rm = TRUE),
    Range_DAP = Max_DAP - Min_DAP
  )
dap_variability = as.data.frame(dap_variability)
View(dap_variability)
data %>%
  summarise(
    mean_DAP = mean(DAP, na.rm = TRUE),
    sd_DAP   = sd(DAP, na.rm = TRUE)
  )

#Brace root mechanical variation is driven by geometry #####
rm(list=ls()) 
cat("\014")
data1 = read.csv("BRW_55Friends.csv", header = TRUE, na.strings = "NA")
head(data1)
data1 = data1 %>%
  group_by(Genotype, Year, Plot.ID, Plant.Number) %>%
  summarise(Brace.Root.Whorls.in.the.Soil = first(Brace.Root.Whorls.in.the.Soil), .groups = "drop")
data = read.csv("InstronData_55&Friends.csv", header = TRUE, na.strings = "NA")
head(data)
colnames(data1)[3] = "PlotNumber"
colnames(data1)[4] = "PlantNumber"
head(data1)
head(data)
data = merge(data, data1, by = c("Genotype","Year","PlotNumber","PlantNumber"))
head(data)
unique(data$Brace.Root.Whorls.in.the.Soil)
str(data)
data$Brace.Root.Whorls.in.the.Soil = as.numeric(data$Brace.Root.Whorls.in.the.Soil)
data$Whorl.or.Internode_Number = as.numeric(data$Whorl.or.Internode_Number)
data = data %>%
  mutate(
    Whorl.or.Internode_Number = as.numeric(as.character(Whorl.or.Internode_Number)),
    WR_ID = case_when(
      # 4 whorls
      Brace.Root.Whorls.in.the.Soil == 4 & Whorl.or.Internode_Number == 4 ~ "BRNode1",
      Brace.Root.Whorls.in.the.Soil == 4 & Whorl.or.Internode_Number == 3 ~ "BRNode2",
      Brace.Root.Whorls.in.the.Soil == 4 & Whorl.or.Internode_Number == 2 ~ "BRNode3",
      Brace.Root.Whorls.in.the.Soil == 4 & Whorl.or.Internode_Number == 1 ~ "BRNode4",
      is.na(Whorl.or.Internode_Number) & Brace.Root.Whorls.in.the.Soil == 4 ~ "BRNode4",
      
      # 3 whorls
      Brace.Root.Whorls.in.the.Soil == 3 & Whorl.or.Internode_Number == 3 ~ "BRNode1",
      Brace.Root.Whorls.in.the.Soil == 3 & Whorl.or.Internode_Number == 2 ~ "BRNode2",
      Brace.Root.Whorls.in.the.Soil == 3 & Whorl.or.Internode_Number == 1 ~ "BRNode3",
      is.na(Whorl.or.Internode_Number) & Brace.Root.Whorls.in.the.Soil == 3 ~ "BRNode3",
      
      # 2 whorls
      Brace.Root.Whorls.in.the.Soil == 2 & Whorl.or.Internode_Number == 2 ~ "BRNode1",
      Brace.Root.Whorls.in.the.Soil == 2 & Whorl.or.Internode_Number == 1 ~ "BRNode2",
      is.na(Whorl.or.Internode_Number) & Brace.Root.Whorls.in.the.Soil == 2 ~ "BRNode2",
      
      # 1 whorl
      Brace.Root.Whorls.in.the.Soil == 1 & Whorl.or.Internode_Number == 1 ~ "BRNode1",
      is.na(Whorl.or.Internode_Number) & Brace.Root.Whorls.in.the.Soil == 1 ~ "BRNode1",
      
      TRUE ~ NA_character_
    )
  )
str(data)
data = data[!is.na(data$WR_ID), ]
colnames(data)
colnames(data)[17] = "K_N.mm"
str(data)
data$a0_mm = (data$Vertical.Specimen.Diameter)/2
data$b0_mm = (data$Horizontal.Specimen.Diameter)/2
data$I_mm4 = (pi/4)*(((data$a0_mm)^3)*(data$b0_mm))
data$E_MPa = data$K_N.mm * ((17.5^3)/(48*(data$I_mm4)))
head(data)
colnames(data)
df = data[,c(1,2,21,17,22:25)]
colnames(df)
unique(df$Year)
genotype_order = readRDS("genotype_order.rds")
df$Genotype = factor(df$Genotype, levels = genotype_order)
df2 = read.csv("SMURFClassifications.csv", header = TRUE, na.strings = "NA")
data_pl = read.csv("LodgingData_55&Friends.csv", header = TRUE, na.strings = "NA")
df = merge(df, data_pl, by = "Genotype")
df = merge(df, df2, by = "Genotype")
head(df)
df$E_GPa = df$E_MPa/1000
colnames(df)
df = df[,c(1:7,12,9:11)]
nrow(df)
df = df %>%
  filter(if_all(K_N.mm:E_GPa, ~ . > 0))
df = subset(df, b0_mm < 10)
nrow(df)
df$SMURF_cat = factor(df$SMURF_cat, levels = c("low","average", "high"))
genotype_order = readRDS("genotype_order.rds")
df$Genotype = factor(df$Genotype, levels = genotype_order)
df$Year = as.character(df$Year)
head(df)
df = subset(df, WR_ID != "BRNode4")
unique(df$WR_ID)
df$WR_ID = factor(df$WR_ID, levels = rev(levels(factor(df$WR_ID))))
head(df)
id_labels = c(
  "BRNode3" = "Node 3",
  "BRNode2" = "Node 2",
  "BRNode1" = "Node 1 (Closest to Soil)"
)
head(df)
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS12.pdf", width=11, height=8.5)
ggplot(df, aes(x = Genotype, y = K_N.mm, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("K (N/mm)")+
  xlab("Maize Inbred Line")+
  ggtitle("A")+
  expand_limits(y = 0)+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
ggplot(df, aes(x = Genotype, y = a0_mm, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("a0 (mm)")+
  xlab("Maize Inbred Line")+
  ggtitle("B")+
  expand_limits(y = 0)+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
ggplot(df, aes(x = Genotype, y = b0_mm, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("b0 (mm)")+
  xlab("Maize Inbred Line")+
  ggtitle("C")+
  expand_limits(y = 0)+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
ggplot(df, aes(x = Genotype, y = I_mm4, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("I (mm4)")+
  xlab("Maize Inbred Line")+
  ggtitle("D")+
  expand_limits(y = 0)+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
ggplot(df, aes(x = Genotype, y = E_GPa, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("E (GPa)")+
  xlab("Maize Inbred Line")+
  ggtitle("E")+
  expand_limits(y = 0)+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
#dev.off()
#Set up data frame to fill in 2
SW = matrix(NA,nrow=2,ncol=2)
rownames(SW) = c("W","pvalue")
SW[1,1] = "W"
SW[2,1] = "pvalue"
Z = "" #placeholder
colnames(df)
df2 = df[,c(4:8,1,3)]
head(df2)
colnames(df2)
for (i in 1:5){
  df3 = df2[,c(6,7,i)]
  a = colnames(df3)[3]
  colnames(df3)[3] = "trait"
  aov = lm(trait ~ Genotype*WR_ID, data = df3)
  hsd = aov(aov)              
  resid = residuals(object = hsd)
  shap = shapiro.test(x=resid)
  SW[1,2] = shap$statistic 
  SW[2,2] = shap$p.value 
  colnames(SW)[2] = a 
  SW = as.data.frame(SW)
  if (shap$p.value > 0.05){  
    print(shapiro.test(x=resid)$p.value)
    aov = lm(trait ~ Genotype*WR_ID, data = df3)
    write.table(a, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS23.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    Z = ""
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS23.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS23.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)  
    hsd_test = HSD.test(hsd, trt = c("WR_ID","Genotype"), console = TRUE)  
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = a 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS24b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS24b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
  #for data that needs normalization
  if (shap$p.value < 0.05){
    print(shapiro.test(x=resid)$p.value)
    norm_col_name = paste0("norm_", a)
    par(mfrow=c(3,1))
    df3$norm_col_name = transformTukey(df3$trait,
                                       start = -10,
                                       end = 10,
                                       int = 0.025,
                                       plotit = TRUE, #can be false 
                                       verbose = FALSE,
                                       quiet = FALSE,
                                       statistic = 1,
                                       returnLambda = FALSE  )
    aov = lm(norm_col_name ~ Genotype*WR_ID, data = df3)
    col = norm_col_name
    write.table(col, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS23.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS23.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS23.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)              
    hsd_test = HSD.test(hsd, trt = c("WR_ID","Genotype"), console = TRUE)
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = col 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS24b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS24b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
}
traits = colnames(df2)[1:5] 
outfile = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS24a.txt"
for (trait in traits) {
  temp_summary = df2 %>%
    group_by(Genotype, WR_ID) %>%
    summarise(
      r = n(),
      Mean = mean(.data[[trait]], na.rm = TRUE),
      SD = sd(.data[[trait]], na.rm = TRUE),
      SE = SD / sqrt(r),
      Min = min(.data[[trait]], na.rm = TRUE),
      Max = max(.data[[trait]], na.rm = TRUE),
      Q25 = quantile(.data[[trait]], 0.25, na.rm = TRUE),
      Q50 = quantile(.data[[trait]], 0.50, na.rm = TRUE),
      Q75 = quantile(.data[[trait]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = WR_ID,
      values_from = c(Mean, SD, r, SE, Min, Max, Q25, Q50, Q75),
      names_glue = "{WR_ID}_{.value}"
    ) %>%
    select(
      Genotype,
      matches("^BRNode1_"),
      matches("^BRNode2_"),
      matches("^BRNode3_"),
      everything()
    )
  
  header = paste0("\n===== Trait: ", trait, " =====\n")
  
  write(header, file = outfile, append = TRUE)
  write.table(temp_summary, file = outfile, append = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
}

#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/FigureS13.pdf", width=6, height=4)
ggplot(df, aes(x = SMURF_cat, y = K_N.mm, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("K (N/mm)")+
  xlab(NULL)+
  ggtitle("A")+
  expand_limits(y = 0)+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(#legend.position = "none",
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
ggplot(df, aes(x = SMURF_cat, y = a0_mm, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("a0 (mm)")+
  xlab(NULL)+
  ggtitle("B")+
  expand_limits(y = 0)+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(#legend.position = "none",
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
ggplot(df, aes(x = SMURF_cat, y = b0_mm, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("b0 (mm)")+
  xlab(NULL)+
  ggtitle("C")+
  expand_limits(y = 0)+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(#legend.position = "none",
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
ggplot(df, aes(x = SMURF_cat, y = I_mm4, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("I (mm4)")+
  xlab(NULL)+
  ggtitle("D")+
  expand_limits(y = 0)+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(#legend.position = "none",
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
ggplot(df, aes(x = SMURF_cat, y = E_GPa, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("E (GPa)")+
  xlab(NULL)+
  ggtitle("E")+
  expand_limits(y = 0)+
  scale_x_discrete(labels = c("low" = "Low", 
                              "average" = "Moderate", 
                              "high" = "High"))+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(#legend.position = "none",
    axis.text.x = element_text(size=8, angle = 60, hjust=1),
    axis.text.y = element_text(size=8),
    axis.text = element_text(size=8),
    axis.title = element_text(face="bold", size=10),
    axis.title.x = element_text(face="bold", size=10),
    axis.title.y = element_text(face="bold", size=10))+
  facet_wrap(~WR_ID, labeller = labeller(WR_ID = id_labels), ncol=1)
#dev.off()

#Set up data frame to fill in 2
SW = matrix(NA,nrow=2,ncol=2)
rownames(SW) = c("W","pvalue")
SW[1,1] = "W"
SW[2,1] = "pvalue"
Z = "" #placeholder
colnames(df)
df2 = df[,c(4:8,11,3)]
head(df2)
colnames(df2)
for (i in 1:5){
  df3 = df2[,c(6,7,i)]
  a = colnames(df3)[3]
  colnames(df3)[3] = "trait"
  aov = lm(trait ~ SMURF_cat*WR_ID, data = df3)
  hsd = aov(aov)              
  resid = residuals(object = hsd)
  shap = shapiro.test(x=resid)
  SW[1,2] = shap$statistic 
  SW[2,2] = shap$p.value 
  colnames(SW)[2] = a 
  SW = as.data.frame(SW)
  if (shap$p.value > 0.05){  
    print(shapiro.test(x=resid)$p.value)
    aov = lm(trait ~ SMURF_cat*WR_ID, data = df3)
    write.table(a, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS25.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    Z = ""
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS25.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS25.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)  
    hsd_test = HSD.test(hsd, trt = c("SMURF_cat"), console = TRUE)  
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = a 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS26b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS26b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    
    hsd_test = HSD.test(hsd, trt = c("WR_ID"), console = TRUE)  
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = a 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS27b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS27b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
  #for data that needs normalization
  if (shap$p.value < 0.05){
    print(shapiro.test(x=resid)$p.value)
    norm_col_name = paste0("norm_", a)
    par(mfrow=c(3,1))
    df3$norm_col_name = transformTukey(df3$trait,
                                       start = -10,
                                       end = 10,
                                       int = 0.025,
                                       plotit = TRUE, #can be false 
                                       verbose = FALSE,
                                       quiet = FALSE,
                                       statistic = 1,
                                       returnLambda = FALSE  )
    aov = lm(norm_col_name ~ SMURF_cat*WR_ID, data = df3)
    col = norm_col_name
    write.table(col, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS25.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    anova = anova(aov)
    anova = as.data.frame(anova)
    write.table(anova, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS25.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS25.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd = aov(aov)              
    hsd_test = HSD.test(hsd, trt = c("SMURF_cat"), console = TRUE)
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = col 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS26b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS26b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
    hsd_test = HSD.test(hsd, trt = c("WR_ID"), console = TRUE)
    groups1 = as.data.frame(hsd_test$groups)
    colnames(groups1)[1] = col 
    write.table(groups1, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS27b.txt", sep = "\t",
                row.names = TRUE, col.names = TRUE, append=TRUE)
    write.table(Z, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS27b.txt", sep = "\t",
                row.names = FALSE, col.names = FALSE, append=TRUE)
  }
}
head(df2)
traits = colnames(df2)[1:5] 
outfile = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS26a.txt"
for (trait in traits) {
  temp_summary = df2 %>%
    group_by(SMURF_cat) %>%
    summarise(
      r = n(),
      Mean = mean(.data[[trait]], na.rm = TRUE),
      SD = sd(.data[[trait]], na.rm = TRUE),
      SE = SD / sqrt(r),
      Min = min(.data[[trait]], na.rm = TRUE),
      Max = max(.data[[trait]], na.rm = TRUE),
      Q25 = quantile(.data[[trait]], 0.25, na.rm = TRUE),
      Q50 = quantile(.data[[trait]], 0.50, na.rm = TRUE),
      Q75 = quantile(.data[[trait]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  header = paste0("\n===== Trait: ", trait, " =====\n")
  
  write(header, file = outfile, append = TRUE)
  write.table(temp_summary, file = outfile, append = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
}
traits = colnames(df2)[1:5] 
outfile = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS27a.txt"
for (trait in traits) {
  temp_summary = df2 %>%
    group_by(WR_ID) %>%
    summarise(
      r = n(),
      Mean = mean(.data[[trait]], na.rm = TRUE),
      SD = sd(.data[[trait]], na.rm = TRUE),
      SE = SD / sqrt(r),
      Min = min(.data[[trait]], na.rm = TRUE),
      Max = max(.data[[trait]], na.rm = TRUE),
      Q25 = quantile(.data[[trait]], 0.25, na.rm = TRUE),
      Q50 = quantile(.data[[trait]], 0.50, na.rm = TRUE),
      Q75 = quantile(.data[[trait]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  header = paste0("\n===== Trait: ", trait, " =====\n")
  
  write(header, file = outfile, append = TRUE)
  write.table(temp_summary, file = outfile, append = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
}

head(df)
df4 = subset(df, WR_ID == "BRNode1")
Figure4A = ggplot(df4, aes(x = SMURF_cat, y = K_N.mm, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("K (N/mm)")+
  xlab("Root System Stiffness Category")+
  scale_x_discrete(labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  expand_limits(y = 0)+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10))
Figure4B = ggplot(df4, aes(x = SMURF_cat, y = I_mm4, fill = SMURF_cat)) +
  geom_boxplot()+
  ylab("I (mm4)")+
  xlab("Root System Stiffness Category")+
  scale_x_discrete(labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  expand_limits(y = 0)+
  scale_fill_manual(values = c("pink4","steelblue2", "burlywood"), labels = c("average" = "Moderate", "high" = "High", "low" = "Low"))+
  labs(fill = "RSS Category")+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10))
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/Figure4.pdf", width=6, height=4)
ggdraw() +
  draw_plot(Figure4A, x = 0, y = 0, width = 0.5, height = 1) +
  draw_plot(Figure4B, x = 0.5, y = 0, width = 0.5, height = 1) +
  draw_plot_label(label = c("A", "B"), 
                  size = 10,
                  x = c(0,0.5), 
                  y = c(1,1))
#dev.off()

#The root system stiffness captures the underground root system mechanics#####
#Model 1
rm(list=ls()) 
par(mfrow=c(1,1))
data = read.csv("RhizovisionData2017_55&Friends.csv", header = TRUE, na.strings = "NA")
head(data)
colnames(data)
data = data[,c(1,23,10,14,22,15,16,21,17,19,18,12,20,13,24,25)]
colnames(data)
df1 = data %>%
  group_by(Genotype) %>%
  summarise(across(1:15, list(mean = ~mean(. , na.rm = TRUE)))
  )
main = as.data.frame(df1)
data = read.csv("RootPhotoData_2019.csv", header = TRUE, na.strings = "NA")
colnames(data)
data = data[,c(1,10,11,12,13,14)]
colnames(data)
df1 = data %>%
  group_by(Genotype) %>%
  summarise(across(1:5, list(mean = ~mean(. , na.rm = TRUE)))
  )
df1 = as.data.frame(df1)
head(df1)
main = merge(main, df1, by = "Genotype")
head(main)

data1 = read.csv("BRW_55Friends.csv", header = TRUE, na.strings = "NA")
head(data1)
data1 = data1 %>%
  group_by(Genotype, Year, Plot.ID, Plant.Number) %>%
  summarise(Brace.Root.Whorls.in.the.Soil = first(Brace.Root.Whorls.in.the.Soil), .groups = "drop")
data = read.csv("InstronData_55&Friends.csv", header = TRUE, na.strings = "NA")
head(data)
colnames(data1)[3] = "PlotNumber"
colnames(data1)[4] = "PlantNumber"
head(data1)
head(data)
data = merge(data, data1, by = c("Genotype","Year","PlotNumber","PlantNumber"))
head(data)
unique(data$Brace.Root.Whorls.in.the.Soil)
str(data)
data$Brace.Root.Whorls.in.the.Soil = as.numeric(data$Brace.Root.Whorls.in.the.Soil)
data$Whorl.or.Internode_Number = as.numeric(data$Whorl.or.Internode_Number)
data = data %>%
  mutate(
    Whorl.or.Internode_Number = as.numeric(as.character(Whorl.or.Internode_Number)),
    WR_ID = case_when(
      # 4 whorls
      Brace.Root.Whorls.in.the.Soil == 4 & Whorl.or.Internode_Number == 4 ~ "BRNode1",
      Brace.Root.Whorls.in.the.Soil == 4 & Whorl.or.Internode_Number == 3 ~ "BRNode2",
      Brace.Root.Whorls.in.the.Soil == 4 & Whorl.or.Internode_Number == 2 ~ "BRNode3",
      Brace.Root.Whorls.in.the.Soil == 4 & Whorl.or.Internode_Number == 1 ~ "BRNode4",
      is.na(Whorl.or.Internode_Number) & Brace.Root.Whorls.in.the.Soil == 4 ~ "BRNode4",
      
      # 3 whorls
      Brace.Root.Whorls.in.the.Soil == 3 & Whorl.or.Internode_Number == 3 ~ "BRNode1",
      Brace.Root.Whorls.in.the.Soil == 3 & Whorl.or.Internode_Number == 2 ~ "BRNode2",
      Brace.Root.Whorls.in.the.Soil == 3 & Whorl.or.Internode_Number == 1 ~ "BRNode3",
      is.na(Whorl.or.Internode_Number) & Brace.Root.Whorls.in.the.Soil == 3 ~ "BRNode3",
      
      # 2 whorls
      Brace.Root.Whorls.in.the.Soil == 2 & Whorl.or.Internode_Number == 2 ~ "BRNode1",
      Brace.Root.Whorls.in.the.Soil == 2 & Whorl.or.Internode_Number == 1 ~ "BRNode2",
      is.na(Whorl.or.Internode_Number) & Brace.Root.Whorls.in.the.Soil == 2 ~ "BRNode2",
      
      # 1 whorl
      Brace.Root.Whorls.in.the.Soil == 1 & Whorl.or.Internode_Number == 1 ~ "BRNode1",
      is.na(Whorl.or.Internode_Number) & Brace.Root.Whorls.in.the.Soil == 1 ~ "BRNode1",
      
      TRUE ~ NA_character_
    )
  )
str(data)
data = data[!is.na(data$WR_ID), ]
colnames(data)
colnames(data)[17] = "K_N.mm"
str(data)
data$a0_mm = (data$Vertical.Specimen.Diameter)/2
data$b0_mm = (data$Horizontal.Specimen.Diameter)/2
data$I_mm4 = (pi/4)*(((data$a0_mm)^3)*(data$b0_mm))
data$E_MPa = data$K_N.mm * ((17.5^3)/(48*(data$I_mm4)))
head(data)
colnames(data)
df = data[,c(1,2,21,17,22:25)]
colnames(df)
unique(df$Year)
df$E_GPa = df$E_MPa/1000
colnames(df)
df = df %>%
  filter(if_all(K_N.mm:E_GPa, ~ . > 0))
df = subset(df, b0_mm < 10)
df = subset(df, WR_ID != "BRNode4")
head(df)
df = df[,c(1,3:9)]
df1 = df %>%
  group_by(Genotype, WR_ID) %>%
  summarise(across(1:6, list(mean = ~mean(. , na.rm = TRUE)))
  )
df1 = as.data.frame(df1)
head(df1)
df1 = df1 %>%
  pivot_wider(
    names_from = WR_ID,
    values_from = K_N.mm_mean:E_GPa_mean
  )
df1 = as.data.frame(df1)
head(df1)
colnames(df1)
df1 = df1[,c(1:4,11:13,17:19)]
main = merge(main, df1, by = "Genotype")
head(main)
colnames(main)
colnames(main) = gsub("_mean", "", colnames(main))
colnames(main) = gsub("_cm3", "", colnames(main))
colnames(main) = gsub("_cm2", "", colnames(main))
colnames(main) = gsub("_cm", "", colnames(main))
colnames(main) = gsub("_deg", "", colnames(main))
colnames(main) = gsub("_cm", "", colnames(main))
colnames(main) = gsub("_N.mm", "", colnames(main))
colnames(main) = gsub("_mm4", "", colnames(main))
colnames(main) = gsub("_mm", "", colnames(main))
colnames(main) = gsub("_GPa", "", colnames(main))
colnames(main)
rm(list = setdiff(ls(), "main"))
data = main
rm(list = setdiff(ls(), "data"))
df2 = read.csv("SMURFClassifications.csv", header = TRUE, na.strings = "NA")
data = merge(data, df2, by = "Genotype")
colnames(data)
data = data[,c(31,1:30)]
colnames(data)
head(data)
plant_data = data[,c(1,3:24,26,27,29,30)]
colnames(plant_data)
plant_data = janitor::clean_names(plant_data)
head(plant_data)
plant_data$smurf_cat = as.factor(plant_data$smurf_cat)
plant_split = initial_split(plant_data, prop = 0.8, strata = "smurf_cat")
plant_train = training(plant_split)
plant_test  = testing(plant_split)
head(plant_train)
rf_recipe = recipe(smurf_cat ~ ., data = plant_train) %>%
  step_string2factor(all_nominal_predictors()) %>%
  step_impute_mean(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors())
plant_folds = vfold_cv(plant_train, v = 5, repeats = 3)
rf_model = rand_forest( 
  trees = 1000,
  mtry = tune(),
  min_n = tune()
) %>%
  set_engine("randomForest") %>%
  set_mode("classification")
rf_wf = workflow() %>% 
  add_recipe(rf_recipe) %>%
  add_model(rf_model)
rf_grid = expand.grid(
  mtry = 3:7,
  min_n = 1:5
)
rf_tune_results = tune_grid( 
  rf_wf,
  resamples = plant_folds,
  grid = rf_grid,
  metrics = metric_set(accuracy)
)
cv_metrics = rf_tune_results %>% 
  collect_metrics() %>%
  arrange(desc(mean))
best_combo = select_best(rf_tune_results, metric = "accuracy") 
best_combo
final_rf = finalize_model( 
  rf_model,
  best_combo
)
final_wf = workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(final_rf)
final_fit = final_wf %>% fit(data = plant_train)
predictions = predict(final_fit, plant_test) %>% 
  bind_cols(plant_test)
predictions
conf_mat(predictions, truth = smurf_cat, estimate = .pred_class)
metrics(predictions, truth = smurf_cat, estimate = .pred_class)
rf_model_obj = final_fit %>% #Variable importance
  extract_fit_parsnip() %>%
  .$fit
var_imp = randomForest::importance(rf_model_obj)
var_imp_df = as.data.frame(var_imp) %>%
  rownames_to_column(var = "variable") %>%
  rename(importance = MeanDecreaseGini) %>%
  arrange(desc(importance))
ggplot(var_imp_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Predictor", y = "Importance")

var_imp_df
var_imp_df$Model = "Model1"
var_imp_df = var_imp_df %>%
  mutate(variable = factor(variable, levels = variable[order(importance)]))

head(var_imp_df)
var_imp_df$DataType = var_imp_df$variable
var_names = levels(var_imp_df$variable)
var_names
var_imp_df = var_imp_df %>%
  mutate(DataType = recode(DataType,
                        "median_diameter" =  "Belowground Architecture", 
                        "spread_width" = "Aboveground Architecture",     
                        "average_diameter" =  "Belowground Architecture",       
                        "root_angle" = "Aboveground Architecture",    
                        "depth" =   "Belowground Architecture",        
                        "i_br_node1" = "Aboveground Mechanics",             
                        "solidity" = "Belowground Architecture",            
                        "e_br_node1" =  "Aboveground Mechanics",          
                        "e_br_node2"  =  "Aboveground Mechanics",          
                        "k_br_node2" = "Aboveground Mechanics",          
                        "root_heightonstalk" ="Aboveground Architecture",      
                        "k_br_node1" =  "Aboveground Mechanics",   
                        "i_br_node2" = "Aboveground Mechanics",           
                        "single_root_width" = "Aboveground Architecture",
                        "perimeter" =  "Belowground Architecture",   
                        "convex_area" = "Belowground Architecture",
                        "stalk_width"  = "Aboveground Architecture",          
                        "max_width" = "Belowground Architecture",         
                        "steep_angle_frequency"  = "Belowground Architecture",          
                        "volume" =  "Belowground Architecture",                        
                        "total_rl" = "Belowground Architecture",                        
                        "lower_root_area" = "Belowground Architecture",                 
                        "max_diameter" = "Belowground Architecture",                    
                        "network_area" = "Belowground Architecture",                   
                        "surface_area" = "Belowground Architecture",                    
                        "medium_angle_frequency" = "Belowground Architecture"))
head(var_imp_df)
var_imp_df = var_imp_df %>%
  mutate(variable = factor(variable, levels = variable))
var_imp_df = var_imp_df %>%
  arrange(importance) %>%
  mutate(y_pos = row_number())
y_labels = setNames(var_imp_df$DataType, var_imp_df$y_pos)
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/Figure5A.pdf", width=5, height=7)
Figure5A = ggplot(var_imp_df, aes(x = Model, y = y_pos, fill = importance)) +
  geom_tile(color = "white", linewidth = 0.2, width = 1.5, height = 0.9) +
  geom_text(aes(label = round(importance, 2)), color = "black", size = 2) +
  scale_fill_gradient2(low = "skyblue1", mid = "white", high = "dodgerblue4",
                       midpoint = 0, na.value = "grey90", name = "Value") +
  scale_y_continuous(
    breaks = var_imp_df$y_pos,
    labels = y_labels
  ) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        panel.grid = element_blank()
  )
Figure5A
#dev.off()

head(var_imp_df)
feature_type_imp = var_imp_df %>%
  group_by(DataType) %>%
  summarise(
    n_traits = n(),
    total_importance = sum(importance),
    mean_importance_per_trait = mean(importance)
  ) %>%
  arrange(desc(mean_importance_per_trait))
feature_type_imp

#Model 2
cat("\014")
rm(list = setdiff(ls(), c("Figure5A")))
par(mfrow=c(1,1))
data = read.csv("RhizovisionData2017_55&Friends.csv", header = TRUE, na.strings = "NA")
head(data)
colnames(data)
data = data[,c(1,23,10,14,22,15,16,21,17,19,12,18,20,13,24,25)]
colnames(data)
df2 = read.csv("SMURFClassifications.csv", header = TRUE, na.strings = "NA")
data = merge(data, df2, by = "Genotype")
colnames(data)
data = data[,c(17,2:16)]
colnames(data)
head(data)
plant_data = data
colnames(plant_data)
plant_data = janitor::clean_names(plant_data)
head(plant_data)
plant_data$smurf_cat = as.factor(plant_data$smurf_cat)
plant_split = initial_split(plant_data, prop = 0.8, strata = "smurf_cat")
plant_train = training(plant_split)
plant_test  = testing(plant_split)
head(plant_train)
rf_recipe = recipe(smurf_cat ~ ., data = plant_train) %>%
  step_string2factor(all_nominal_predictors()) %>%
  step_impute_mean(all_numeric_predictors()) %>%
  step_impute_mode(all_nominal_predictors())
plant_folds = vfold_cv(plant_train, v = 5, repeats = 3)
rf_model = rand_forest( 
  trees = 1000,
  mtry = tune(),
  min_n = tune()
) %>%
  set_engine("randomForest") %>%
  set_mode("classification")
rf_wf = workflow() %>% 
  add_recipe(rf_recipe) %>%
  add_model(rf_model)
rf_grid = expand.grid(
  mtry = 3:7,
  min_n = 1:5
)
rf_tune_results = tune_grid( 
  rf_wf,
  resamples = plant_folds,
  grid = rf_grid,
  metrics = metric_set(accuracy, roc_auc)
)
cv_metrics = rf_tune_results %>% 
  collect_metrics() %>%
  arrange(desc(mean))
best_combo = select_best(rf_tune_results, metric = "accuracy") 
best_combo
final_rf = finalize_model( 
  rf_model,
  best_combo
)
final_wf = workflow() %>%
  add_recipe(rf_recipe) %>%
  add_model(final_rf)
final_fit = final_wf %>% fit(data = plant_train)
predictions = predict(final_fit, plant_test) %>% 
  bind_cols(plant_test)
predictions
conf_mat(predictions, truth = smurf_cat, estimate = .pred_class)
metrics(predictions, truth = smurf_cat, estimate = .pred_class)
rf_model_obj = final_fit %>% #Variable importance
  extract_fit_parsnip() %>%
  .$fit
var_imp = randomForest::importance(rf_model_obj)
var_imp_df = as.data.frame(var_imp) %>%
  rownames_to_column(var = "variable") %>%
  rename(importance = MeanDecreaseGini) %>%
  arrange(desc(importance))
ggplot(var_imp_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(x = "Predictor", y = "Importance")

head(var_imp_df)
unique(var_imp_df$variable)
var_imp_df = var_imp_df %>%
  mutate(variable = recode(variable,
                           "solidity" = "Solidity",
                           "max_diameter_cm" = "Max Diameter",
                           "steep_angle_frequency" = "Steep Angle Freq.",
                           "surface_area_cm2" = "Surface Area",
                           "convex_area_cm2" = "Convex Area",
                           "network_area_cm2" = "Network Area",
                           "lower_root_area_cm2" = "Lower Root Area",
                           "max_width_cm" = "Max Width",          
                           "medium_angle_frequency" = "Medium Angle Freq.",
                           "total_rl_cm" = "Total Root Length",
                           "depth_cm" = "Depth",
                           "volume_cm3" = "Volume",
                           "perimeter_cm" = "Perimeter",
                           "median_diameter_cm" = "Median Diameter",
                           "average_diameter_cm" = "Average Diameter"))
head(var_imp_df)
print(var_imp_df)
trait_description = tibble(
  variable = c(
    # Distribution
    "Total Root Length", "Solidity", "Lower Root Area",
    
    # Shape
    "Width:Depth",
    
    # Extent
    "Depth", "Max Width", "Network Area", "Convex Area",
    
    # Size
    "Perimeter", "Surface Area",
    "Average Diameter", "Median Diameter", "Max Diameter",
    "Volume","Steep Angle Freq.","Medium Angle Freq." 
  ),
  Description = c(
    rep("Distribution", 3),
    "Shape",
    rep("Extent", 4),
    rep("Size", 8)
  )
)
var_imp_df = var_imp_df %>%
  left_join(trait_description, by = "variable")
print(var_imp_df)
Figure5B = ggplot(var_imp_df, aes(x = reorder(variable, importance), y = importance, fill = Description)) +
  geom_col() +
  scale_fill_manual(values = c("skyblue4","dodgerblue3","skyblue1"))+
  geom_text(aes(label = round(importance, 2)), 
            hjust = -0.1,        # small offset to the right
            size = 3) +          # text size
  coord_flip() +
  labs(x = "", y = "Mean Decrease in Gini") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  expand_limits(y = max(var_imp_df$importance) * 1.1)  # make room for text
Figure5B
#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/Figure5.pdf", width=7, height=5)
ggdraw() +
  draw_plot(Figure5A, x = 0, y = 0, width = 0.25, height = 1) +
  draw_plot(Figure5B, x = 0.25, y = 0, width = 0.75, height = 1) +
  draw_plot_label(label = c("A", "B"), 
                  size = 10,
                  x = c(0,0.25), 
                  y = c(1,1))
#dev.off()



#Root hairs positively impact root system stiffness ####
par(mfrow=c(1,1))
cat("\014")
rm(list=ls()) 
data = read.csv("RhizovisionData2017_55&Friends.csv", header = TRUE, na.strings = "NA")
head(data)
unique(data$Genotype)
unique(data$Year)
colnames(data)
data = data[,c(1,16)]
head(data)
data$MaxWidth_mm = data$MaxWidth_cm *10
data = data[,c(1,3)]
df1 = data %>%
  group_by(Genotype) %>%
  summarise(across(1, list(mean = ~mean(. , na.rm = TRUE)))
  )
df1 = as.data.frame(df1)
head(df1)
colnames(df1)[2] = "width_mm"
data = read.csv("SMURFDatabase.csv", header = TRUE, na.strings = "NA")
data = subset(data, Year == "2021")
data = data[,c(1,15)]
head(data)
data = data %>%
  group_by(Genotype) %>%
  summarise(across(1, list(mean = ~mean(. , na.rm = TRUE)))) 
head(data)
colnames(data)[2] = "RSS"
df1 = merge(df1, data, by = "Genotype")
head(df1)
data = read.csv("SMURF_RootHairless.csv", header = TRUE, na.strings = "NA")
head(data)
unique(data$Genotype)
data$Genotype = factor(data$Genotype, levels = c("wt","mt"))
Figure6A = ggplot(data, aes(x = Genotype, y = line_raw_slope_N.m, fill = Genotype)) +
  geom_boxplot()+
  ylab("Root System Stiffness (N/m)")+
  xlab(NULL)+
  scale_y_continuous(limits=c(0,12000), breaks=seq(0,12000,1000))+
  scale_x_discrete(labels = c("wt" = "Wild Type", "mt" = "rth3 mutant"))+
  scale_fill_manual(values = c("seagreen4","lightgreen"))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8, angle = 60, hjust=1),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10))
Figure6A
attach(data)
head(data)
lm_x = lm(line_raw_slope_N.m ~ Genotype)
anova(lm_x)
lm_x_aov=aov(lm_x) 
X = HSD.test(lm_x_aov, trt = c("Genotype"), console = TRUE)
Means = X$means
#write.csv(Means, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS30a.csv", row.names = TRUE)
Groups = X$groups
Groups
#write.csv(Groups, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS30b.csv", row.names = TRUE)
par(mfrow=c(2,2))
plot(lm_x)
par(mfrow=c(2,1))
plot(data$line_raw_slope_N.m)
hist(data$line_raw_slope_N.m)
resid = residuals(object = lm_x_aov)
shapiro.test(x=resid)
detach(data)


head(data)
colnames(data)
data = data[,c(4,15)]
data = data %>%
  group_by(Genotype) %>%
  summarise(across(1, list(mean = ~mean(. , na.rm = TRUE)))) 
colnames(data)[2] = "RSS"
df = read.csv("rth_RhizovisionData_2025.csv", header = TRUE, na.strings = "NA")
head(df)
unique(df$Geno)
head(df)
df$Geno = factor(df$Gen, levels = c("WT","Mut"))
Figure6B = ggplot(df, aes(x = Geno, y = Width_mm, fill = Geno)) +
  geom_boxplot()+
  ylab("Max Width (mm)")+
  xlab(NULL)+
  scale_y_continuous(limits=c(0,325), breaks=seq(0,325,25))+
  scale_x_discrete(labels = c("WT" = "Wild Type", "Mut" = "rth3 mutant"))+
  scale_fill_manual(values = c("seagreen4","lightgreen"))+
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(size=8, angle = 60, hjust=1),
        axis.text.y = element_text(size=8),
        axis.text = element_text(size=8),
        axis.title = element_text(face="bold", size=10),
        axis.title.x = element_text(face="bold", size=10),
        axis.title.y = element_text(face="bold", size=10))
Figure6B

attach(df)
head(df)
lm_x = lm(Width_mm ~ Geno)
anova(lm_x)
lm_x_aov=aov(lm_x) 
X = HSD.test(lm_x_aov, trt = c("Geno"), console = TRUE)
Means = X$means
Means
#write.csv(Means, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS32a.csv", row.names = TRUE)
Groups = X$groups
#write.csv(Groups, file = "/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/TableS32b.csv", row.names = TRUE)
Groups
par(mfrow=c(2,2))
plot(lm_x)
par(mfrow=c(2,1))
plot(df$Width_mm)
hist(df$Width_mm)
resid = residuals(object = lm_x_aov)
shapiro.test(x=resid)
detach(df)

head(df)
df = df[,c(3,4)]
df = df %>%
  group_by(Geno) %>%
  summarise(across(1, list(mean = ~mean(. , na.rm = TRUE)))) 
head(df)
colnames(df)[1] = "Genotype"
colnames(df)[2] = "Width_mm"
head(data)
data = data %>%
  mutate(Genotype = recode(as.character(Genotype),
                           "wt" = "WT",
                           "mt" = "Mut"))
head(data)
head(df)
data = merge(data, df, by = "Genotype")
head(data)
head(df1)

df1_plot = df1 %>% mutate(source = "original")
head(df1_plot)
data2 = data %>% rename(width_mm = Width_mm) %>% mutate(source = "extra")
head(data2)
fit = lm(RSS ~ width_mm, data = df1_plot)
fit
s = summary(fit)
s
intercept = coef(fit)[1]
intercept
slope = coef(fit)[2]
slope
r2 = s$r.squared
r2
plot_df = bind_rows(df1_plot, data2)
head(plot_df) #full dataframe, original + extra

get_val = function(df, geno, col){
  r = df[df$Genotype == geno, , drop = FALSE]
  if(nrow(r) == 0) stop(paste0("Genotype '", geno, "' not found in dataframe"))
  as.numeric(r[[col]][1])
}

b73_w = get_val(df1, "B73", "width_mm")
b73_r = get_val(df1, "B73", "RSS")

wt_w  = get_val(data2, "WT", "width_mm")
wt_r  = get_val(data2, "WT", "RSS")

mut_w = get_val(data2, "Mut", "width_mm")
mut_r = get_val(data2, "Mut", "RSS")


scale_w = wt_w / b73_w
scale_r = wt_r / b73_r


mut_to_B73_w = mut_w / scale_w
mut_to_B73_r = mut_r / scale_r
mut_to_B73_w
mut_to_B73_r

WT_to_B73_w = wt_w / scale_w
WT_to_B73_r = wt_r / scale_r
WT_to_B73_w
WT_to_B73_r


mutE_row = data.frame(
  Genotype = "Mut_Pred",
  width_mm = mut_to_B73_w,
  RSS      = mut_to_B73_r,
  source   = "predicted")
plot_df = rbind(plot_df, mutE_row)
WTE_row = data.frame(
  Genotype = "WT_Pred",
  width_mm = WT_to_B73_w,
  RSS      = WT_to_B73_r,
  source   = "predicted")
plot_df = rbind(plot_df, WTE_row)
data3 = rbind(WTE_row, mutE_row)

head(df1_plot) #geno averages; 43&friends
head(data2) #rth3 averages 
head(data3) #rth3 predicted
Figure6C = ggplot(plot_df, aes(x = width_mm, y = RSS)) +
  geom_point(data = df1_plot, color = "black", size = 2) +
  geom_smooth(data = df1_plot, method = "lm", se = TRUE, color = "black", fullrange = TRUE) +
  geom_point(data = data3, aes(color = Genotype, shape = Genotype), size = 3) +
  geom_text(data = data3, aes(color = Genotype, label = Genotype),
            vjust = -1, fontface = "bold",show.legend = FALSE ) +
  scale_color_manual(values = c("Mut_Pred" = "lightgreen", "WT_Pred" = "seagreen4")) +
  scale_shape_manual(values = c("Mut_Pred" = 17, "WT_Pred" = 15)) +
  labs(x = "Width (mm)", y = "Root System Stiffness (N/m)",
       color = "Genotype", shape = "Genotype") +
  scale_y_continuous(limits = c(0, 10000), breaks = seq(0, 10000, 1000)) +
  scale_x_continuous(limits = c(0, 400), breaks = seq(0, 400, 25)) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8, angle = 60, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 10),
    axis.title.y = element_text(face = "bold", size = 10)
  )
Figure6C


#pdf(file="/Users/ashley/Desktop/Hostetler_Pierce_et_al_2025/Figures/Figure6.pdf", width=7, height=5)
ggdraw() +
  draw_plot(Figure6A, x = 0, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(Figure6B, x = 0.5, y = 0.5, width = 0.5, height = 0.5) +
  draw_plot(Figure6C, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B","C"), 
                  size = 10,
                  x = c(0,0.5,0), 
                  y = c(1,1,0.5))
#dev.off()
