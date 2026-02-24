####stats barplots
#collapse GHG on plot level 
library(dplyr)
library(readr)
library(ggpubr)
library(ggplot2)
library(gridExtra)

library(ggpattern)

d <- read_csv("metadata231124.csv")

library(dplyr)

d_core <- d %>%
  filter(core == "yes") %>%
  select(site,
         m_graminoid,
         m_mean_qPCR_ITS,
         rel_ab_ascomycota,
         rel_ab_basidiomycota,
         m_CO2fluxdark,
         CO2fluxlight) %>%
  mutate(site = factor(site))
#graminoids
ggplot(d_core, aes(site, m_graminoid)) +
  geom_boxplot() +
  theme_bw()

gram_aov <- aov(m_graminoid ~ site, data = d_core)

# Normality of residuals
shapiro.test(residuals(gram_aov))

# Homogeneity of variances
bartlett.test(m_graminoid ~ site, data = d_core)

###is not noramlly distributed and is not homogeneous! thus kruskal walis
kruskal.test(m_graminoid ~ site, data = d_core)

# Kruskal-Wallis rank sum test
# 
# data:  m_graminoid by site
# Kruskal-Wallis chi-squared = 30.719, df = 2, p-value = 2.135e-07

pairwise.wilcox.test(
  d_core$m_graminoid,
  d_core$site,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  d_core$m_graminoid and d_core$site 
# 
# ADT     BJB    
# BJB 4.6e-05 -      
#   LOD 4.6e-05 4.6e-05
d_core %>%
  group_by(site) %>%
  summarise(
    median_gram = median(m_graminoid, na.rm = TRUE),
    IQR_gram = IQR(m_graminoid, na.rm = TRUE)
  )
# # A tibble: 3 × 3
# site  median_gram IQR_gram
# <fct>       <dbl>    <dbl>
#   1 ADT           3.5      3  
# 2 BJB          25       27.5
# 3 LOD          85       12.5


#co2dark
ggplot(d_core, aes(site, m_CO2fluxdark)) +
  geom_boxplot() +
  theme_bw()

codark_aov <- aov(m_CO2fluxdark ~ site, data = d_core)

# Normality of residuals
shapiro.test(residuals(codark_aov))

# Homogeneity of variances
bartlett.test(m_CO2fluxdark ~ site, data = d_core)

###is not noramlly distributed and is not homogeneous! thus kruskal walis
kruskal.test(m_CO2fluxdark ~ site, data = d_core)

# Kruskal-Wallis rank sum test
# 
# data:  m_CO2fluxdark by site
# Kruskal-Wallis chi-squared = 17.622, df = 2, p-value = 0.0001491
pairwise.wilcox.test(
  d_core$m_CO2fluxdark,
  d_core$site,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  d_core$m_CO2fluxdark and d_core$site 
# 
# ADT     BJB    
# BJB 0.00052 -      
#   LOD 0.00052 0.72059
# 
# P value adjustment method: BH 
d_core %>%
  group_by(site) %>%
  summarise(
    median_codark = median(m_CO2fluxdark, na.rm = TRUE),
    IQR_codark = IQR(m_CO2fluxdark, na.rm = TRUE)
  )
# # A tibble: 3 × 3
# site  median_codark IQR_codark
# <fct>         <dbl>      <dbl>
#   1 ADT         0.00178   0.000726
# 2 BJB         0.00377   0.00115 
# 3 LOD         0.00351   0.00442 

#co2light
ggplot(d_core, aes(site, CO2fluxlight)) +
  geom_boxplot() +
  theme_bw()

col_aov <- aov(CO2fluxlight ~ site, data = d_core)

# Normality of residuals
shapiro.test(residuals(col_aov))

# Homogeneity of variances
bartlett.test(CO2fluxlight ~ site, data = d_core)

###is not noramlly distributed and is not homogeneous! thus kruskal walis
kruskal.test(CO2fluxlight ~ site, data = d_core)
# 
# Kruskal-Wallis rank sum test
# 
# data:  CO2fluxlight by site
# Kruskal-Wallis chi-squared = 21.089, df = 2, p-value = 2.634e-05
pairwise.wilcox.test(
  d_core$CO2fluxlight,
  d_core$site,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  d_core$CO2fluxlight and d_core$site 
# 
# ADT     BJB    
# BJB 0.00052 -      
#   LOD 0.00052 0.00557
# 
# P value adjustment method: BH 
d_core %>%
  group_by(site) %>%
  summarise(
    median_codark = median(CO2fluxlight, na.rm = TRUE),
    IQR_codark = IQR(CO2fluxlight, na.rm = TRUE)
  )
# A tibble: 3 × 3
# site  median_codark IQR_codark
# <fct>         <dbl>      <dbl>
#   1 ADT        0.000995   0.000917
# 2 BJB        0.00403    0.00110 
# 3 LOD        0.00273    0.00154 

#itsabundance
ggplot(d_core, aes(site, m_mean_qPCR_ITS)) +
  geom_boxplot() +
  theme_bw()

qI_aov <- aov(m_mean_qPCR_ITS ~ site, data = d_core)

# Normality of residuals
shapiro.test(residuals(qI_aov))

# Homogeneity of variances
bartlett.test(m_mean_qPCR_ITS ~ site, data = d_core)

###is not noramlly distributed and is not homogeneous! thus kruskal walis
kruskal.test(m_mean_qPCR_ITS ~ site, data = d_core)
# 
# Kruskal-Wallis rank sum test
# 
# data:  m_mean_qPCR_ITS by site
# Kruskal-Wallis chi-squared = 11.946, df = 2, p-value = 0.002547
pairwise.wilcox.test(
  d_core$m_mean_qPCR_ITS,
  d_core$site,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  d_core$m_mean_qPCR_ITS and d_core$site 
# 
# ADT    BJB   
# BJB 0.0503 -     
#   LOD 0.0051 0.0362
# 
# P value adjustment method: BH 
d_core %>%
  group_by(site) %>%
  summarise(
    median_codark = median(m_mean_qPCR_ITS, na.rm = TRUE),
    IQR_codark = IQR(m_mean_qPCR_ITS, na.rm = TRUE)
  )
# # A tibble: 3 × 3
# site  median_codark IQR_codark
# <fct>         <dbl>      <dbl>
#   1 ADT      1420216899 661431723.
# 2 BJB      1957805191 964240207.
# 3 LOD      2980002020 945971074.

#ascomycota
ggplot(d_core, aes(site, rel_ab_ascomycota)) +
  geom_boxplot() +
  theme_bw()

Asco_aov <- aov(rel_ab_ascomycota ~ site, data = d_core)

# Normality of residuals
shapiro.test(residuals(Asco_aov))

# Homogeneity of variances
bartlett.test(rel_ab_ascomycota ~ site, data = d_core)

###is not noramlly distributed and is not homogeneous! thus kruskal walis
kruskal.test(rel_ab_ascomycota ~ site, data = d_core)
# 
# Kruskal-Wallis rank sum test
# 
# data:  rel_ab_ascomycota by site
# Kruskal-Wallis chi-squared = 18.988, df = 2, p-value = 7.532e-05
pairwise.wilcox.test(
  d_core$rel_ab_ascomycota,
  d_core$site,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  d_core$rel_ab_ascomycota and d_core$site 
# 
# ADT     BJB    
# BJB 0.063   -      
#   LOD 6.2e-05 6.2e-05
# 
# P value adjustment method: BH  
d_core %>%
  group_by(site) %>%
  summarise(
    median_codark = median(rel_ab_ascomycota, na.rm = TRUE),
    IQR_codark = IQR(rel_ab_ascomycota, na.rm = TRUE)
  )
# # A tibble: 3 × 3
# site  median_codark IQR_codark
# <fct>         <dbl>      <dbl>
#   1 ADT           0.177     0.505 
# 2 BJB           0.626     0.169 
# 3 LOD           0.861     0.0897

#basisdiomycota
ggplot(d_core, aes(site, rel_ab_basidiomycota)) +
  geom_boxplot() +
  theme_bw()

baso_aov <- aov(rel_ab_basidiomycota ~ site, data = d_core)

# Normality of residuals
shapiro.test(residuals(baso_aov))

# Homogeneity of variances
bartlett.test(rel_ab_basidiomycota ~ site, data = d_core)

###is not noramlly distributed and is not homogeneous! thus kruskal walis
kruskal.test(rel_ab_basidiomycota ~ site, data = d_core)
# 
# Kruskal-Wallis rank sum test
# 
# data:  rel_ab_basidiomycota by site
# Kruskal-Wallis chi-squared = 16.511, df = 2, p-value = 0.0002598
pairwise.wilcox.test(
  d_core$rel_ab_basidiomycota,
  d_core$site,
  p.adjust.method = "BH"
)
# Pairwise comparisons using Wilcoxon rank sum exact test 
# 
# data:  d_core$rel_ab_basidiomycota and d_core$site 
# 
# ADT     BJB    
# BJB 1.00000 -      
#   LOD 0.00012 0.00012
# 
# P value adjustment method: BH 
d_core %>%
  group_by(site) %>%
  summarise(
    median_codark = median(rel_ab_basidiomycota, na.rm = TRUE),
    IQR_codark = IQR(rel_ab_basidiomycota, na.rm = TRUE)
  )
# # A tibble: 3 × 3
# site  median_codark IQR_codark
# <fct>         <dbl>      <dbl>
#   1 ADT          0.212      0.286 
# 2 BJB          0.285      0.245 
# 3 LOD          0.0421     0.0452











ggplot(d_core, aes(site, m_graminoid)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  theme_bw()

d_core <- d_core %>%
  mutate(
    log_ITS = log10(m_mean_qPCR_ITS + 1),
    asin_Asco = asin(sqrt(rel_ab_ascomycota))
  )

cor.test(d_core$m_graminoid, d_core$m_mean_qPCR_ITS, method = "spearman")
cor.test(d_core$m_graminoid, d_core$rel_ab_ascomycota, method = "spearman")
cor.test(d_core$m_mean_qPCR_ITS, d_core$rel_ab_ascomycota, method = "spearman")

library(lme4)
library(lmerTest)
mod_ITS <- lmer(
  log_ITS ~ m_graminoid + (1 | site),
  data = d_core
)

summary(mod_ITS)

mod_Asco <- lmer(
  asin_Asco ~ m_graminoid + (1 | site),
  data = d_core
)

summary(mod_Asco)

mod_ITS_Asco <- lmer(
  asin_Asco ~ log_ITS + (1 | site),
  data = d_core
)

summary(mod_ITS_Asco)

par(mfrow = c(2, 2))
plot(mod_ITS)

library(sjPlot)

plot_model(mod_ITS, type = "pred", terms = "m_graminoid")







#####
#ok plot boxplot with grasses?
d <- read_csv("metadata231124.csv")

dm <- d %>%
  group_by(site) %>%
  summarise(
    mean_CO2 = mean(m_CO2fluxdark, na.rm = TRUE),
    sd_CO2 = sd(m_CO2fluxdark, na.rm = TRUE),
  )


ca <- ggplot(dm, aes(x = site, y = mean_CO2, fill=site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + # Side by side bars for CO2
  geom_errorbar(aes(ymin = mean_CO2 - sd_CO2, ymax = mean_CO2 + sd_CO2), position = position_dodge(width = 0.7), width = 0.2) + # Error bars
  labs(title = "Mean CO2 Flux",
       x = "Site",
       y = "Mean CO2 Flux [mol m-3 s-1]") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() + 
  scale_fill_manual(values = c("#ff7f0e","blue","#2ca02c")) # Customize colors for dark and light
ca

#"#ff7f0e","orange","#d62728", "#FEAFAB","#2ca02c", "light green"



dg <- d %>%
  filter(core == "yes") %>%
  group_by(site) %>%
  summarise(
    mean_gram = mean(m_graminoid, na.rm = TRUE),
    sd_gram = sd(m_graminoid, na.rm = TRUE),
  )

gra <- ggplot(dg, aes(x = site, y = mean_gram, fill=site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + # Side by side bars for CO2
  geom_errorbar(aes(ymin = mean_gram - sd_gram, ymax = mean_gram + sd_gram), position = position_dodge(width = 0.7), width = 0.2) + # Error bars
  labs(title = "Mean Graminoid cover",
       x = "Site",
       y = "% Gramioid Cover") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() + 
  scale_fill_manual(values = c("#ff7f0e","blue","#2ca02c")) # Customize colors for dark and light
gra

dI <- d %>%
  filter(core == "yes") %>%
  group_by(site) %>%
  summarise(
    mean_sha = mean(m_mean_qPCR_ITS, na.rm = TRUE),
    sd_sha = sd(m_mean_qPCR_ITS, na.rm = TRUE),
  )

ITSq <- ggplot(dI, aes(x = site, y = mean_sha, fill=site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + # Side by side bars for CO2
  geom_errorbar(aes(ymin = mean_sha - sd_sha, ymax = mean_sha + sd_sha), position = position_dodge(width = 0.7), width = 0.2) + # Error bars
  labs(title = "Mean ITS abundance",
       x = "Site",
       y = "ITS copy numbers") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() + 
  scale_fill_manual(values = c("#ff7f0e","blue","#2ca02c")) # Customize colors for dark and light
ITSq
dp <- d %>%
  filter(core == "yes") %>%
  group_by(site) %>%
  summarise(
    mean_sha = mean(m_mean_qPCR_16s, na.rm = TRUE),
    sd_sha = sd(m_mean_qPCR_16s, na.rm = TRUE),
  )

pq <- ggplot(dp, aes(x = site, y = mean_sha, fill=site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + # Side by side bars for CO2
  geom_errorbar(aes(ymin = mean_sha - sd_sha, ymax = mean_sha + sd_sha), position = position_dodge(width = 0.7), width = 0.2) + # Error bars
  labs(title = "Mean 16S abundance",
       x = "Site",
       y = "16S copy numbers") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() + 
  scale_fill_manual(values = c("#ff7f0e","blue","#2ca02c")) # Customize colors for dark and light
pq

dA <- d %>%
  filter(core == "yes") %>%
  group_by(site) %>%
  summarise(
    mean_sha = mean(rel_ab_ascomycota, na.rm = TRUE),
    sd_sha = sd(rel_ab_ascomycota, na.rm = TRUE),
  )
pA <- ggplot(dA, aes(x = site, y = mean_sha, fill=site)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + # Side by side bars for CO2
  geom_errorbar(aes(ymin = mean_sha - sd_sha, ymax = mean_sha + sd_sha), position = position_dodge(width = 0.7), width = 0.2) + # Error bars
  labs(title = "Mean relative abundance Ascomycota",
       x = "Site",
       y = "Relative Abundance Ascomycota") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme_bw() + 
  scale_fill_manual(values = c("#ff7f0e","blue","#2ca02c")) # Customize colors for dark and light
pA

