#'##############################################################################
#' Run analysis of exposures PCs vs gene expression
#'##############################################################################

# Load libraries and data ####
library(SummarizedExperiment)
library(limma)
library(tidyverse)
library(RColorBrewer )

load("results/gene_expression/gexpr.hipathia.pathways.Rdata")
load("results/exposures/exposomePC.RData")
load("data/exposome.RData")

## Add covariates to colData
rownames(covariates) <- covariates$ID
covariates <- covariates[ , -1]
colData(path_vals) <- cbind(colData(path_vals), covariates[colnames(path_vals), c("h_cohort", "h_edumc_None")])

rownames(phenotype) <- phenotype$ID
phenotype <- phenotype[ , -1]
colData(path_vals) <- cbind(colData(path_vals), phenotype[colnames(path_vals), "hs_zbmi_who", drop = FALSE ])

## Add exposures to colData
colData(path_vals) <- cbind(colData(path_vals), exposomePC[paste0('sample_', colnames(path_vals)), ])


exps <- colnames(exposomePC)

runExpo <- function(expname){
  
  model <- model.matrix(formula(paste("~", expname, "+ e3_sex + h_cohort + h_edumc_None + hs_zbmi_who + age_sample_years")),
                        colData(path_vals))
  lmF <- lmFit(assay(path_vals), model)
  lmFe <- eBayes(lmF)
  tab <- topTable(lmFe, number = Inf, coef = 2)
  tab$exposure <- expname
  tab$feat.ID <- rownames(tab)
  tab
}
expoRes <- lapply(exps, runExpo)
expoResdf <- Reduce(rbind, expoRes) %>%
  tibble() %>%
  arrange(P.Value) %>%
  mutate(adj.P.Value.all = p.adjust(P.Value, method = "BH")) %>%
  left_join(as.data.frame(rowData(path_vals)) %>% select(-decomposed), by = "feat.ID") %>%
  left_join(mutate(exposomePC_fam, exposure = var) %>% select(-var), by = "exposure")

## Heatmap FC ####
a <- expoResdf %>%
  select(logFC, feat.name, exposure) %>%
  spread(feat.name, logFC) 

logmat <- data.matrix(a[, -1])
rownames(logmat) <- a$exposure

cols <- brewer.pal(n = 9, name = "Paired")
names(cols) <- unique(expoResdf$family)
colPlot <- cols[exposomePC_fam[rownames(logmat), "family"]]

heatmap(t(logmat), ColSideColors = colPlot, main = "FC plot")

## Heatmap pvalues ####
b <- expoResdf %>%
  mutate(logP = -log10(P.Value)) %>%
  select(logP, feat.name, exposure) %>%
  spread(feat.name, logP) 

pmat <- data.matrix(b[, -1])
rownames(pmat) <- b$exposure

colPlotp <- cols[exposomePC_fam[rownames(pmat), "family"]]

heatmap(t(pmat), ColSideColors = colPlotp, main = "-log10 plot")

expoRes.sel <- expoResdf %>%
  filter(P.Value < 1e-3) %>%
  arrange(family)

## Heatmap FC filtered ####
a <- expoResdf %>%
  filter(P.Value < 1e-3) %>%
  select(logFC, feat.name, exposure) %>%
  spread(feat.name, logFC) 

logmat <- data.matrix(a[, -1])
rownames(logmat) <- a$exposure
logmat[is.na(logmat)] <- 0

cols <- brewer.pal(n = 9, name = "Paired")
names(cols) <- unique(exposomePC_fam[rownames(logmat), "family"])
colPlot <- cols[exposomePC_fam[rownames(logmat), "family"]]

heatmap(t(logmat), ColSideColors = colPlot, main = "FC plot")

## Heatmap pvalues ####
b <- expoResdf %>%
  filter(P.Value < 1e-3) %>%
  mutate(logP = -log10(P.Value)) %>%
  select(logP, feat.name, exposure) %>%
  spread(feat.name, logP) 

pmat <- data.matrix(b[, -1])
rownames(pmat) <- b$exposure
pmat[is.na(pmat)] <- 2.9

colPlotp <- cols[exposomePC_fam[rownames(pmat), "family"]]

heatmap(t(pmat), ColSideColors = colPlotp, main = "-log10 plot")



expoRes.sel %>%
  group_by(exposure) %>%
  summarize(n = n())
tab <- tab[, colSums(tab) > 0]

data.frame(path = as.vector(assay(path_vals["P-hsa04015-95", ])), 
                 colData(path_vals)[,  c("pbde_pc1", "h_cohort", "e3_sex")]) %>%
  ggplot(aes(x = pbde_pc1, y = path, color = h_cohort, shape = e3_sex)) +
  geom_point() +
  geom_smooth(aes(x = pbde_pc1, y = path), method = lm, se = FALSE, inherit.aes = FALSE) +
  theme_bw()

data.frame(path = as.vector(assay(path_vals["P-hsa04010-65", ])), 
           colData(path_vals)[,  c("airpollution_pc1", "h_cohort", "e3_sex")]) %>%
  ggplot(aes(x = airpollution_pc1, y = path, color = h_cohort, shape = e3_sex)) +
  geom_point() +
  geom_smooth(aes(x = airpollution_pc1, y = path), method = lm, se = FALSE, inherit.aes = FALSE) +
  theme_bw()

data.frame(path = as.vector(assay(path_vals["P-hsa05205-451", ])), 
           colData(path_vals)[,  c("pbde_pc1", "h_cohort", "e3_sex")]) %>%
  ggplot(aes(x = pbde_pc1, y = path, color = h_cohort, shape = e3_sex)) +
  geom_point() +
  geom_smooth(aes(x = pbde_pc1, y = path), method = lm, se = FALSE, inherit.aes = FALSE) +
  theme_bw()


data.frame(path = as.vector(assay(path_vals["P-hsa04015-95", ])), 
           colData(path_vals)[,  c("hs_pbde153_cadj_Log2", "h_cohort", "e3_sex")]) %>%
  ggplot(aes(x = hs_pbde153_cadj_Log2, y = path, color = h_cohort, shape = e3_sex)) +
  geom_point() +
  geom_smooth(aes(x = hs_pbde153_cadj_Log2, y = path), method = lm, se = FALSE, inherit.aes = FALSE) +
  theme_bw()
