#'##############################################################################
#' Collapse gene expression data into pathways with Hipathia
#'##############################################################################

# Load libraries and data ####
library(hipathia)
library(Biobase)
library(SummarizedExperiment)

load("data/genexpr.Rdata")

# Compute pathways activation ####
## Map gene IDs
gexp_SE <- makeSummarizedExperimentFromExpressionSet(genexpr)
translated <- translate_data(gexp_SE, "hsa")

## Normalize data
exp_data <- normalize_data(translated)

## Apply hipathia
pathways <- load_pathways(species = "hsa")
gexp_path <- hipathia(exp_data, pathways, decompose = FALSE, verbose = TRUE)
path_vals <- get_paths_data(gexp_path)
save(gexp_path, file = "results/gene_expression/gexpr.hipathia.full.Rdata")
save(path_vals, file = "results/gene_expression/gexpr.hipathia.pathways.Rdata")
