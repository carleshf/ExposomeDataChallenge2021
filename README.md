The following is a project proposal for the [Exposome Data Challenge 2021](https://www.isglobal.org/-/exposome-data-analysis-challenge) organised by the [Barcelona Institute of Global Heath](https://www.isglobal.org).

---

__Timeline__ (from the oficial webpage):

  * February 18: Exposome dataset available [here](https://github.com/isglobal-exposomeHub/ExposomeDataChallenge2021/blob/main/README.md).
  * March 4: _We submitted this proposal_.
  * March 8: Proposal submission deadline.
  * March 15: Proposal selection results.
  * March-April:Time for the selected participants to apply their method on the Exposome dataset provided.
  * April 28: Selected participants apply their methods on the Exposome dataset provided, submit the results and codes on the event Github page.
  * April 28-30: Showcase of the results and discussion. Selected participants present their methods and results on the Exposome dataset provided. The results and codes will be submitted on the event Github page.

---

# Dimensionality reduction in exposome and omic data to uncover molecular pathologic mechanisms 

__Authors__:
 * Carlos Ruiz-Arenas ([linkedin](https://es.linkedin.com/in/cruizarenas))
 * Carles Hernandez-Ferrer ([linkedin](https://es.linkedin.com/in/carleshf))

__Keyword__: _exposome_, _transcriptome_, _metabolome_, _dimensionality reduction_, _pathway association analyses_, _integration analyses_. 

__Background__: The exposome represents the sum of all environmental exposures over a lifetime and, at the intersection of public health and toxicology, is a key tool to look forward in the understanding of the disease development machinery, completing the knowledge obtained from genetics (1). The traditional genome-wide association studies inspired a range of methods to investigate the underlying links in exposome with multiple health outcomes (2). The exposome-wide association approach consists of a series of linear and logistic regression models to assess the association between single exposures and outcomes (3). Current exposome studies have included omic measurements (aka, trascriptomics, proteomics, and metabolomics) in the study design. Omic measurements have provided new insights for disease physiology and have been proved as valuable biomarkers for disease (4,5). However, omic measurements greatly increase the data dimensionality, requiring enormous datasets to achieve statistical power to detect associations using single exposure-omic measurements tests. Besides, omic features have not typically a clear biological meaning (6), so exposure-omic feature associations reaching statistical significance are commonly difficult to interpret and translate to biological knowledge. 

__Goal__: We aim to define a new approach to elucidate the relationships between the exposome and a series of molecular signatures (aka, trascriptomics, and metabolomics) in a generalized fashion. To improve the interpretability of the results, we will reduce the data dimensionality by mapping the exposures and molecular signatures to biochemical and biological pathways. 

__Methods__: We will collapse the exposome into biochemical pathways using an exposure-to-exposure mapping from the [Human Metabolome Database](https://hmdb.ca) (10) and the [Toxin and Toxin-Target Database](www.t3db.ca) (11). The same methodology will be applied to both urine and serum metabolome datasets. We will consider merging urine and serum datasets, after collapsing them to pathway level. We will collapse gene expression in pathways activation using the [`hipatia`](https://bioconductor.org/packages/release/bioc/html/hipathia.html) (4) algorithm. Then, the state of the art exposome-wide association analysis will be performed but at a pathway-to-pathway level. Finally, [`MOFA+`](https://bioconductor.org/packages/release/bioc/html/MOFA2.html) (7,9) and other integrative methodologies (5,8,12–15) may be considered to corroborate the previous results. 

__Expected results__: The first set of results comprehends the relation between each feature (exposome and omic data) with the pathway(s) it belongs to (results of dimensionality reduction) and the cross-association between the exposome and the molecular signatures (results from the exposome-wide association study at pathway level). If considered, the second set of results, from the integration analysis, provide corroboration of the relationship exposure-to-molecular signature in an N-to-N consideration. 

__Challenges__: 1) Combined effects of exposures, thanks to dimensionality reduction on the exposome side; 2) Using omics data to improve inference on the link between exposome and health, due to the exposome-wide association study at pathway level between the reduced exposome and reduced gene expression and metabolome; and 3) Multi-omics analysis, due to the methods for data integration that will include the exposome, the gene expression, and the metabolome. 
