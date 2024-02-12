
---

***This repository contains instructions and code to reproduce the analyses presented in***

# Systematic Analysis of Biological Processes Reveals Gene Co-expression Modules Driving Pathway Dysregulation in Alzheimer’s Disease
---

Here, we employ single-nucleus RNA-sequencing data from various brain regions to perform a comprehensive analysis of Alzheimer's disease (AD). We goes beyond traditional differential gene expression (DEG) analysis by integrating pathway activity analysis with weighted gene co-expression networks. This approach allows for a detailed mapping of gene interconnectivity and the identification of region- and cell-type-specific drivers of AD-related biological processes. The findings reveal significant heterogeneity in modular organization and functional disruptions in both neuronal and glial cells, with extended involvement of astrocytes and microglia in various processes beyond neuroinflammation. Our study also observes limited overlap of DEGs within perturbed pathways, indicating that DEGs might not fully represent the disease's complexity. Further, we show distinct dynamics of hub DEGs in neuronal versus glial modules, suggesting a greater impact of DEGs on neurons than on glial cells. These results underscore the importance of a systems-oriented approach, combining pathway enrichment and co-expression methods, for a comprehensive understanding of AD-related biological processes.

#### Paper
[Access the Paper Here](...)

<p align="center">
  <img src="Figure 1.png" width="800" title="hover text">
</p>

#### snRNAseq Data Availability

For processing `raw counts matrices, associated metadata, fully-annotated and QC-ed data`:

- [Leng et al., 2021](https://doi.org/10.1038/s41593-020-00764-7) files available at AD Knowledge Portal, Synapse ID: [syn21788402](https://www.synapse.org/#!Synapse:syn21788402).
- [Gabitto et al., 2023](https://www.biorxiv.org/content/10.1101/2023.05.08.539485v2) data at the [SEA-AD Documentation Page](https://portal.brain-map.org/explore/seattle-alzheimers-disease/seattle-alzheimers-disease-brain-cell-atlas-download?edit&language=en).

#### Other Data Availability

- Subject-level metadata with pathology grouping available [here](...).

### Reproduce Analyses and Plots

Follow these instructions to replicate analyses and plots from the paper.

#### Clone the Repository:

```bash
git clone https://github.com/TemiLeke/systematic_ad_analysis.git
```

##### Set Up Environment:

```
conda env create -f systematic_ad_analysis.yml
conda activate systematic_ad_analysis
```

##### Steps to Reproduce Analyses

1. Download data into /data directory. See table below for data files and their sources.
2. For each study/brain region, create `/data/{study_name}/`. Save .h5ad files and CSV metadata accordingly.
2. Ensure sample-level pathological status is stratified as `no-`, `early-`, and `late-pathology`
3. Run notebooks. 
- `pathway_analysis.ipynb` ***(note study-specific differences)***
- `DEG_analysis.ipynb` ***(note study-specific differences)***
- `pathway_meta_analysis.ipynb` ***(note study-specific differences)***

##### Data File Overview

| Data File                                                     | Description / Origin                                                                                                                                                                                                    |       
|---------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| pathway_databases/GO_Biological_Process_2018.txt              | from mayaan lab  [here](https://maayanlab.cloud/Enrichr/#libraries)                                                                                                                                                     |
| pathway_databases/MONDO_0004975-associated-diseases.tsv       | from Open Targets Platform [here](https://platform.opentargets.org/disease/MONDO_0004975/associations)                                                                                                                                                      |
| pathway_databases/AD_genes.csv       | from Harmonizome (Mayaanlab) [here](https://maayanlab.cloud/Harmonizome/gene_set/Alzheimer+Disease/dbGAP+Gene-Trait+Association s)                                                                                                                                                      |
| pathway_databases/KEGG_2019_Human.txt                         | from mayaan lab [here](https://maayanlab.cloud/Enrichr/#libraries)                                                                                                                                                      |
| pathway_databases/go_terms.obo                         | from Gene Ontology [here](https://geneontology.org/docs/download-ontology/)                                                                                                                                                      |


### Citation

If you use the methods and scripts provided here, please cite:

`Adeoye, T, Shah, S.I., Ullah, G. Systematic Analysis of Biological Processes Reveals Gene Co-expression Modules Driving Pathway Dysregulation in Alzheimer’s Disease. (2023).`


### Contact

For questions, inconsistencies, help requests, or ideas for future collaborations, reach out at tadeoye@usf.edu

If you have any questions, notice any inconsistencies, need help, or would like to brainstorm future collaborations and ideas, please don't hesitate to reach out: `tadeoye@usf.edu``
