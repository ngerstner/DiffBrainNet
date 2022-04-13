# Brain region specific gene regulation after treatment with glucocorticoids

![<img src="./06_Shiny/www/mousejavi_reversebrain.png" width="150"/>](./06_Shiny/www/mousejavi_reversebrain.png) TODO: General text about project, similar to paper abstract.

## 1. Dataset

### 1.1 Experimental setup

TODO: add data overview from slides

### 1.2 Data availability

TODO: add here link to Shiny App where data can be downloaded

## 2. Data analysis

### 2.1 Differential gene expression analysis

Differential gene expression (DE) analysis between mice treated with a vehicle and mice treated with dexamethasone was performed with DESeq2 [cite?] in R. Analysis scripts for the DE analysis can be found in the folder [01_DiffExp](01_DiffExp/).

### 2.2 Differential network analysis

Differential expression networks were calculated with KiMONo [cite?] and DiffGRN [cite?]. As a prior for the network analysis, we used [Funcoup](https://funcoup5.scilifelab.se/search/) which includes functional associations between genes and their respective proteins from various evidence types. Scripts that were used to calculate the expression networks on control and DEX treated samples respectively, can be found in the folder [02_CoExp_Kimono](02_CoExp_Kimono/). 

During the post-processing of the networks inferred by KiMONo, we filtered out interactions with beta or r-squared values below a certain threshold and inferred a *differential gene expression network* for each brain region in our dataset using the DiffGRN [cite?] approach. Scripts used for the post-processing can be found in the folder [03_CoExp_Analysis](03_CoExp_Analysis/).

## 3. Manuscript plots

The plots in our manuscript [link?] were generated with the scripts in the folder (04_PlotsManuscript)[04_PlotsManuscript/].


