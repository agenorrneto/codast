# Introduction

CODAST (**CO**nstraint and **D**omain-based r**A**re variant a**S**sociation **T**est) is framework for rare variant association tests that integrates constraint and three dimensional data.
In this repository, users can find a brief (and still under construction and review) step-by-step guide to perform this framework.


# Step by Step

## Overview

CODAST can be summarized in two main steps:

1. Data mining and,
2. Association tests.

Most part of the coded provided in this repository is part of Step 1, since this framework relies heavily in the processing of secondary data from biological databases.
In the folder `scripts`, users can find scripts used for automatization of both steps. The folder `notebooks` is intended for storing Python notebook for steps not yet automatized or for statistical analysis
and image generation.

### Data Mining
CODAST needs protein domain data coordinates for data annotation and integration. Here, based on the importance of Protein-Protein Interactions (PPIs) for cell biology, we suggest users to first retrieve manually data from [STRING](https://string-db.org/) and then using [InterProScan](https://www.ebi.ac.uk/interpro/search/sequence/) for domain annotation.

With domain annotation available, users can use `scripts/gcoord_translator.py` as a translator from protein coordinates to genomic coordinates (GRCh38). You can use this functionality with:

```
python3 gcoord_translator.py -ip [InterProScan  annotation result] \
                             -db [Database of choice to retrieve the annotations (e.g. Pfam)]
                             -outdir [Directory for output]    
```

This script's output file can used as an input to the other scripts. We propose integration with the following metrics:

| Type of metric  |  |
| ------------- | ------------- |
| Intraspecies constraint  | LOEUF, z-score missense, pLI, CCRs, PER  |
| Interspecies constraint  | GERP++, ConSurf |
| 3D data  | RSA, secondary structure |

In `notebooks/ccrs_integration.ipynb`, you can find the codes for integration and analysis of CCR (constrained coding regions) and gnomAD contraint metrics and the domain genomic coordinates. We plan to automatize this step.

`scripts/gerp_metrics_integrator.py` is intended for integration with GERP++ data. You can use it with:

```
python3 gcoord_translator.py -gc [Protein domain genomic coordinates] \
                             -out_dir [Directory for output]    
```

`scripts/3d_metrics_integrator.py` is intended for integration with three-dimensional data. You can use it with:
```
python3 gcoord_translator.py -gc [Protein domain genomic coordinates] \
                             -out_dir [Directory for output]    
```
### Association Test

Our framework uses [rvtest]() to perform rare variant association test. The biological unit proposed in the tests is the domain. 
The data generated in the previous test will be used to classify the domains and help in the interpretation of results. We are still studying how to make this classification and integrate in a better way this data into the tests.



## 💭 Feedback and Contributing

