# EpiCon
## Introduction
EpiCon is a computational method, proposing an epigenetic controllability score to quantify the ability of TFs to open up closed chromatins for single cells and cell types, taking sc-multiome data as inputs. EpiCon first constructs a TF-RE regulatory network by integrating the single cell multiome input data and motif binding affinity using a correlation-based approach. The epi-controllability score for the single cell for a given TF is then defined as the weighted sum of regulatory strength of TF-RE pairs over REs where the weight is the cell specific chromatin accessibility of REs. For cell type-specific scores, we incorporate motif enrichment analysis, removing non-enriched motifs before network construction.
<div style="text-align: right">
  ![Image](Fig1_small.png)
</div>
EpiCon is validated using experimental datasets including chromatin immunoprecipitation sequencing (ChIP-seq) and TF knockdown data. EpiCon distinguished TFs even when they belong to the same family and share identical motifs. EpiCon discovers experimentally validated driver regulators. Our approach is broadly applicable to any sc-multiome data.
## Install the packages
git 
