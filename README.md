# A pan-cancer genome-wide analysis reveals tumour dependencies by induction of nonsense-mediated decay

Zhiyuan Hu1,2,3, Christopher Yau3,4\* and Ahmed Ashour Ahmed1,2\*

1 Weatherall Institute of Molecular Medicine, University of Oxford, Oxford, OX3 9DU, UK  

2 Nuffield Department of Obstetrics and Gynaecology, University of Oxford, Oxford, OX3 9DU, UK

3 Wellcome Trust Centre for Human Genetics, University of Oxford, Oxford, OX3 7BN, UK

4 Centre for Computational Biology, Institute of Cancer and Genomic Sciences, University of Birmingham, Birmingham, B15 2TT, UK

*Corresponding authors

## Background

Herein the repo contains the raw data and coding to reproduce the results in the manusript:
Z. Hu, C. Yau and A. Ahmed (2017) *A pan-cancer genome-wide analysis reveals tumour dependencies by induction of non-sense mediated decay* (Under review).

While some mutations can introduce premature termination codons (PTCs) into genes, nonsense mediated decay (NMD) can detect these PTCs and then eliminate the abnormal transcripts. PTCs and NMD play important roles in genetic diseases and cancers.

Herein, three rules are used to predict whether a PTC-generating mutation is NMD-elicit or NMD-escape:

1. PTC is more than 50-54bp upstream of the last-exon-exon junction.

2. targeted gene is not intronless.

3. PTC is more than 200bp downstream of the start codon.

Using these rules, we can predict whether a called mutation will elicit NMD on the mRNA from the mutated gene, i.e. the NMD-elicit mutations.

## Related R package masonmd

To predict the NMD-elicit mutations in your own dataset, please try our R package [masonmd] (https://github.com/ZYBunnyHu/masonmd)


In R or Rstudio, use the following codes to install the masonmd package directly from Github:

```{R}
install.packages(“devtools”)
devtools::install_github("ZYBunnyHu/masonmd")
```

## Explanation of files


#### Raw data

The raw data is in the `data/` folder. Before start, combine the splitted RNA-seq data by putting the following commands into the terminal.

```{bash}
cd NMD-paper
cd data
cat PANCAN_HiSeqV2a*.zip > PANCAN_HiSeqV2.zip
```

#### R scripts

In total, eight R scripts are included in the `R` subdirectory. You can run them by the order from `01_PANCAN_data.R`.

1. `00_PANCAN_function.R` contains the functions that are used in the prediction or analysis.

2. `01_PANCAN_data.R` is about loading data and preprocessing.

3. `02_PANCAN_classify.R` is about prediction of NMD-elicit mutations.

4. `06_PANCAN_figures_manuscript.R` generates the figures in the manuscript

5. `07_PANCAN_supplementary_figures.R` generates the supplementary figures.


## Citation

A manuscript detailing our work has been submitted for publication and a bioRxiv submission is in preparation:

Z. Hu, C. Yau and A. Ahmed (2017) *A pan-cancer genome-wide analysis reveals tumour dependencies by induction of non-sense mediated decay* (Under review).
