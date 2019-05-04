# Comparison between MaxQuant and Comet/PAW pipelines (v2)

## Comparison between edgeR statistical testing and two-sample t-test

#### Notebooks:

[Comet/PAW](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_PAW.html) - Updated notebook with analysis using Comet/PAW results

> Thorough workup of Comet/PAW data with QC checks. Testing is done with edgeR.

[Comet/PAW edgeR vs t-test](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_PAW_t-test.html) - Comparison of edgeR to t-test

> Comparison of edgeR to a two-sample t-test data workup.

[Comet/PAW edgeR vs limma](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_PAW_limma.html) - Comparison of edgeR to limma

> Comparison of edgeR to limma. limma with its normal data distribution may be more appropriate for intensity data.

[Comet/PAW edgeR vs limma-voom](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_PAW_limma-voom.html) - Comparison of edgeR to limma with voom

> The group that developed limma and edgeR has a recommended way to combine the best of both packages for analyzing RNA-seq data. We follow that outline here and see how it compares.

[MaxQuant](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_MQ.html) - updated notebook with analysis using MaxQuant results

> Thorough workup of MaxQuant data with QC checks. Testing is done with edgeR.

---

[old notebook](https://pwilmart.github.io/TMT_analysis_examples/KUR1502_MQ_PAW.html) - older notebook with first analysis

> Older MaxQuant version was used. Both PAW and MaxQuant results are in the same notebook.

## Data is from this publication:
> Huan, J., Hornick, N.I., Goloviznina, N.A., Kamimae-Lanning, A.N., David, L.L., Wilmarth, P.A., Mori, T., Chevillet, J.R., Narla, A., Roberts Jr, C.T. and Loriaux, M.M., 2015. Coordinate regulation of residual bone marrow function by paracrine trafficking of AML exosomes. Leukemia, 29(12), p.2285.

This was a mouse bone marrow cell culture experiment with controls (n=3) and leukemia exosome-dosed cells (n=4). The data was collected on a Thermo Fusion tribrid instrument using the SPS MS3 method. Two free data analysis options were explored: the [Comet/PAW pipeline](https://github.com/pwilmart/PAW_pipeline.git) or [MaxQuant](https://www.maxquant.org). A canonical UniProt mouse reference protein database was used in both analyses. More details on the search parameters are in the respective notebooks. There are two folders of results files: one for PAW and one for MaxQuant. Each folder has separate Jupyter notebooks for the TMT analysis. The notebooks are very similar in layout and are designed to be opened in side-by-side browser windows for a head-to-head comparison. Some of the markdown cells (mostly in the MaxQuant notebook) differ to illustrate important points. The notebooks are pre-rendered as HTML files and have links at the top of this file.

A widely-used genomics statistical package (edgeR) is used in the analysis of results from both MaxQuant (`KUR1502_MQ`) and PAW (`KUR1502_PAW`). Comparisons of edgeR to other statistical analyses are in separate notebooks (t-test:`KUR1502_PAW_t-test`, limma:`KUR1502_PAW_limma`, and limma-voom:`KUR1502_PAW_limma-voom`). There are many papers that have demonstrated the inadequacy of t-tests when replicate numbers are small in these types of omics results. limma and edgeR, with designs for these types of datasets, outperform t-tests in sensitivity and accuracy.

The notebooks in this 2019 update have nicer R code cell contents compared to the 2018 notebook. There is better use of functions and some simplifications in the workflow. The notebooks and R scripts in the repository can serve as starting templates for other analyses.

Further information:

- [Jupyter Notebooks](https://jupyter.org/)
- [R kernel for notebooks](https://irkernel.github.io/)
- [modern R](https://r4ds.had.co.nz/)
- [RStudio](https://www.rstudio.com/)
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html)
- [limma](http://bioconductor.org/packages/release/bioc/html/limma.html)


updated May 4, 2019
