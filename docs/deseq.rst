***************
DESeq/DEBrowser
***************

This guide contains a breif discription of DESeq2 used within the DEBrowser


Introduction
============

Differential gene expression analysis has become an increasingly popular tool
in determining and viewing up and/or down experssed genes between two sets of
samples.  The goal of Differential gene expression analysis is to find genes
or transcripts whose difference in expression, when accounting for the
variance within condition, is higher than expected by chance.  `DESeq2
<https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ is an R
package available via Bioconductor and is designed to normalize count data
from high-throughput sequencing assays such as RNA-Seq and test for
differential expression (Love et al. 2014).  For more information on the
DESeq2 algorithm, you can visit `this website <https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf>`_  With multiple parameters such as
padjust values, log fold changes, and plot styles, altering plots
created with your DE data can be a hassle as well as time consuming.  The
Differential Expression Browser uses DESeq2 coupled with shiny to produce
real-time changes within your plot queries and allows for interactive browsing
of your DESeq results. In addition to DESeq analysis, DEBrowser also offers
a variety of other plots and analysis tools to help visualize your data
even further.


Getting Started
===============

In order to conduct differential expression analysis, we first need data to analyze.  In order to call DESeq2,
we're going to need gene quantifications and expected counts for those genes.
To obtain these quantifications, we typically use `RSEM <http://deweylab.github.io/RSEM/>`_,
however there are other ways to obtain this data.

The TSV files used to describe the quantification counts are similar to this:

IE:

=====  =====  =====  =====  =====  =====
gene   trans   exp1   exp2  cont1  cont2
=====  =====  =====  =====  =====  =====
DQ714  uc007   0.00   0.00   0.00   0.00
DQ554  uc008   0.00   0.00   0.00   0.00
AK028  uc011   2.00   1.29   0.00   0.00
=====  =====  =====  =====  =====  =====

Where the gene column represent the gene name, the transcript column represents the transcript(s) name (comma separated for multiple),
and the rest of the columns are the raw counts for your samples.

DESeq2
=========

For the details please check the user guide.
`DESeq2 userguide <https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf>`_

DESeq2 performs multiple steps in order to analyze the data you've provided for it.
The first step is to indicate the condition that each column (experiment) in the table represent.
You can group multiple samples into one condition column.
DESeq2 will compute the probability that a gene is differentially expressed (DE) for ALL genes in the table. It outputs
both a nominal and a multiple hypothesis corrected p-value (padj) using a negative binomial distribution.

Un-normalized counts
====================
DESeq2 rquires count data as input obtained from RNA-Seq or another high-thorughput sequencing experiment in the form of matrix values. Here we convert un-integer values to integer to be able to run DESeq2. The matrix values should be un-normalized, since DESeq2 model internally corrects for library size. So, transformed or normalized values such as counts scaled by library size should not be used as input. Please use edgeR or limma for normalized counts.

Used parameters for DESeq2
==========================
  - fitType:
     either "parametric", "local", or "mean" for the type 
     of fitting of dispersions to the mean intensity. 
     See estimateDispersions for description.

  - betaPrior: 
     whether or not to put a zero-mean normal prior
     on the non-intercept coefficients See nbinomWaldTest for 
     description of the calculation of the beta prior. By default, 
     the beta prior is used only for the Wald test, but can also be 
     specified for the likelihood ratio test.

  - testType: 
     either "Wald" or "LRT", which will then use either 
     Wald significance tests (defined by nbinomWaldTest), or the 
     likelihood ratio test on the difference in deviance between a 
     full and reduced model formula (defined by nbinomLRT)

  - rowsum.filter: 
     regions/genes/isoforms with total count (across all samples) below this value will be filtered out

EdgeR
========
For the details please check the user guide.
`EdgeR userguide <https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf>`_.

Used parameters for EdgeR
=========================

  - Normalization:
     Calculate normalization factors to scale the raw 
     library sizes. Values can be "TMM","RLE","upperquartile","none".

  - Dispersion:
     either a numeric vector of dispersions or a character 
     string indicating that dispersions should be taken from the data 
     object. If a numeric vector, then can be either of length one or 
     of length equal to the number of genes. Allowable character 
     values are "common", "trended", "tagwise" or "auto". 
     Default behavior ("auto" is to use most complex dispersions 
     found in data object.

  - testType: 
     exactTest or glmLRT. exactTest: Computes p-values for differential 
     abundance for each gene between two digital libraries, conditioning 
     on the total count for each gene. The counts in each group as a 
     proportion of the whole are assumed to follow a binomial distribution. 
     glmLRT: Fit a negative binomial generalized log-linear model to the read 
     counts for each gene. Conduct genewise statistical tests for a given 
     coefficient or coefficient contrast.
     
  - rowsum.filter: 
     regions/genes/isoforms with total count (across all samples) below this value will be filtered out
  
Limma
========
For the details please check the user guide.
`Limma userguide <https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf>`_.

Limma is a package to analyse of microarray or RNA-Seq data. If data is normalized with spike-in or any other scaling, tranforamtion or normalization method, Limma can be ideal. In that case, prefer limma rather than DESeq2 or EdgeR.

Used parameters for Limma
=========================

  - Normalization: 
     Calculate normalization factors to scale the raw library sizes. Values can be "TMM","RLE","upperquartile","none".
  
  - Fit Type: 
     fitting method; "ls" for least squares or "robust" for robust regression
  
  - Norm. Bet. Arrays: 
     Normalization Between Arrays; Normalizes expression intensities so that the    
     intensities or log-ratios have similar distributions across a set of arrays.

  - rowsum.filter:
      regions/genes/isoforms with total count (across all samples) below this value will be filtered out

DEBrowser
=========

DEBrowser utilizes `Shiny <http://shiny.rstudio.com/>`_, a R based application development tool that creates a wonderful interactive user interface (UI)
combinded with all of the computing prowess of R.  After the user has selected the data to analyze and has used the shiny
UI to run DESeq2, the results are then input to DEBrowser.  DEBrowser manipulates your results in a way that allows for
interactive plotting by which changing padj or fold change limits also changes the displayed graph(s).
For more details about these plots and tables, please visit our quickstart guide for some helpful tutorials.

For comparisons against other popular data visualization tools, see the table below.

.. image:: debrowser_pics/comparison_table.png
	:align: center

References
==========

1. Love MI, Huber W and Anders S (2014). Moderated estimation of fold change and
    dispersion for RNA-seq data with DESeq2.  Genome Biology, 15, pp. 550.
    http://doi.org/10.1186/s13059-014-0550-8.
2. Robinson, MD, and Smyth, GK (2008). Small sample estimation of negative binomial dispersion,
    with applications to SAGE data. Biostatistics 9, 321–332.
3. Ritchie, ME, Phipson, B, Wu, D, Hu, Y, Law, CW, Shi, W, and Smyth, GK (2015).
    limma powers differential expression analyses for RNA-sequencing and microarray studies.
    Nucleic Acids Research 43(7), e47.
