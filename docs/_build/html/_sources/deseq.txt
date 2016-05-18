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

	gene	transcript	exper_rep1	exper_rep2	exper_rep3	control_rep1	control_rep2	control_rep3

	DQ714826	uc007tfl.1	0.00	0.00	0.00	0.00	0.00	0.00

	DQ551521	uc008bml.1	0.00	0.00	0.00	0.00	0.00	0.00

	AK028549	uc011wpi.1	2.00	1.29	0.00	0.00	0.00	1.40


Where the gene column represent the gene name, the transcript column represents the transcript(s) name (comma separated for multiple),
and the rest of the columns are the raw counts for your samples.

DESeq2 Workflow
===============

DESeq2 performs multiple steps in order to analyze the data you've provided for it.
The first step is to indicate the condition that each column (experiment) in the table represent.
You can group multiple samples into one condition column.
DESeq2 will compute the probability that a gene is differentially expressed (DE) for ALL genes in the table. It outputs
both a nominal and a multiple hypothesis corrected p-value (padj) using a negative binomial distribution.

DEBrowser
=========

DEBrowser utilizes `Shiny <http://shiny.rstudio.com/>`_, a R based application development tool that creates a wonderful interactive user interface (UI)
combinded with all of the computing prowess of R.  After the user has selected the data to analyze and has used the shiny
UI to run DESeq2, the results are then input to DEBrowser.  DEBrowser manipulates your results in a way that allows for
interactive plotting by which changing padj or fold change limits also changes the displayed graph(s).
For more details about these plots and tables, please visit our quickstart guide for some helpful tutorials.

References
==========

1. Love MI, Huber W and Anders S (2014). Moderated estimation of fold change and
    dispersion for RNA-seq data with DESeq2.  Genome Biology, 15, pp. 550.
    http://doi.org/10.1186/s13059-014-0550-8.
