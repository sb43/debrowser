#' getQCPanel, conditional panel for QC plots 
#'
#'
#' @note \code{getQCSection}
#' @return the panel for QC plots;
#'
#' @examples
#'     x <- getQCPanel()
#'
#' @export
#'

getQCPanel <- function() {
    a <- list( conditionalPanel( (condition <- "input.qcplot=='all2all' ||
                                input.qcplot=='heatmap' ||
                                input.qcplot=='pca'"),
            column(2,
                sliderInput("width", "width",
                        min = 100, max = 2000, step = 10, value = 700)),
            column(2,
                sliderInput("height", "height",
                        min = 100, max = 2000, step = 10, value = 500)),
        conditionalPanel( (condition <- "input.qcplot=='all2all'"),
            column(2, sliderInput("cex", "corr font size",
                                  min = 0.1, max = 10,
                                  step = 0.1, value = 2))),
        conditionalPanel( (condition <- "input.qcplot=='heatmap'"),
                column(3, selectInput("clustering_method", "Clustering Method:",
                    choices <- c("complete", "ward.D2", "single", "average",
                        "mcquitty", "median", "centroid"))),
                column(3, selectInput("distance_method", "Distance Method:",
                    choices <- c("cor", "euclidean", "maximum", "manhattan",
                        "canberra", "binary", "minkowski"))),
            column(1, actionButton("startQCPlot", "Submit"))),
            column(1, downloadButton("downloadPlot", "")),
            uiOutput("plotarea")))
}

#'Left menu for QC plots
#'
#' @note \code{getLeftMenu}
#' @return returns the left menu according to the selected tab;
#'
#' @examples
#'     x <- getLeftMenu()
#'
#' @export
#'
getLeftMenu <- function() {
    a <- list( conditionalPanel( (condition <- "input.methodtabs=='panel1'"),
            wellPanel(radioButtons("mainplot", paste("Main Plots:", sep = ""),
            c(Scatter = "scatter", VolcanoPlot = "volcano",
                MAPlot = "maplot")))),
        conditionalPanel( (condition <- "input.methodtabs=='panel2'"),
            wellPanel(radioButtons("qcplot",
                paste("QC Plots:", sep = ""),
            c(All2All = "all2all", Heatmap = "heatmap", PCA = "pca")))),
        conditionalPanel( (condition <- "input.methodtabs=='panel3'"),
            wellPanel(radioButtons("goplot", paste("Go Plots:", sep = ""),
            c(enrichGO = "enrichGO", enrichKEGG = "enrichKEGG",
        Disease = "disease", compareClusters = "compare")))),
        tags$small("Note: Please don't forget to choose appropriate
            dataset to visualize it in the QC plots."))
}

#' getQCPlots, for quality checks 
#'
#'
#' @note \code{getQCPlots}
#' @return the panel for QC plots;
#' @param dataset, the dataset to use
#' @param datasetname, name of the dataset
#' @param qcplot, type of plot to add
#' @param metadata, coupled samples and conditions
#' @param clustering_method, clustering method used
#' @param distance_method, distance method used
#' @param cex, font size
#' @examples
#'     x <- getQCPlots(mtcars)
#'
#' @export
#'
getQCPlots <- function(dataset, datasetname = "Up", qcplot = "heatmap",
                        metadata = NULL,
                        clustering_method = "complete",
                        distance_method = "cor", cex = 2) {
    a <- NULL
    if (nrow(dataset) > 0) {
        if (qcplot == "all2all") {
            a <- all2all(dataset, cex)
        } else if (qcplot == "heatmap") {
            a <- runHeatmap(dataset, title = paste("Dataset:", datasetname),
                clustering_method = clustering_method,
                distance_method = distance_method)
        } else if (qcplot == "pca") {
            colnames(metadata) <- c("samples", "conditions")
            pca_data <- run_pca(getNormalizedMatrix(dataset))
            a <- plot_pca(pca_data$PCs, explained = pca_data$explained,
                metadata = metadata, color = "samples",
                size = 5, shape = "conditions",
                factors = c("samples", "conditions"))
        }
    }
    a
}
