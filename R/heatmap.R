#' runHeatmap
#'
#' Creates a heatmap based on the user selected parameters within shiny.
#'
#' @param data, a matrixthat includes expression values
#' @param title, title of the heatmap
#' @param dend, dendogram
#' @param names, a flag to show the rownames
#' @param clustering_method = c('complete', 'ward.D2', 'single', 'average',
#' 'mcquitty', 'median' , 'centroid')
#' @param distance_method = c('cor','euclidean', 'maximum', 'manhattan',
#' 'canberra', 'binary' ,'minkowski')
#' @return heatmap.2 plot
#'
#' @examples
#'     x <- runHeatmap(mtcars)
#'
#' @export
#' @import gplots
#' @import RColorBrewer
#'
runHeatmap <- function(data, title="Title", dend = "both",
    names = FALSE,
    clustering_method = c("ward.D2", "complete", "single",
        "average", "mcquitty", "median", "centroid"),
    distance_method = c("euclidean", "cor", "maximum",
        "manhattan", "canberra", "binary", "minkowski")) {
    ld <- log2(data + 0.1)
    cldt <- scale(t(ld), center = TRUE, scale = TRUE)
    cld <- t(cldt)
    hclust2 <- function(x, ...) hclust(x, method = clustering_method)
    dist2 <- function(x, ...) {
        if (distance_method != "cor") {
            return(dist(x, method = distance_method))
        } else {
            return(as.dist(1 - cor(t(x))))
        }
    }
    m <- heatmap.2(cld, Rowv = TRUE, main = title, dendrogram = dend,
        Colv = TRUE, col = redblue(256), labRow = names,
        distfun = dist2, hclustfun = hclust2, density.info = "none",
        trace = "none", margins = c(10, 10))
    m
}

#' cellInfo
#'
#' hover info in heatmap
#'
#' @param x, data
#' @return data
#' @export
#'
#' @examples
#'     x <- cellInfo()
#'
cellInfo <- function(x = NULL) {
    if(is.null(x)) return(NULL)
    sprintf("%s</br>%s", x$Genes, x$Samples)
}

#' getIntHeatmap
#'
#' getIntHeatmap
#'
#' @param heatData, heatData
#' @param inputQCPlot, input params for QC
#' @param linked brush object for selection
#' @return plot
#' @export
#'
#' @examples
#'     load(system.file("extdata", "demo", "demodata.Rda", package = "debrowser"))
#'     varied <- getMostVariedList(demodata, c(2:7), 50, 100)[,c(2:7)]
#'     heatdat<- runHeatmap(varied)
#'     lbheat <- linked_brush(keys = rownames(t(heatdat$carpet)))
#'     getIntHeatmap(t(heatdat$carpet), lbheat)
#'
getIntHeatmap <- function(heatdat = NULL, count = NULL, 
    init_data = NULL, lbheat = NULL) {
    if (is.null(heatdat)) return (NULL)
    graphheight = "auto"
    if (count > 50)
        graphheight <- count + 200
    a <- heatdat %>%
        ggvis(~Samples, ~Genes, fill = ~Values) %>%
        layer_points( size := 200, size.hover := 200,
            fillOpacity := 0.0, fillOpacity.hover := 0.0,
            fill.brush := "red", opacity := 0.0,
            stroke = ~Samples, key := ~ID) %>%
        layer_rects(width = band(), height = band(), strokeWidth := 0) %>%
        lbheat$input() %>%
        add_tooltip(cellInfo,"hover") %>%
        scale_nominal("x", padding = 0, points=FALSE) %>%
        scale_nominal("y", padding = 0, points=FALSE,  reverse = TRUE) %>%
        add_axis("x", properties=axis_props(labels=list(angle=270)),
            tick_padding = 50, title_offset=120) %>%
        set_options(width = "auto", height = graphheight)  %>%
        scale_numeric("fill", range = c("darkblue", "white")) %>%
        hide_legend("stroke")
    if (count > 50)
        a <- a %>% hide_axis("y")
    else
        a <- a %>% add_axis("y", title_offset=80)
    a %>% bind_shiny(paste0("heatmapplot"))
}
