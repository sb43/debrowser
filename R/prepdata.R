#' getSamples
#'
#' Gathers the sample names to be used within DEBrowser.
#'
#' @param cnames, names of the  samples
#' @param index, starting column in a tab separated file
#' @return choices
#' @export
#'
#' @examples
#'     x <- getSamples()
#'
getSamples <- function (cnames = NULL, index = 2) { 
    m <- NULL
    if (!is.null(cnames)) {
        cn <- cnames[index:length(cnames)]
        m <- as.list(NULL)
        for (i in seq(cn)) {
            m[i] <- cn[i]
        }
    }
    m
}

#' prepDataContainer
#'
#' Prepares the data container that stores values used within DESeq.
#'
#' @param data, loaded dataset
#' @param counter, the number of comparisons
#' @param input, input parameters
#' @param session, session var
#' @return data
#' @export
#'
#' @examples
#'     x <- prepDataContainer()
#'
prepDataContainer <- function(data = NULL, counter=NULL, 
    input = NULL, session=NULL) {
    if (is.null(input$goButton)) return(NULL)
    if (input$goButton[1]==0) return(NULL)
    inputconds <- reactiveValues(fittype=NULL, conds = list())
    inputconds <- eventReactive(input$goButton, {
    m <- c()
    updateTabsetPanel(session, "methodtabs", selected = "panel1")
    hide(selector = "#methodtabs li a[data-value=panel0]")
    shiny::validate(need(input$condition1, "Condition1 has to be selected"),
        need(input$condition2, "Condition2 has to be selected"))
    m$conds <- list()
    for (cnt in seq(1:(2*counter)))
    {
        m$conds[cnt] <- list(input[[paste0("condition",cnt)]])
    }
    shinyjs::disable("resetsamples")
    m$fittype <- input$fittype
    shinyjs::disable("goButton")
    m
    })
    if (is.null(input$condition1)) return(NULL)
    dclist<-list()
    for (i in seq(1:counter))
    {
        conds <- c(rep(paste0("Cond", 2*i-1), 
        length(inputconds()$conds[[2*i-1]])), 
        rep(paste0("Cond", 2*i), length(inputconds()$conds[[2*i]])))
        cols <- c(paste(inputconds()$conds[[2*i-1]]), 
        paste(inputconds()$conds[[2*i]]))
        m<-prepDESeqOutput(data, cols, conds, inputconds(), i)
        m<-list(conds = conds, cols = cols, init_data=m)
        dclist[[i]] <- m
    }
    togglePanels(1, c(1:10), session)
    return(dclist)
}

#' getMean
#'
#' Gathers the mean for selected condition.
#'
#' @param norm_data, loaded dataset
#' @param de_res, de results
#' @param inputconds, input parameters
#' @param colnum, colnum
#' @return data
#' @export
#'
#' @examples
#'     x <- getMean()
#'
getMean<-function(norm_data = NULL, de_res = NULL, 
    inputconds = NULL, colnum = NULL) {
    if (is.null(norm_data)) return (NULL)
    mean_cond<-NULL
    if (length(inputconds$conds[[colnum]]) > 1)
        mean_cond <-list(rowMeans( norm_data[ rownames( de_res ),
            paste( inputconds$conds[[colnum]] )] ))
    else
        mean_cond <-list(norm_data[ rownames( de_res ),
            paste( inputconds$conds[[colnum]] )])
    mean_cond
}

#' prepDESeqOutput
#'
#' Prepares the output data from DESeq to be used within
#' DEBrowser
#'
#' @param data, loaded dataset
#' @param cols, columns
#' @param conds, conds
#' @param inputconds, inputconds
#' @param i, selected comparison number
#' @return data
#' @export
#'
#' @examples
#'     x <- prepDESeqOutput()
#'
prepDESeqOutput <- function(data = NULL, cols = NULL, 
    conds = NULL, inputconds=NULL, i=NULL) {
    if (is.null(data)) return (NULL)
    if (length(cols) == length(conds))
        de_res <- runDESeq(data, cols, conds, inputconds$fittype,
            non_expressed_cutoff = 10)
    de_res <- data.frame(de_res)
    norm_data <- getNormalizedMatrix(data[, cols])
    mean_cond <- c()
    mean_cond_first <- getMean(norm_data, de_res, 
        inputconds, 2*i-1)
    mean_cond_second <- getMean(norm_data, de_res, 
        inputconds, 2*i)
    m <- cbind(rownames(de_res), norm_data[rownames(de_res), cols],
        log10(unlist(mean_cond_first) + 0.1),
        log10(unlist(mean_cond_second) + 0.1),
        de_res[rownames(de_res),
        c("padj", "log2FoldChange", "pvalue")], 
        2 ^ de_res[rownames(de_res),
        "log2FoldChange"],
        -1 * log10(de_res[rownames(de_res), "padj"]))
    colnames(m) <- c("ID", cols, "x", "y",
        "padj", "log2FoldChange", "pvalue",
        "foldChange", "log10padj")
    m <- as.data.frame(m)
    m$padj[is.na(m[paste0("padj")])] <- 1
    m
}

#' applyFilters
#'
#' Applies filters based on user selected parameters to be
#' displayed within the DEBrowser.
#'
#' @param filt_data, loaded dataset
#' @param cols, selected samples
#' @param input, input parameters
#' @return data
#' @export
#'
#' @examples
#'     x <- applyFilters()
#'
applyFilters <- function(filt_data = NULL, cols = NULL, 
    input = NULL){
    if (is.null(input$padjtxt) || is.null(input$foldChangetxt) 
        || is.null(filt_data)) return(NULL)
    padj_cutoff <- as.numeric(input$padjtxt)
    foldChange_cutoff <- as.numeric(input$foldChangetxt)
    m <- filt_data
    # Add column which says whether a gene significant or not
    m$Legend <- character(nrow(m))
    m$Size <- character(nrow(m))
    m[, "Size"] <- "40"
    m$Legend <- "NS"
    if (input$dataset == "up" || input$dataset == "up+down") 
        m$Legend[m$foldChange >= foldChange_cutoff &
               m$padj <= padj_cutoff] <- "Up"
    if (input$dataset == "down" || input$dataset == "up+down")
        m$Legend[m$foldChange <= (1 / foldChange_cutoff) &
               m$padj <= padj_cutoff] <- "Down"
    if (input$dataset == "most-varied" && !is.null(cols)) {
        most_varied <- getMostVariedList(m, cols, input$topn, input$mincount)
        m[rownames(most_varied), c("Legend")] <- "MV"
    }
    
    if (!is.null(input$genesetarea) && input$genesetarea != "") {
        genelist <- getGeneSetData(m, c(input$genesetarea))
        m[rownames(genelist), "Legend"] <- "GS"
        m[rownames(genelist), "Size"] <- "100"
        m <- m[rev(order(m$Legend)),]
    }
    m
}

#' getSelectedDatasetInput
#'
#' Gathers the user selected dataset output to be displayed.
#'
#' @param rdata, filtered dataset
#' @param getSelected, selected data
#' @param getMostVaried, most varied data
#' @param getGeneSet, given gene set
#' @param getSelectedDatasetInput, merged comparison data
#' @param input, input parameters
#' @return data
#' @export
#'
#' @examples
#'     x <- getSelectedDatasetInput()
#'
getSelectedDatasetInput<-function(rdata = NULL, getSelected = NULL, 
    getMostVaried = NULL, getMergedComparison = NULL, 
    input = NULL) {
    if (is.null(rdata)) return (NULL)
    m <- rdata
    if (input$dataset == "up") {
        m <- rdata[rdata[, "Legend"] == "Up", ]
    } else if (input$dataset == "down") {
        m <- rdata[rdata[, "Legend"] == "Down", ]
    } else if (input$dataset == "up+down") {
        m <- rdata[which(rdata[, "Legend"] == "Down" |
            rdata[, "Legend"] == "Up"), ]
    } else if (input$dataset == "alldetected") {
        m <- rdata
    } else if (input$dataset == "selected") {
        m <- getSelected
    } else if (input$dataset == "most-varied") {
        m <- getMostVaried
    } else if (input$dataset == "comparisons") {
        m <- getMergedComparison
    }
    m
}

#' prepDataForQC
#'
#' Prepares selected data for QC plots.
#'
#' @param dataset, loaded dataset
#' @return data
#' @export
#'
#' @examples
#'     x <- prepDataForQC()
#'
prepDataForQC<-function(dataset = NULL){
    if (is.null(dataset)) return (NULL)
    columns <-colnames(dataset)
    dataset <- data.frame(dataset[,columns])
    dataset[, columns] <- apply(dataset[, columns], 2,
        function(x) as.integer(x))
    dataset1 <- rowSums(dataset[,1:ncol(dataset)])
    filtd <- data.frame(subset(dataset, rowSums(dataset[,1:ncol(dataset)]) > 10))
    norm_data <- getNormalizedMatrix(filtd)
    return(norm_data)
}

#' getMostVariedList
#'
#' Calculates the most varied genes to be used for specific plots
#' within the DEBrowser.
#'
#' @param datavar, loaded dataset
#' @param cols, selected columns
#' @param topn, most varied records
#' @param mincount, total min read count for selected samples
#' @return data
#' @export
#'
#' @examples
#'     x <- getMostVariedList()
#'
getMostVariedList <- function(datavar = NULL, cols = NULL,
    topn = 500, mincount = 10){
    if (is.null(datavar)) return (NULL)
    topn <- as.integer(as.numeric(topn))
    mincount <- as.integer(as.numeric(mincount))
    norm_data_var <- getNormalizedMatrix(
        datavar[rowSums(datavar[,cols])>mincount,cols])  
    cv<-cbind(apply(norm_data_var, 1, function(x) 
        (sd(x,na.rm=TRUE)/mean(x,na.rm=TRUE))), 1)
    colnames(cv)<-c("coeff", "a")
    cvsort<-cv[order(cv[,1],decreasing=TRUE),]
    topindex<-nrow(cvsort)
    if (topindex > topn) topindex <- topn
    cvsort_top <- head(cvsort, topindex)
    selected_var <- data.frame(datavar[rownames(cvsort_top),])
}

#' getSearchData
#'
#' search the geneset in the tables and return it
#'
#' @param dat, table data
#' @param input, input params
#' @return data
#' @export
#'
#' @examples
#'     x <- getSearchData()
#'
getSearchData <- function(dat = NULL, input = NULL)
{
  if (is.null(dat)) return(NULL)
  if (input$genesetarea != ""){
    dat <- getGeneSetData(dat, c(input$genesetarea))
  }
  dat
}

#' getGeneSetData
#'
#' Gathers the specified gene set list to be used within the
#' DEBrowser.
#'
#' @param data, loaded dataset
#' @param geneset, given gene set
#' @return data
#' @export
#'
#' @examples
#'     x <- getGeneSetData()
#'
getGeneSetData <- function(data = NULL, geneset = NULL) {
    if (is.null(data)) return (NULL)
    geneset1 <- unique(unlist(strsplit(geneset, split="[:;, \t\n\t]")))
    geneset2 <- geneset1[geneset1 != ""]
    dat1 <- data.frame(data)
    if(!("ID" %in% names(dat1)))
        dat1 <- addID(dat1)

    geneset4 <- unique(as.vector(unlist(lapply(toupper(geneset2), 
        function(x){ dat1[(grepl(x, toupper(dat1[,"ID"]))), "ID"] }))))
    retset <- data.frame(dat1[geneset4, ])
    retset
}

#' addID
#'
#' Adds an id to the data frame being used.
#'
#' @param data, loaded dataset
#' @return data
#' @export
#'
#' @examples
#'     x <- addID()
#'
addID <- function(data = NULL) {
    if (is.null(data)) return (NULL)
    dat1 <- data.frame(data)
    dat1 <- cbind(rownames(data), data)
    colnames(dat1) <- c("ID", colnames(data))
    dat1
}

#' getMergedComparison
#'
#' Gathers the merged comparison data to be used within the
#' DEBrowser.
#'
#' @param dc, data container
#' @param nc, the number of comparisons
#' @param input, input params
#' @return data
#' @export
#'
#' @examples
#'     x <- getMergedComparison()
#'
getMergedComparison <- function(dc = NULL, nc = NULL, input = NULL){
    merged <- c()
    if (is.null(dc)) return (NULL)
    
    padj_cutoff <- as.numeric(input$padjtxt)
    foldChange_cutoff <- as.numeric(input$foldChangetxt)
    for ( ni in seq(1:nc)) {
        tmp <- dc[[ni]]$init_data[,c("foldChange", "padj")]
        tt <- paste0("C", (2*ni-1),".vs.C",(2*ni))
        colnames(tmp) <- c(paste0("foldChange.", tt),  
            paste0("padj", tt))
        if (ni==1){
            merged <- tmp
        }
        else{
            merged <- merge(merged, tmp, all = TRUE, by=0)
            rownames(merged) <- merged$Row.names
            merged$Row.names <- NULL
        }
        if (is.null(merged$Legend)){
          merged$Legend <- character(nrow(merged))
          merged$Legend <- "NS"
        }
        merged[which(merged[,c(paste0("foldChange.", tt))] >= foldChange_cutoff &
            merged[,c(paste0("padj", tt))] <= padj_cutoff), "Legend"] <- "Sig"
        merged[which(merged[,c(paste0("foldChange.", tt))] <= 1/foldChange_cutoff &
            merged[,c(paste0("padj", tt))] <= padj_cutoff), "Legend"] <- "Sig"
    }
    merged <- merged[merged$Legend == "Sig", ]
    merged[,c("Legend")]<- NULL
    merged
}

#' removeCols
#'
#' remove unnecessary columns
#'
#' @param cols, columns that are going to be removed from data frame
#' @param dat, data
#' @return data
#' @export
#'
#' @examples
#'     x <- removeCols()
#'
removeCols <- function( cols = NULL, dat = NULL) {
    for (colnum in seq(1:length(cols))){
         if (cols[colnum] %in% colnames(dat) )
              dat[, cols[colnum]]<- NULL
    }
    dat
}