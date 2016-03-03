#' shinyUI to be able to run interactively
#' 
#' @note \code{deUI}
#' @return the panel for main plots;
#'
#' @examples
#'    x<-deUI()
#'    
#' @export
#'

deUI<- function() {
    shinyUI(fluidPage(
    titlePanel("DE Browser"),
    textOutput("text"),
    sidebarLayout(
        sidebarPanel( tags$head(tags$style(type="text/css", "
            #loadmessage {
                position: fixed;
                top: 0px;
                left: 200px;
                width: 70%;
                height: 100;
                padding: 5px 0px 5px 0px;
                text-align: center;
                font-weight: bold;
                font-size: 100%;
                color: #000000;
                opacity: 0.8;
                background-color: #FFFFFFF;
                z-index: 100;
            }
            ")),
            conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                tags$div("Please wait! Loading...",id="loadmessage",
                    tags$img(
                        src = "http://www.ajaxload.info/images/exemples/34.gif"
            ))),
            width=2,
            conditionalPanel(condition="!output.fileUploaded && !input.demo",
                actionLink("demo", "Load Demo!"),
                fileInput("file1", "Choose TSV File",
                    accept=c("text/tsv",
                        "text/comma-separated-values,text/plain",
                        ".tsv"))),
            conditionalPanel(condition="(output.fileUploaded || input.demo) && !input.goButton",
                uiOutput("condition1Selector"),
                uiOutput("condition2Selector"),
                selectInput("fittype",
                    label ="Fit Type:",
                    choices = list("parametric"="parametric",
                    "local"="local",
                    "mean"="mean"),
                    selected = 1,
                    multiple = FALSE
                ),
                actionButton("goButton", "Run DESeq!")),
            conditionalPanel(condition="input.goButton",
                h4("Filter"),
                sliderInput("padj", "padj value cut off",
                    0, 1, 0.01, step = 0.001),
            sliderInput("foldChange", "Fold Change cut off",
                    1, 10, 2, step = 0.1),
                    uiOutput("downloadSection"),
                    uiOutput("leftMenu"))
    ),
        mainPanel(
            tags$head(
                tags$style(type="text/css",
                    "#methodtabs.nav-tabs {font-size: 14px} ")),
            tabsetPanel(id="methodtabs",type = "tabs",
                tabPanel(title="Main Plots", value="panel1",
                  conditionalPanel(condition="input.demo || output.fileUploaded",
                    uiOutput("mainpanel")),
                  conditionalPanel(condition="!input.demo && !output.fileUploaded",
                    uiOutput("startup")),
                  conditionalPanel(condition="(input.demo || output.fileUploaded) && !input.goButton",
                    uiOutput("afterload"))),
                tabPanel(title="Additional Plots",value="panel2",
                    uiOutput("addpanel")),
                tabPanel(title="GO Term",  value="panel3", 
                    uiOutput("gopanel")),
                tabPanel(title="All Detected", value="panel4",
                    DT::dataTableOutput("table")),
                tabPanel(title="Up", value="panel5",
                    DT::dataTableOutput("up")),
                tabPanel(title="Down", value="panel6",
                    DT::dataTableOutput("down")),
                tabPanel(title="Selected", value="panel7",
                    DT::dataTableOutput("selected"))
            )
        )
    )
))
}
