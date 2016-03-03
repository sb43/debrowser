#' startDEBrowser to be able to run interactively
#' 
#' @note \code{startDEBrowser}
#' @return NULL
#'
#' @examples
#'    startDEBrowser()
#'    
#' @export
#'

startDEBrowser<-function()
{
  if (interactive()) {
    #the upload file size limit is 30MB
    options( shiny.maxRequestSize = 30 * 1024 ^ 2)
    addResourcePath(prefix="demo", directoryPath=system.file("extdata/demo", package="debrowser"))
    
    environment(deServer) <- environment() 
    
    app <- list( ui = deUI,
      server = deServer)
    
    runApp(app)
  }
}