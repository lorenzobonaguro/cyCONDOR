#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  welcome_msg <- paste("Welcome to", packageDescription("condor")$Package, "version", packageDescription("condor")$Version, "!")
  ascii_art <- "

                                                                   dddddddd
                                                                   d::::::d
                                                                   d::::::d
                                                                   d::::::d
                                                                   d:::::d
    cccccccccccccccc   ooooooooooo   nnnn  nnnnnnnn        ddddddddd:::::d    ooooooooooo   rrrrr   rrrrrrrrr
  cc:::::::::::::::c oo:::::::::::oo n:::nn::::::::nn    dd::::::::::::::d  oo:::::::::::oo r::::rrr:::::::::r
 c:::::::::::::::::co:::::::::::::::on::::::::::::::nn  d::::::::::::::::d o:::::::::::::::or:::::::::::::::::r
c:::::::cccccc:::::co:::::ooooo:::::onn:::::::::::::::nd:::::::ddddd:::::d o:::::ooooo:::::orr::::::rrrrr::::::r
c::::::c     ccccccco::::o     o::::o  n:::::nnnn:::::nd::::::d    d:::::d o::::o     o::::o r:::::r     r:::::r
c:::::c             o::::o     o::::o  n::::n    n::::nd:::::d     d:::::d o::::o     o::::o r:::::r     rrrrrrr
c:::::c             o::::o     o::::o  n::::n    n::::nd:::::d     d:::::d o::::o     o::::o r:::::r
c::::::c     ccccccco::::o     o::::o  n::::n    n::::nd:::::d     d:::::d o::::o     o::::o r:::::r
c:::::::cccccc:::::co:::::ooooo:::::o  n::::n    n::::nd::::::ddddd::::::ddo:::::ooooo:::::o r:::::r
 c:::::::::::::::::co:::::::::::::::o  n::::n    n::::n d:::::::::::::::::do:::::::::::::::o r:::::r
  cc:::::::::::::::c oo:::::::::::oo   n::::n    n::::n  d:::::::::ddd::::d oo:::::::::::oo  r:::::r
    cccccccccccccccc   ooooooooooo     nnnnnn    nnnnnn   ddddddddd   ddddd   ooooooooooo    rrrrrrr

  (o o)   (o o)
 (  V  ) (  V  )
/--m-m- /--m-m-

  "

  # Use packageStartupMessage to display the messages
  packageStartupMessage(welcome_msg)
  packageStartupMessage(ascii_art)

}

#' Print Welcome Message and ASCII Art upon Loading `ggvolc`
#'
#' This function is executed when the `ggvolc` package is attached to the R session.
#' It prints a welcome message and an ASCII representation related to the package.
#'
#' @name ggvolc-onAttach
#' @keywords internal
#' @seealso \code{\link[base]{library}}
#' @examples
#' # The function is automatically called when you use:
#' # library(ggvolc)
#` This function was inspired by loukesio/ggvolc

# This code is to inform R that the listed names are intentionally global variables
# to prevent 'no visible binding for global variable' warnings during R CMD check
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("log2FoldChange", "pvalue", "threshold","size_aes","genes" ,"x.start","x.end","y.start","y.end")) # Add other variables as needed
}
