#' Print Welcome Message and ASCII Art upon Loading `ggvolc`
#'
#' This function is executed when the `condor` package is attached to the R session.
#' It prints a welcome message and an ASCII representation related to the package.
#'
#' @name condor-onAttach
#' @keywords internal
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

# This code is to inform R that the listed names are intentionally global variables
# to prevent 'no visible binding for global variable' warnings during R CMD check
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("prcomp", "loadings", "marker", "Freq", "groups", "poi",
                           "quantile", "IQR", "UMAP1", "UMAP2", "DM_1", "DM_2",
                           "...level...", "poi", "tSNE1", "tSNE2", "model.matrix",
                           "read.delim", "separator_fc_csv", "read.delim", "nodelist",
                           "gs", "filename_col", "umap_name", "prcomp", "input_scale",
                           "colorRampPalette", "model.matrix", "trainControl",
                           "train", "Accuracy", "AccuracySD", "varImp", "..level..",
                           "DC_1", "DC_2")) # Add other variables as needed
}
