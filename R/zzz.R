#' Print Welcome Message and ASCII Art upon Loading this function is inspired by the developers of the `ggvolc` package
#'
#' This function is executed when the `cyCONDOR` package is attached to the R session.
#' It prints a welcome message and an ASCII representation related to the package.
#'
#' @name condor-onAttach
#' @keywords internal
#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  welcome_msg <- paste0("Welcome to ", packageDescription("cyCONDOR")$Package, " version ", packageDescription("cyCONDOR")$Version, "!")
  ascii_art <- "

:'######::'##:::'##::'######:::'#######::'##::: ##:'########:::'#######::'########::
'##... ##:. ##:'##::'##... ##:'##.... ##: ###:: ##: ##.... ##:'##.... ##: ##.... ##:
 ##:::..:::. ####::: ##:::..:: ##:::: ##: ####: ##: ##:::: ##: ##:::: ##: ##:::: ##:
 ##:::::::::. ##:::: ##::::::: ##:::: ##: ## ## ##: ##:::: ##: ##:::: ##: ########::
 ##:::::::::: ##:::: ##::::::: ##:::: ##: ##. ####: ##:::: ##: ##:::: ##: ##.. ##:::
 ##::: ##:::: ##:::: ##::: ##: ##:::: ##: ##:. ###: ##:::: ##: ##:::: ##: ##::. ##::
. ######::::: ##::::. ######::. #######:: ##::. ##: ########::. #######:: ##:::. ##:
:......::::::..::::::......::::.......:::..::::..::........::::.......:::..:::::..::

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
                           "DC_1", "DC_2", "PC1", "PC2", "X", "Y", "color", "comparison_result",
                           "variable", "value", "condor", "expfcs_filename", ".", "level",
                           "data_path", "max_cell", "seed", "p_val", "p_adj", "p.adj.signif",
                           "group", "value1", "value2", "y.base", "y.step", "cluster_id")) # Add other variables as needed
}
