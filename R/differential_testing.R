#' frequency_anova_test
#'
#' @title Anova test to compare cell population frequencies of independent samples
#' @description \code{frequency_anova_test()} performs an independent measures Anova to compare cell population frequencies of three of more groups. Optionally, post-hoc testing is performed using Emmeans test.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable in cell_anno that should be used to group samples in sample_var. group_var must have three or more levels.
#' @param sample_var string indicating variable in cell_anno that defines sample IDs to be used.
#' @param anova_p.adjust.method p-value adjustment method to use for multiple test correction of Anova tests, e.g "bonferroni"(default) or "BH" (Benjamini-Hochberg). All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.
#' @param post_hoc_test string specifying which post hoc test to perform. Select \code{"tukey"} (balanced data is required to perform a Tukey HSD test) or \code{"emmeans"}. By default no post-hoc test is performed.
#' @param post_hoc_p.adjust.method p-value adjustment method to use for post-hoc testing using \code{"emmeans"}, e.g "bonferroni" (default). All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.
#' @param anova_sig_threshold significance threshold of the Anova test. For all Anova tests with an adjusted p-value equal or smaller than the threshold, post-hoc tests are performed (default 0.05)
#' @param numeric logical, if TRUE numeric levels in cluster_var are ordered in ascending order and "Cluster_" is pasted before number, if FALSE alphabetical ordering is applied.
#' @details \code{frequency_anova_test()} is a wrapper function around \code{\link[rstatix]{anova_test}}, \code{\link[rstatix]{tukey_hsd}} and \code{\link[rstatix]{emmeans_test}} implemented in the package \code{rstatix}.
#' The function first calculates cell population frequencies for each sample in sample_var. Then a independent measures, one-way Anova test is performed for each cell population followed by p-value adjustment. If \code{post_hoc = T}, post-hoc testing with pairwise emmeans tests and p-value correction is performed for each significant Anova test.
#' @returns \code{frequency_anova_test()} returns a list of two data frames, "anova_test" and "emmeans_test". "anova_test" comprises results produced by \code{\link[rstatix]{anova_test}} and "emmeans_test" contains results obtained by \code{\link[rstatix]{emmeans_test}}. Both data frames have one additional columns, "cluster", containing the information, which cell population was tested.
#' @import rstatix
#' @import dplyr
#' @import reshape2
#'
#' @export
frequency_anova_test<-function(fcd,
                               cluster_slot,
                               cluster_var,
                               sample_var,
                               group_var,
                               anova_p.adjust.method = "bonferroni",
                               numeric = F,
                               post_hoc_test = NULL,
                               post_hoc_p.adjust.method = "bonferroni",
                               anova_sig_threshold = 0.05)
{

  #### check slots, cell IDs und variables
  checkInput(fcd = fcd,
             check_cluster_slot = T,
             check_cell_anno = T,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var,
             sample_var = sample_var
  )

  if (!(is.null(post_hoc_test) || post_hoc_test %in% c("emmeans", "tukey"))){
    stop('argument "post_hoc_test needs to be a valid method or NULL.')
  }
  if(!is.character(anova_p.adjust.method)){
    stop('argument "anova_p.adjust.method" needs to be a character string.')
  }
  if(!is.character(post_hoc_p.adjust.method)){
    stop('argument "post_hoc_p.adjust.method" needs to be a character string.')
  }

  #### check if group_var has more than two levels
  if(length(unique(fcd$anno$cell_anno[[group_var]])) <= 2){
    stop('group_var has two or less unique levels.')
  }

  #### prepare data
  data <- data.frame(cellID = rownames(fcd$clustering[[cluster_slot]]),
                     cluster = fcd$clustering[[cluster_slot]][[cluster_var]],
                     group_var = fcd$anno$cell_anno[[group_var]],
                     sample_var = fcd$anno$cell_anno[[sample_var]])

  if(is.numeric(data$sample_var)){
    data$sample_var <- as.character(data$sample_var)}
  if(is.numeric(data$group_var)){
    data$group_var <- as.character(data$group_var)}

  #### prepare percentage table
  tmp <- cyCONDOR::confusionMatrix(paste0(data$sample_var), paste0(data$cluster))
  tmp <- as.matrix(tmp)
  if (numeric == TRUE) {
    tmp <- tmp[order(rownames(tmp)), order(as.numeric(colnames(tmp)))]
    colnames(tmp) <- paste0("Cluster_", colnames(tmp))
  } else if (numeric == FALSE) {
    tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
  } else {
    stop('argument "numeric" needs to be set to TRUE or FALSE.')
  }

  tmp <- tmp/rowSums(tmp) * 100
  tmp <- as.data.frame(tmp)
  tmp$sample_var <- rownames(tmp)

  #### perform independent measures Anova
  tmp <- dplyr::left_join(tmp, unique(data[,c("group_var","sample_var")]), by = "sample_var")
  tmp <- reshape2::melt(tmp, id.vars=c("group_var","sample_var"))

  ##get levels from original condor object or change to factor.
  if(is.factor(data$group_var)){
    tmp$group_var <- factor(tmp$group_var, levels = levels(data$group_var))
  }else{
    tmp$group_var <- factor(tmp$group_var)
  }

  ## perform statistical test
  results <- tmp %>% dplyr::group_by(variable) %>% rstatix::anova_test(data =., value ~ group_var, detailed=F)
  results <- rstatix::adjust_pvalue(data = results, p.col = "p", method = anova_p.adjust.method)
  results$p.adj_method <- anova_p.adjust.method

  names(results)[names(results) == "variable"] <- "cluster"

  results.list<-list()
  results.list$anova_test<-results


  if(!is.null(post_hoc_test)){

    ##get variables with significant Anova test
    if(!is.numeric(anova_sig_threshold)){
      stop('argument "anova_sig_threshold" needs to be a numeric.')
    }
    cluster_keep <- unique(as.character(results[which(results$p.adj <= anova_sig_threshold),]$cluster))
    if(length(cluster_keep) == 0){
      warning('None of the ANOVAs has an adjusted p-value below the provided anova_sig_threshold '
              ,anova_sig_threshold,'.')
    }else{
      tmp_filtered <- tmp[tmp$variable %in% cluster_keep,]

      #perform tukey hsd test
      if(post_hoc_test == "tukey"){
        ##perform tukey test
        results_pht <- tmp_filtered %>%    group_by(variable) %>%  tukey_hsd(value ~ group_var,detailed = F)
        results_pht$info<-paste("Tukey HSD test for Anova tests with p.adj <=", anova_sig_threshold,".", collapse = "")

      }else if(post_hoc_test == "emmeans"){

       ##perform emmeans test
        # Check if tmp_filtered has only one unique variable value
        if(length(unique(tmp_filtered$variable)) == 1){
          results_pht  <- tmp_filtered %>%
            rstatix::emmeans_test(value ~ group_var, detailed = FALSE)
          results_pht$variable <- unique(tmp_filtered$variable)
        } else {
          # Proceed with multiple variables
          results_pht  <- bind_rows(
            lapply(unique(tmp_filtered$variable), function(i) {
              emmeans_tmp <- tmp_filtered %>%
                filter(variable == i) %>%
                rstatix::emmeans_test(value ~ group_var, detailed = FALSE)
            })
          )
        }
        ##p-value adjustment for comparisons of same variable (cluster), same adj.p as when provided during dunn_test.
        results_pht <- results_pht %>% dplyr::group_by(variable) %>%
          rstatix::adjust_pvalue(data = ., p.col = "p", method = post_hoc_p.adjust.method) %>%
          rstatix::add_significance("p.adj")
        results_pht$p.adj_method <- post_hoc_p.adjust.method
        results_pht$info<-paste("Emmeans test for Anova tests with p.adj <=",anova_sig_threshold,".", collapse = "")
      }

      names(results_pht)[names(results_pht) == "variable"] <- "cluster"
      results.list$post_hoc_test <- results_pht
    }
  }

  return(results.list)
}


#' frequency_friedman_test
#'
#' @title Friedman Rank Sum test to compare cell population frequencies of paired samples
#' @description \code{frequency_friedman_test()} performs a Friedman Rank Sum test to compare cell population frequencies of three of more groups. Optionally, post-hoc testing is performed using Wilcoxon Rank Sum test.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable in cell_anno that should be used to group samples in sample_var. group_var must have three or more levels.
#' @param sample_var string indicating variable in cell_anno that defines sample IDs to be used.
#' @param friedman_p.adjust.method p-value adjustment method to use for multiple comparisons of Friedman Rank Sum test, e.g "bonferroni" (default) or "BH" (Benjamini-Hochberg). All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.
#' @param post_hoc_test logical, whether to perform post-hoc testing (TRUE, default) or not (FALSE).
#' @param post_hoc_p.adjust.method p-value adjustment method to use for post-hoc testing, e.g "bonferroni" (default). All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.
#' @param friedman_sig_threshold significance threshold Friedman Rank Sum test. For all Friedman Rank Sum comparisons with an adjusted p-value equal or smaller than the threshold, post-hoc tests are performed (default 0.05)
#' @param numeric logical, if TRUE numeric levels in cluster_var are ordered in ascending order and "Cluster_" is pasted before number, if FALSE alphabetical ordering is applied.
#' @param pair_var string indicating variable in cell_anno that should be used to pair the samples.
#' @details \code{frequency_friedman_test()} is a wrapper function around \code{\link[rstatix]{friedman_test}},  \code{\link[rstatix]{friedman_effsize}} and  \code{\link[rstatix]{wilcox_test}} implemented in the package \code{rstatix}. The function first calculates cell population frequencies for each sample in sample_var. Then a Friedman Rank Sum test is performed for each cell population followed by p-value adjustment. If \code{post_hoc = T}, post-hoc testing with pairwise Wilcoxon Rank Sum Tests and p-value correction is performed for each significant Friedman Rank Sum test comparison.
#' @returns \code{frequency_friedman_test} returns a list of two data frames, "friedman_test" and "wilcox_test". "friedman_test" comprises results produced by \code{\link[rstatix]{friedman_test}} and \code{\link[rstatix]{friedman_effsize}} and "wilcox_test" contains results obtained by \code{\link[rstatix]{wilcox_test}}. Both data frames have one additional columns, "cluster", containing the information, which cell population was tested.
#' @import rstatix
#' @import dplyr
#' @import reshape2
#'
#' @export
frequency_friedman_test<-function(fcd,
                                  cluster_slot,
                                  cluster_var,
                                  sample_var,
                                  group_var,
                                  pair_var,
                                  friedman_p.adjust.method = "bonferroni",
                                  numeric = F,
                                  post_hoc_test = F,
                                  post_hoc_p.adjust.method = "bonferroni",
                                  friedman_sig_threshold = 0.05)
{

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_cluster_slot = T,
             check_cell_anno = T,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var,
             sample_var = sample_var,
             pair_var = pair_var
  )
  if(!(isTRUE(post_hoc_test) || isFALSE(post_hoc_test))){
    stop('argument "post_hoc_test needs to be TRUE or FALSE.')
  }
  if(!is.character(friedman_p.adjust.method)){
    stop('argument "friedman_p.adjust.method" needs to be a character string.')
  }
  if(!is.character(post_hoc_p.adjust.method)){
    stop('argument "post_hoc_p.adjust.method" needs to be a character string.')
  }


  #### check if group_var has more than two levels
  if(length(unique(fcd$anno$cell_anno[[group_var]])) <= 2){
    stop('group_var has two or less unique levels.')
  }else{
    n_groups<-length(unique(fcd$anno$cell_anno[[group_var]]))
  }

  #### prepare data
  data <- data.frame(cellID = rownames(fcd$clustering[[cluster_slot]]),
                     cluster = fcd$clustering[[cluster_slot]][[cluster_var]],
                     group_var = fcd$anno$cell_anno[[group_var]],
                     sample_var = fcd$anno$cell_anno[[sample_var]],
                     pair_var = fcd$anno$cell_anno[[pair_var]])


  if(is.numeric(data$sample_var)){
    data$sample_var <- as.character(data$sample_var)}

  #### make sure that all data pairs are fine.
  pairing<-unique(data[,c("sample_var","pair_var","group_var")])

  if(sum(is.na(pairing))>0){
    stop('Paired data connot be identified.
         NAs were found in sample_var, pair_var or group_var. Please remove samples with missing information.')
  }
  if(length(unique(pairing$sample_var)) < nrow(pairing)){
    stop('Paired data cannot be identified.
         Make sure that each sample ID in sample_var is uniquely assigned to one donor+group combination.')
  }else{
    tmp <- data.frame(base::table(pairing$pair_var))
    if(!sum(tmp$Freq==n_groups)==nrow(tmp)){
      stop('Paired data cannot be identified.
           Make sure that each donor has exactly one sample for each group_var.(1)')
    }else{
      tmp <- data.frame(base::table(pairing$pair_var, pairing$group_var))
      if(!sum(tmp$Freq==1)==nrow(tmp)){
        stop('Paired data cannot be identified.
             Make sure that each donor has exactly one sample for each group_var.(2)')
      }
    }
  }


  #### prepare percentage table
  tmp <- cyCONDOR::confusionMatrix(paste0(data$sample_var), paste0(data$cluster))
  tmp <- as.matrix(tmp)
  if (numeric == TRUE) {
    tmp <- tmp[order(rownames(tmp)), order(as.numeric(colnames(tmp)))]
    colnames(tmp) <- paste0("Cluster_", colnames(tmp))
  } else if (numeric == FALSE) {
    tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
  } else {
    #tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
    stop('argument "numeric" needs to be set to TRUE or FALSE.')
  }

  tmp <- tmp/rowSums(tmp) * 100
  tmp <- as.data.frame(tmp)
  tmp$sample_var <- rownames(tmp)

  #### perform friedman test
  tmp <- dplyr::left_join(tmp, unique(data[,c("group_var","sample_var","pair_var")]), by = "sample_var")
  tmp <- reshape2::melt(tmp, id.vars=c("group_var","sample_var","pair_var"))


  ## ensure correct order of samples and check pairs
  tmp <- tmp %>% dplyr::group_by(variable) %>% dplyr::arrange(pair_var, .by_group = TRUE)
  # if(!identical(tmp[tmp$group_var==unique(tmp$group_var)[1],]$pair_var,
  #               tmp[tmp$group_var==unique(tmp$group_var)[2],]$pair_var)){
  #   stop('Something unexpected is wrong in pairing of sample_var and pair_var.') #should not be needed, but better safe than sorry.
  # }

  ##get levels from original condor object or change to factor.
  if(is.factor(data$group_var)){
    tmp$group_var <- factor(tmp$group_var, levels = levels(data$group_var))
  }else{
    tmp$group_var <- factor(tmp$group_var)
  }

  ## perform statistical test (paired)
  results <- tmp %>% dplyr::group_by(variable) %>% rstatix::friedman_test(data =., value ~ group_var | pair_var)
  results <- rstatix::adjust_pvalue(data = results, p.col = "p", method = friedman_p.adjust.method)
  results$p.adj_method<-friedman_p.adjust.method

  ## calculate effect size
  effect_size <- tmp %>% dplyr::group_by(variable) %>% rstatix::friedman_effsize(data =., value ~ group_var | pair_var)
  effect_size <- effect_size[,c("variable","effsize","n","method","magnitude")]
  colnames(effect_size)[3:5] <- paste0("effsize_",colnames(effect_size)[3:5])
  results <- dplyr::left_join(results, effect_size, by = "variable")

  names(results)[names(results) == "variable"] <- "cluster"


  results.list<-list()
  results.list$friedman_test<-results


  if(post_hoc_test == T){
    #### post hoc test (wilcox test)

    ##get variables with significant Friedman test
    if(!is.numeric(friedman_sig_threshold)){
      stop('argument "friedman_sig_threshold" needs to be a numeric.')
    }
    cluster_keep <- unique(as.character(results[which(results$p.adj <= friedman_sig_threshold),]$cluster))
    if(length(cluster_keep) == 0){
      warning('None of the Friedman tests has an adjusted p-value below the provided friedman_sig_threshold '
              ,friedman_sig_threshold,'.')
    }else{
      tmp_filtered <- tmp[tmp$variable %in% cluster_keep,]

      ##perform post hoc test (wilcox)
      results_wilcox <- tmp_filtered %>% dplyr::group_by(variable) %>% rstatix::wilcox_test(data =., value ~ group_var, detailed = T, paired = T)
      ##p-value adjustment for comparisons of same variable (cluster), same adj.p as when provided during dunn_test.
      results_wilcox <- results_wilcox %>% dplyr::group_by(variable) %>%
        rstatix::adjust_pvalue(data = ., p.col = "p", method = post_hoc_p.adjust.method) %>%
        rstatix::add_significance("p.adj")
      results_wilcox$p.adj_method <- post_hoc_p.adjust.method
      results_wilcox$info<-paste("Wilcox test for Friedman tests with p.adj <=",friedman_sig_threshold,".", collapse = "")
      # ##p-value adjustment for all comparisons
      # results_wilcox <- rstatix::adjust_pvalue(data = results_wilcox, p.col = "p", method = post_hoc_p.adjust.method) %>%
      #   rstatix::add_significance("p.adj")

      names(results_wilcox)[names(results_wilcox) == "variable"] <- "cluster"
      results.list$wilcox_test <- results_wilcox
    }
  }

  return(results.list)
}


#' frequency_kruskal_test
#'
#' @title Kruskal-Wallis test to compare cell population frequencies
#' @description \code{frequency_kruskal_test()} performs a Kruskal-Wallis Rank Sum test to compare cell population frequencies of three of more groups. Optionally, post-hoc testing is performed using Dunne's test.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable in cell_anno that should be used to group samples in sample_var. group_var must have three or more levels.
#' @param sample_var string indicating variable in cell_anno that defines sample IDs to be used.
#' @param kruskal_p.adjust.method p-value adjustment method to use for multiple comparisons of Kruskal-Wallis test, e.g "bonferroni" (default) or "BH" (Benjamini-Hochberg). All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.
#' @param post_hoc_test logical, whether to perform post-hoc testing (TRUE, default) or not (FALSE).
#' @param post_hoc_p.adjust.method p-value adjustment method to use for post-hoc testing, e.g "bonferroni" (default). All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.
#' @param kruskal_sig_threshold significance threshold for Kruskal-Wallis test. For all Kruskal-Wallis comparisons with an adjusted p-value equal or smaller than the threshold, post-hoc tests are performed (default 0.05)
#' @param numeric logical, if TRUE numeric levels in cluster_var are ordered in ascending order and "Cluster_" is pasted before number, if FALSE alphabetical ordering is applied.
#' @details \code{frequency_kruskal_test()} is a wrapper function around \code{\link[rstatix]{kruskal_test}},  \code{\link[rstatix]{kruskal_effsize}} and  \code{\link[rstatix]{dunn_test}} implemented in the package *rstatix*. The function first calculates cell population frequencies for each sample in sample_var. Then a Kruskal-Wallis rank sum test is performed for each cell population followed by p-value adjustment. If \code{post_hoc = T}, post-hoc testing with Dunne's Test and p-value correction is performed for each significant Kruskal-Wallis comparison.
#' @returns \code{frequency_kruskal_test()} returns a list of two data frames, "kruskal_test" and "dunn_test". "kruskal_test" comprises results produced by \code{\link[rstatix]{kruskal_test}} and \code{\link[rstatix]{kruskal_effsize}} and "dunn_test" contains results obtained by \code{\link[rstatix]{dunn_test}}. Both data frames have one additional columns, "cluster", containing the information, which cell population was tested.
#' @import rstatix
#' @import dplyr
#' @import reshape2
#'
#' @export
frequency_kruskal_test<-function(fcd,
                                 cluster_slot,
                                 cluster_var,
                                 sample_var,
                                 group_var,
                                 kruskal_p.adjust.method = "bonferroni",
                                 post_hoc_test = T,
                                 post_hoc_p.adjust.method = "bonferroni",
                                 kruskal_sig_threshold = 0.05,
                                 numeric = F)
{

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_cluster_slot = T,
             check_cell_anno = T,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var,
             sample_var = sample_var
  )

  if(!(isTRUE(post_hoc_test) || isFALSE(post_hoc_test))){
    stop('argument "post_hoc_test needs to be TRUE or FALSE.')
  }
  if(!is.character(kruskal_p.adjust.method)){
    stop('argument "kruskal_p.adjust.method" needs to be a character string.')
  }
  if(!is.character(post_hoc_p.adjust.method)){
    stop('argument "post_hoc_p.adjust.method" needs to be a character string.')
  }

  #### check if group_var has more than two levels
  if(length(unique(fcd$anno$cell_anno[[group_var]])) <= 2){
    stop('group_var has two or less unique levels.')
  }

  #### prepare data
  data <- data.frame(cellID = rownames(fcd$clustering[[cluster_slot]]),
                     cluster = fcd$clustering[[cluster_slot]][[cluster_var]],
                     group_var = fcd$anno$cell_anno[[group_var]],
                     sample_var = fcd$anno$cell_anno[[sample_var]])

  if(is.numeric(data$sample_var)){
    data$sample_var <- as.character(data$sample_var)}

  #### prepare percentage table
  tmp <- cyCONDOR::confusionMatrix(paste0(data$sample_var), paste0(data$cluster))
  tmp <- as.matrix(tmp)
  if (numeric == TRUE) {
    tmp <- tmp[order(rownames(tmp)), order(as.numeric(colnames(tmp)))]
    colnames(tmp) <- paste0("Cluster_", colnames(tmp))
  } else if (numeric == FALSE) {
    tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
  } else {
    #tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
    stop('argument "numeric" needs to be set to TRUE or FALSE.')
  }

  tmp <- tmp/rowSums(tmp) * 100
  tmp <- as.data.frame(tmp)
  tmp$sample_var <- rownames(tmp)


  #### perform Kruskal-Wallis Rank Sum Test
  tmp <- dplyr::left_join(tmp, unique(data[,c("group_var","sample_var")]), by = "sample_var")
  tmp <- reshape2::melt(tmp, id.vars=c("group_var","sample_var"))

  ## get levels from original condor object or change to factor
  if(is.factor(data$group_var)){
    tmp$group_var <- factor(tmp$group_var, levels = levels(data$group_var))
  }else{
    tmp$group_var <- factor(tmp$group_var)
  }

  ## perform statistical test
  results <- tmp %>% dplyr::group_by(variable) %>% rstatix::kruskal_test(data =., value ~ group_var)
  results <- rstatix::adjust_pvalue(data = results, p.col = "p", method = kruskal_p.adjust.method)
  results$p.adj_method <- kruskal_p.adjust.method

  ## calculate effect size
  effect_size <- tmp %>% dplyr::group_by(variable) %>% rstatix::kruskal_effsize(data =., value ~ group_var)
  effect_size <- effect_size[,c("variable","effsize","n","method","magnitude")]
  colnames(effect_size)[3:5] <- paste0("effsize_",colnames(effect_size)[3:5])
  results <- dplyr::left_join(results, effect_size, by = "variable")

  names(results)[names(results) == "variable"] <- "cluster"

  results.list<-list()
  results.list$kruskal_test<-results


  if(post_hoc_test == T){
    #### post hoc test (dunn test)

    ## get variables with significant Kruskal-Wallis test
    if(!is.numeric(kruskal_sig_threshold)){
      stop('argument "kruskal_sig_threshold" needs to be a numeric.')
    }
    cluster_keep <- unique(as.character(results[which(results$p.adj <= kruskal_sig_threshold),]$cluster))
    if(length(cluster_keep) == 0){
      warning('None of the Kruskal-Wallis tests has an adjusted p-value below the provided kruskal_sig_threshold '
              ,kruskal_sig_threshold,'.')
    }else{
      tmp_filtered <- tmp[tmp$variable %in% cluster_keep,]

      ##perform dunn test
      results_dunn <- tmp_filtered %>% dplyr::group_by(variable) %>% rstatix::dunn_test(data =., value ~ group_var, detailed = F)

      ##p-value adjustment for comparisons of same variable (cluster), same adj.p as when provided during dunn_test.
      results_dunn <- results_dunn %>% dplyr::group_by(variable) %>%
        rstatix::adjust_pvalue(data = ., p.col = "p", method = post_hoc_p.adjust.method) %>%
        rstatix::add_significance("p.adj")
      results_dunn$p.adj_method <- post_hoc_p.adjust.method
      results_dunn$info<-paste("Dunn test for Kruskal-Wallis tests with p.adj <=",kruskal_sig_threshold,".", collapse = "")
      # ##p-value adjustment for all comparisons
      # results_dunn <- rstatix::adjust_pvalue(data = results_dunn, p.col = "p", method = post_hoc_p.adjust.method) %>%
      #   rstatix::add_significance("p.adj")

      names(results_dunn)[names(results_dunn) == "variable"] <- "cluster"
      results.list$dunn_test <- results_dunn
    }
  }

  return(results.list)
}


#' frequency_t_test
#'
#' @title t-test to compare cell population frequencies
#' @description \code{frequency_t_test()} performs a two-sided, two sample t-test on cell type frequencies.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable in cell_anno that should be used to group samples in sample_var. group_var must have exactly two levels.
#' @param sample_var string indicating variable in cell_anno that defines sample IDs to be used.
#' @param pair_var string indicating variable in cell_anno that defines pairing of samples, e.g. donor ID, that should be used if \code{paired_test = T}.
#' @param paired_test logical, indicating if a paired (TRUE) or unpaired test (FALSE, default) should be performed.
#' @param var.equal logical, indicating whether variance in both groups should be treated as equal (TRUE) or not (FALSE, default). TRUE uses pooled variance, FALSE the Welch approximation, see documentation of \code{\link[rstatix]{t_test}} .
#' @param detailed logical if detailed output from \code{\link[rstatix]{t_test}} should be reported.
#' @param p.adjust.method p-value adjustment method to use for multiple comparison testing, e.g "bonferroni" (default) or "BH" (Benjamini-Hochberg). All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.
#' @param numeric logical, if TRUE numeric levels in cluster_var are ordered in ascending order and "Cluster_" is pasted before number, if FALSE alphabetical ordering is applied.
#' @param print_results Logical, indicating if the test results are printed to the console (TRUE) or not (FALSE).
#' @details \code{frequency_t_test()} is a wrapper function around the \code{t_test()} implemented in the package \code{rstatix}.
#' The function first calculates cell population frequencies for each sample in sample_var. Then, a two-sided, two sample t-test is performed between two groups defined in group_var.
#' The test can either be run unpaired (two independent groups) or paired. Afterwards,  p-value adjustment will be performed across all comparisons that were made.
#' @returns \code{frequency_t_test()} returns the fcd containing a data frame produced by \code{\link[rstatix]{t_test}} with two additional columns, "cluster" containing the information, which cell population was tested, and "applied_test", indicating which test was used. Results are stored in the fcd under extras$statistics.
#' @import rstatix
#' @import dplyr
#' @import reshape2
#'
#' @export
frequency_t_test<-function(fcd,
                           cluster_slot,
                           cluster_var,
                           sample_var,
                           group_var,
                           pair_var = NULL,
                           paired_test = F,
                           var.equal = F,
                           detailed = F,
                           p.adjust.method = "bonferroni",
                           numeric = F,
                           print_results = T)
{

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_cluster_slot = T,
             check_cell_anno = T,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var,
             sample_var = sample_var,
             pair_var = pair_var
  )

  if(!(isTRUE(paired_test) || isFALSE(paired_test))){
    stop('argument "paired_test" needs to be TRUE or FALSE.')
  }
  if(!is.character(p.adjust.method)){
    stop('argument "p.adjust.method" needs to be a character string.')
  }

  #### check if group_var has exactly two levels
  if(length(unique(fcd$anno$cell_anno[[group_var]])) > 2){
    stop('group_var has more than two unique levels.')
  }

  #### prepare data
  data <- data.frame(cellID = rownames(fcd$clustering[[cluster_slot]]),
                     cluster = fcd$clustering[[cluster_slot]][[cluster_var]],
                     group_var = fcd$anno$cell_anno[[group_var]],
                     sample_var = fcd$anno$cell_anno[[sample_var]])

  if(is.numeric(data$sample_var)){
    data$sample_var <- as.character(data$sample_var)}

  #### take care of pairing
  if(paired_test == F & !is.null(pair_var)){
    warning('The argument pair_var will be ignored, since argument "paired_test" is set to FALSE. Please note that the test is run without pairing!')
  }
  if(paired_test == T){
    if(is.null(pair_var)){
      stop('paired_test is set to T. Please provide a variable for pairing via argument pair_var')

    }else{
      data$pair_var <- as.character(fcd$anno$cell_anno[[pair_var]])

      #### make sure that all data pairs are fine.
      pairing<-unique(data[,c("sample_var","pair_var","group_var")])

      if(sum(is.na(pairing)) > 0){
        stop('Paired data connot be identified.
         NAs were found in sample_var, pair_var or group_var. Please remove samples with missing information.')
      }
      if(length(unique(pairing$sample_var)) < nrow(pairing)){
        stop('Paired data cannot be identified.
         Make sure that each sample ID in sample_var is uniquely assigned to one donor+group combination.')
      }

      tmp <- data.frame(base::table(pairing$pair_var))
      if(!sum(tmp$Freq==2) == nrow(tmp)){
        stop('Paired data cannot be identified.
           Make sure that each donor has exactly one sample for each group_var.(1)')
      }
      tmp <- data.frame(base::table(pairing$pair_var, pairing$group_var))
      if(!sum(tmp$Freq == 1) == nrow(tmp)){
        stop('Paired data cannot be identified.
             Make sure that each donor has exactly one sample for each group_var.(2)')
      }
    }
  }


  #### prepare percentage table
  tmp <- cyCONDOR::confusionMatrix(paste0(data$sample_var), paste0(data$cluster))
  tmp <- as.matrix(tmp)
  if (numeric == TRUE) {
    tmp <- tmp[order(rownames(tmp)), order(as.numeric(colnames(tmp)))]
    colnames(tmp) <- paste0("Cluster_", colnames(tmp))
  } else if (numeric == FALSE) {
    tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
  } else {
    #tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
    stop('argument "numeric" needs to be set to TRUE or FALSE.')
  }

  tmp <- tmp/rowSums(tmp) * 100
  tmp <- as.data.frame(tmp)
  tmp$sample_var <- rownames(tmp)


  #### perform wilcoxon test
  if(paired_test == T){

    tmp <- dplyr::left_join(tmp, unique(data[,c("group_var","sample_var","pair_var")]), by = "sample_var")
    tmp <- reshape2::melt(tmp, id.vars=c("group_var","sample_var","pair_var"))

    ## ensure correct order of samples and check pairs
    tmp <- tmp %>% dplyr::group_by(variable) %>% dplyr::arrange(pair_var, .by_group = TRUE)
    if(!identical(tmp[tmp$group_var==unique(tmp$group_var)[1],]$pair_var,
                  tmp[tmp$group_var==unique(tmp$group_var)[2],]$pair_var)){
      stop('Something unexpected is wrong in pairing of sample_var and pair_var.') #should not be needed, but better safe than sorry.
    }

    ## get levels from original condor object. Allows user to adjust group1 and group2 via levels of group_var
    if(is.factor(data$group_var)){
      tmp$group_var <- factor(tmp$group_var, levels = levels(data$group_var))
    }else{
      tmp$group_var <- factor(tmp$group_var)
    }

    ## perform statistical test (paired)
    results <- tmp %>% dplyr::group_by(variable) %>% rstatix::t_test(data =., value ~ group_var, detailed = detailed, paired = T, var.equal = var.equal)
    results <- rstatix::adjust_pvalue(data = results, p.col = "p", method = p.adjust.method) %>%
      rstatix::add_significance("p.adj")

    names(results)[names(results) == "variable"] <- "cluster"
    results$p.adj_method<-p.adjust.method
    results$applied_test<-"paired t test"

  }else{

    tmp <- dplyr::left_join(tmp, unique(data[,c("group_var","sample_var")]), by = "sample_var")
    tmp <- reshape2::melt(tmp, id.vars=c("group_var","sample_var"))

    ## get levels from original condor object. Allows user to adjust group1 and group2 via levels of group_var
    if(is.factor(data$group_var)){
      tmp$group_var <- factor(tmp$group_var, levels = levels(data$group_var))
    }else{
      tmp$group_var <- factor(tmp$group_var)
    }

    ## perform statistical test (unpaired)
    results <- tmp %>% dplyr::group_by(variable) %>% rstatix::t_test(data =., value ~ group_var, detailed = detailed, paired = F, var.equal = var.equal)
    results <- rstatix::adjust_pvalue(data = results, p.col = "p", method = p.adjust.method) %>%
      rstatix::add_significance("p.adj")

    names(results)[names(results) == "variable"] <- "cluster"
    results$p.adj_method<-p.adjust.method
    results$applied_test<-"unpaired t test"
  }

  #### store the results in the fcd
  ## check if the statistics slot exists. If not, initiate a list for storage.
  if(is.null(fcd$extras$statistics)){
    fcd$extras$statistics <- list()
  }

  ## assing the results to their unique slot
  fcd$extras$statistics[["t_test"]] <- results
  message("Statistic results were saved in the fcd under extras$statistics")

  #### print the results to the console if desired
  if(print_results == T){
    print(results)
  }

  return(fcd)
}


#' frequency_wilcox_test
#'
#' @title Wilcoxon test to compare cell population frequencies
#' @description \code{frequency_wilcox_test()} performs a two-sided, two-sample Wilcoxon Rank Sum test (independent) or Wilcoxon Signed Rank test (paired data) on cell type frequencies.
#' @param fcd flow cytometry data set, that has been subjected to clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable  in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable  in cell_anno that should be used to group samples in sample_var. group_var must have exactly two levels.
#' @param sample_var string indicating variable  in cell_anno that defines sample IDs to be used.
#' @param pair_var string indicating variable  in cell_anno that defines pairing of samples, e.g. donor ID, that should be used if \code{paired_test = T}.
#' @param paired_test logical, indicating if a paired (TRUE) or unpaired test (FALSE, default) should be performed.
#' @param p.adjust.method p-value adjustment method to use for multiple comparison testing, e.g "bonferroni" (default) or "BH" (Benjamini-Hochberg). All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.
#' @param detailed logical if detailed output from \code{\link[rstatix]{wilcox_test}} should be reported.
#' @param numeric logical, if TRUE numeric levels in cluster_var are ordered in ascending order and "Cluster_" is pasted before number, if FALSE alphabetical ordering is applied.
#' @param print_results Logical, indicating if the test results are printed to the console (TRUE) or not (FALSE).
#' @details \code{frequency_wilcox_test()} is a wrapper function around the \code{wilcox_test()} function implemented in the package \code{rstatix}.
#' The function first calculates cell population frequencies for each sample in sample_var. Then a two-sided, two sample wilcoxon test is performed between two groups defined in group_var.
#' The test can either be run unpaired (two independent groups) or paired. Afterwards,  p-value adjustment will be performed across all comparisons that were made.
#' @returns \code{frequency_wilcox_test} returns the fcd containing a data frame produced by \code{\link[rstatix]{wilcox_test}} with two additional columns, "cluster" containing the information, which cell population was tested, and "applied_test", indicating which test was used. Results are stored in the fcd under extras$statistics.
#' @import rstatix
#' @import dplyr
#' @import reshape2
#'
#' @export
frequency_wilcox_test<-function(fcd,
                                cluster_slot,
                                cluster_var,
                                sample_var,
                                group_var,
                                pair_var = NULL,
                                paired_test = F,
                                p.adjust.method = "bonferroni",
                                detailed = F,
                                numeric = F,
                                print_results = T)
{

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_cluster_slot = T,
             check_cell_anno = T,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var,
             sample_var = sample_var,
             pair_var = pair_var
  )

  if(!(isTRUE(paired_test) || isFALSE(paired_test))){
    stop('argument "paired_test" needs to be TRUE or FALSE.')
  }
  if(!is.character(p.adjust.method)){
    stop('argument "p.adjust.method" needs to be a character string')
  }

  #### check if group_var has exactly two levels
  if(length(unique(fcd$anno$cell_anno[[group_var]])) > 2){
    stop('group_var has more than two unique levels.')
  }

  #### prepare data
  data <- data.frame(cellID = rownames(fcd$clustering[[cluster_slot]]),
                     cluster = fcd$clustering[[cluster_slot]][[cluster_var]],
                     group_var = fcd$anno$cell_anno[[group_var]],
                     sample_var = fcd$anno$cell_anno[[sample_var]])

  if(is.numeric(data$sample_var)){
    data$sample_var <- as.character(data$sample_var)}

  #### take care of pairing
  if(paired_test == F & !is.null(pair_var)){
    warning('The argument pair_var will be ignored, since argument "paired_test" is set to FALSE. Please note that the test is run without pairing!')
  }
  if(paired_test == T){
    if(is.null(pair_var)){
      stop('paired_test is set to T. Please provide a variable for pairing via argument pair_var')

    }else{
      data$pair_var <- as.character(fcd$anno$cell_anno[[pair_var]])

      #### make sure that all data pairs are fine.
      pairing<-unique(data[,c("sample_var","pair_var","group_var")])

      if(sum(is.na(pairing)) > 0){
        stop('Paired data connot be identified.
         NAs were found in sample_var, pair_var or group_var. Please remove samples with missing information.')
      }
      if(length(unique(pairing$sample_var)) < nrow(pairing)){
        stop('Paired data cannot be identified.
         Make sure that each sample ID in sample_var is uniquely assigned to one donor+group combination.')
      }

      tmp <- data.frame(base::table(pairing$pair_var))
      if(!sum(tmp$Freq==2) == nrow(tmp)){
        stop('Paired data cannot be identified.
           Make sure that each donor has exactly one sample for each group_var.(1)')
      }
      tmp <- data.frame(base::table(pairing$pair_var, pairing$group_var))
      if(!sum(tmp$Freq == 1) == nrow(tmp)){
        stop('Paired data cannot be identified.
             Make sure that each donor has exactly one sample for each group_var.(2)')
      }
    }
  }


  #### prepare percentage table
  tmp <- cyCONDOR::confusionMatrix(paste0(data$sample_var), paste0(data$cluster))
  tmp <- as.matrix(tmp)
  if (numeric == TRUE) {
    tmp <- tmp[order(rownames(tmp)), order(as.numeric(colnames(tmp)))]
    colnames(tmp) <- paste0("Cluster_", colnames(tmp))
  } else if (numeric == FALSE) {
    tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
  } else {
    #tmp <- tmp[order(rownames(tmp)), order(colnames(tmp))]
    stop('argument "numeric" needs to be set to TRUE or FALSE.')
  }

  tmp <- tmp/rowSums(tmp) * 100
  tmp <- as.data.frame(tmp)
  tmp$sample_var <- rownames(tmp)

  #### perform wilcoxon test
  if(paired_test == T){

    tmp <- dplyr::left_join(tmp, unique(data[,c("group_var","sample_var","pair_var")]), by = "sample_var")
    tmp <- reshape2::melt(tmp, id.vars=c("group_var","sample_var","pair_var"))

    ## ensure correct order of samples and check pairs
    tmp <- tmp %>% dplyr::group_by(variable) %>% dplyr::arrange(pair_var, .by_group = TRUE)
    if(!identical(tmp[tmp$group_var==unique(tmp$group_var)[1],]$pair_var,
                  tmp[tmp$group_var==unique(tmp$group_var)[2],]$pair_var)){
      stop('Something unexpected is wrong in pairing of sample_var and pair_var.') #should not be needed, but better safe than sorry.
    }

    ## get levels from original condor object. Allows user to adjust group1 and group2 via levels of group_var
    if(is.factor(data$group_var)){
      tmp$group_var <- factor(tmp$group_var, levels = levels(data$group_var))
    }else{
      tmp$group_var <- factor(tmp$group_var)
    }

    ## perform statistical test (paired)
    results <- tmp %>% dplyr::group_by(variable) %>% rstatix::wilcox_test(data =., value ~ group_var, detailed = detailed, paired = T)
    results <- rstatix::adjust_pvalue(data = results, p.col = "p", method = p.adjust.method) %>%
      rstatix::add_significance("p.adj")

    names(results)[names(results) == "variable"] <- "cluster"
    results$p.adj_method<-p.adjust.method
    results$applied_test<-"paired wilcox test"

  }else{

    tmp <- dplyr::left_join(tmp, unique(data[,c("group_var","sample_var")]), by = "sample_var")
    tmp <- reshape2::melt(tmp, id.vars=c("group_var","sample_var"))

    ##get levels from original condor object. Allows user to adjust group1 and group2 via levels of group_var
    if(is.factor(data$group_var)){
      tmp$group_var <- factor(tmp$group_var, levels = levels(data$group_var))
    }else{
      tmp$group_var <- factor(tmp$group_var)
    }

    ## perform statistical test (unpaired)
    results <- tmp %>% dplyr::group_by(variable) %>% rstatix::wilcox_test(data =., value ~ group_var, detailed = detailed, paired = F)
    results <- rstatix::adjust_pvalue(data = results, p.col = "p", method = p.adjust.method) %>%
      rstatix::add_significance("p.adj")

    names(results)[names(results) == "variable"] <- "cluster"
    results$p.adj_method<-p.adjust.method
    results$applied_test<-"unpaired wilcox test"
  }

  #### store the results in the fcd
  ## check if the statistics slot exists. If not, initiate a list for storage.
  if(is.null(fcd$extras$statistics)){
    fcd$extras$statistics <- list()
  }

  ## assing the results to their unique slot
  fcd$extras$statistics[["wilcox"]] <- results
  message("Statistic results were saved in the fcd under extras$statistics")

  #### print the results to the console if desired
  if(print_results == T){
    print(results)
  }

  return(fcd)
}

#' prepInputDiffcyt
#'
#' @title Convert a condor object to se object
#' @description \code{prepInputDiffcyt()} converts a fcd object into a SummarizedExperiment object compatible with \code{diffcyt} functions \code{\link[diffcyt]{calcCounts}} and \code{\link[diffcyt]{calcMedians}}.
#' @param fcd flow cytometry dataset, that has been subjected to clustering or cell type label prediction with cyCONDOR before
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels)
#' @param sample_var Charlotte
#' @param meta_vars vector of variables in cell_anno, which contain sample level metadata, which means that each sample ID is associated with exactly one level per variable. All variables that the user wants to use in the test design need to be listed, e.g. group, donor_id. Variables with names "sample_id" and "cluster_id" are not allowed, since these names have designated purposes in diffcyt workflow.
#' @param marker_state vector of marker names that should get the marker_class "state". If no markers are provided in marker_state and marker_type all available markers and features in expr data will be set as "type".
#' @param marker_type vector of marker names available in expr data, that should get the marker_class "type". If no markers are provided in marker_state and marker_type all available markers and features in expr data will get the marker_class "type".
#' @details The function will carry over the original transformed expression. The flexible experimental design of  diffcyt's testing functions allows to include batch variables.
#' @import diffcyt
#' @import SummarizedExperiment
#' @returns
#' A SummarizedExperiment object suitable to be used as input for \code{diffcyt} functions \code{\link[diffcyt]{calcCounts}} and \code{\link[diffcyt]{calcMedians}}. The object contains the following components:
#' \itemize{
#'  \item{metadata "experiment_info"} : sample-level metadata table, containing all variables provided in meta_vars and sample_var, whereby sample_var is renamed to "sample_id"
#'  \item{metadata "n_cells"} : number of cells per sample_id
#'  \item{assay "exprs"} : contains expression data
#'  \item{rowData} : cell-level information, containing all variables provided in meta_vars, sample_var and cluster_var, whereby sample_var is renamed to "sample_id" and cluster_var tp "cluster_id" to be compatible with diffcyt workflow.
#'  \item{colData} : marker information required for diffcyt differential state analysis.
#' }
#'
#' @export
prepInputDiffcyt<-function(fcd,
                           cluster_slot,
                           cluster_var,
                           sample_var,
                           meta_vars,
                           marker_state = NULL,
                           marker_type = NULL){

  expr_slot <- "orig"

  #### check slots, cell IDs and variables
  checkInput(fcd = fcd,
             check_expr_slot = T,
             check_cluster_slot = T,
             check_cell_anno = T,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             sample_var = sample_var)

  if(is.null(sample_var)){
    stop("sample_var needs to be set to run this function.")
  }
  if(is.null(cluster_var)){
    stop("cluster_var needs to be set to run this function.")
  }


  #### prepare marker_info data frame
  if(length(intersect(marker_type,marker_state)) >= 1){
    stop("Error: Specified markers in marker_state and marker_type are overlapping.")
  }

  ## type marker
  if(!is.null(marker_type)){
    marker_type<-unique(marker_type)

    ##check if markers are present in expr slot
    marker_type_present<-marker_type[marker_type %in% colnames(fcd$expr[[expr_slot]])]
    if(length(marker_type_present) == 0){
      stop('None of the provided type marker are present in expr slot "', expr_slot, '".')
    }
    if(length(marker_type_present) < length(marker_type)){
      warning('The following type marker could not be found in expr slot "', expr_slot, '": ',
              paste(marker_type[!marker_type %in% marker_type_present], collapse = ","))
    }
  }else{
    marker_type_present<-colnames(fcd$expr[[expr_slot]])
  }

  ## state marker
  if(!is.null(marker_state)){
    marker_state<-unique(marker_state)

    ##check if markers are present in expr slot
    marker_state_present<-marker_state[marker_state %in% colnames(fcd$expr[[expr_slot]])]
    if(length(marker_state_present) == 0){
      stop('None of the provided state marker are present in expr slot "', expr_slot, '".')
    }
    if(length(marker_state_present) < length(marker_state)){
      warning('The following state marker could not be found in expr slot "', expr_slot, '": ',
              paste(marker_state[!marker_state %in% marker_state_present], collapse = ","))
    }
  }else{
    marker_state_present<-NULL
  }

  ## prepare data frame marker_info
  marker_info<-data.frame(channel_name = colnames(fcd$expr[[expr_slot]]),
                          marker_name = colnames(fcd$expr[[expr_slot]]),
                          marker_class = c(rep("none",length(colnames(fcd$expr[[expr_slot]])))))
  marker_info[marker_info$marker_name %in% marker_type_present,]$marker_class<-"type"

  if(!is.null(marker_state)){
    marker_info[marker_info$marker_name %in% marker_state_present,]$marker_class<-"state"
  }


  #### prepare experiment_info meta data
  if(missing(meta_vars)){
    stop('meta_vars is not definded.
         At least one variable in cell_anno grouping the samples for testing musst be provided.')
  }

  meta_vars_present<-meta_vars[meta_vars %in% colnames(fcd$anno$cell_anno)]
  if(is.null(meta_vars_present) | length(meta_vars_present) < 1){
    stop('None of the variables in meta_vars was found in cell_anno.')
  }
  if(length(meta_vars_present) < length(meta_vars)){
    stop('The following variable listed in meta_vars could not be found in cell anno: ',
         paste(meta_vars[!meta_vars %in% meta_vars_present], collapse=","))
  }
  if(sum(meta_vars %in% c("sample_id")) > 0){
    stop('meta_vars cannot be named "sample_id",
         since this name will be given to the sample_var in the new object.')
  }
  if(sum(meta_vars %in% c("cluster_id")) > 0){
    stop('meta_vars cannot be named "cluster_id",
         since this name will be given to the cluster_var in the new object.')
  }

  experiment_info<-unique(fcd$anno$cell_anno[,c(sample_var,meta_vars)])
  names(experiment_info)[names(experiment_info) == sample_var] <- "sample_id"

  if(!nrow(experiment_info) == length(unique(experiment_info$sample_id))){
    stop('At least one ID in sample_var is associated with more than one level in at least one variable in meta_vars).
         Please ensure, that variables used in meta_vars have only one value for each sample identifier in sample_var.')
  }
  rownames(experiment_info)<-experiment_info$sample_id

  #### prepare n_cells
  n_cells <- table(fcd$anno$cell_anno[[sample_var]])
  n_cells <-n_cells[match(experiment_info$sample_id,names(n_cells))]


  #### prepare rowData
  rowData <- fcd$anno$cell_anno[,c(sample_var, meta_vars)]
  rowData$cluster_id <- fcd$clustering[[cluster_slot]][[cluster_var]]
  names(rowData)[names(rowData) == sample_var] <- "sample_id"


  #### prepare expression data
  data <- as.matrix(fcd$expr[[expr_slot]])


  #### prepare SummarizedExperiment object
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(exprs=data),
                                                   rowData = rowData,
                                                   colData = marker_info,
                                                   metadata = list(experiment_info = experiment_info,
                                                                   n_cells = n_cells))

  return(se)
}

#' marker_wilcox_test
#'
#' @title Differential expression testing on cell level
#' @description
#' EXPERIEMTAL FEATURE, we advise against using it yet and suggest you have a look at the Vignette "Differential Analysis" for other options.
#' \code{marker_wilcox_test()} performes a Wilcoxon Rank Sum Test for two groups of cells for each marker and cell population combination.
#' @param fcd flow cytometry data set, that has been subjected to the clustering or cell type label prediction with cyCONDOR
#' @param cluster_slot string specifying which clustering slot to use to find variable specified in cluster_var
#' @param cluster_var string specifying variable name in cluster_slot that identifies cell population labels to be used (e.g. clusters, metaclusters or predicted labels).
#' @param group_var string indicating variable name in cell_anno that should be used as grouping variable. The grouping variable needs to have two or more groupss
#' @param group1 string indicating group level in group_var that should be used to select cells for group 1
#' @param group2 string indicating group level in group_var that should be used to select cells for group 2
#' @param p.adjust.method p-value adjustment method to use for multiple comparison testing, e.g "BH" (Benjamini-Hochberg, default) or "bonferroni". All available options can be checked in the documentation of the \code{\link[rstatix]{adjust_pvalue}} function from the package \code{rstatix}.
#' @param marker (optional) vector of character strings indicating which features in the expression expr_slot should be considered during testing. by default, all features are tested.
#' @param min_cells_per_group Minimum number of cells per group required to include a cell population for differential testing.
#' @returns
#' A data frame containing test results for each marker and cell population combination - one combination per row.
#' \itemize{
#'  \item{cluster} : cell population that was tested
#'  \item{marker} : feature from expression data that was tested
#'  \item{group1} : group name (level) of group 1
#'  \item{group2} : group name (level) of group 2
#'  \item{n1 / n2} : absolute cell counts in group 1 (n1) and group 2 (n2)
#'  \item{mean1 / mean2} : mean marker expression in group 1 (mean1) and group 2 (mean2)
#'  \item{p} : p-value
#'  \item{p.adj} : adjusted p-value.
#'  \item{delta_mean} : delta of mean1 and mean2
#' }
#' @details
#' The function \code{marker_wilcox_test()} compares two groups of cells for each marker-cell population combination. Expression values will be extracted from expr_slot "orig", containing the transformed data.
#' In case group_var has more than two levels, the dataset will be subsetted to the two levels specified in group1 and group2.
#' Wilcoxon Rank Sum Test is performed using the \code{\link[stats]{wilcox.test}} implemented in the \code{stats} package. Afterwards p-value adjustment is performed considering all comparisons that were made.
#' @import tidyr
#' @import reshape2
#' @import dplyr
#' @import stats
#' @import rstatix
#'
#' @export
marker_wilcox_test<-function(fcd,
                             cluster_slot,
                             cluster_var,
                             group_var,
                             group1,
                             group2,
                             p.adjust.method = "BH",
                             marker = NULL,
                             min_cells_per_group = 10
){

  warning('HOLD UP! This is an experimental feature - we advise against using it yet and suggest you have a look at the Vignette "Differential Analysis" for other options.')

  expr_slot <- "orig"

  #### check slots, cellIDs und varibles
  checkInput(fcd = fcd,
             check_expr_slot = T,
             check_cluster_slot = T,
             check_cell_anno = T,
             expr_slot = expr_slot,
             cluster_slot = cluster_slot,
             cluster_var = cluster_var,
             group_var = group_var)

  #### get markers of interest
  if(!is.null(marker)){
    #remove potential duplicates
    marker <- unique(marker)

    ##check if markers are present in expr slot
    marker_present <- marker[marker %in% colnames(fcd$expr[[expr_slot]])]
    if(length(marker_present) == 0){
      stop('None of the provided markers are present in data slot "',expr_slot,'".')
    }
    if(length(marker_present)<length(marker)){
      warning('The following markers could not be found in data slot ',expr_slot,': ',
              paste(marker[!marker %in% marker_present], collapse = ","))
    }

  }else{
    marker_present <- colnames(fcd$expr[[expr_slot]])
  }


  #### prepare data
  data <- fcd$expr[[expr_slot]]
  data$cluster <- fcd$clustering[[cluster_slot]][[cluster_var]]
  data$group_var <- fcd$anno$cell_anno[[group_var]]

  if(is.numeric(data$cluster)){
    data$cluster <- as.character(data$cluster)}

  ##subset to marker of interest
  data <- as.data.frame(data[,c(marker_present,"cluster","group_var")])

  ##subset to groups of interested in case group_var has more than two entries
  groups <- unique(data$group_var)
  if(!group1 %in% groups){
    stop('group "',group1,'" could not be found in provided group_var.')
  }
  if(!group2 %in% groups){
    stop('group "',group2,'" could not be found in provided group_var.')
  }
  if(length(groups) > 2){
    data <- data[data$group_var %in% groups,]
  }

  ##subset to cluster of interest by cell counts per group
  clusters <- unique(as.character(data$cluster))

  cell_counts <- as.matrix(cyCONDOR::confusionMatrix(paste0(data$cluster),
                                           paste0(data$group_var)))

  cluster_oi <- rownames(cell_counts)[cell_counts[,1] >= min_cells_per_group & cell_counts[,2] >= min_cells_per_group]
  if(is.null(cluster_oi)){
    stop("No cluster passed the min_cells_per_group threshold.")
  }
  if(length(cluster_oi) < length(clusters)){
    warning("The following cluster(s) did not pass the min cells_per_group threshold: ",
            paste(clusters[!clusters %in% cluster_oi],collapse=", "),
            ". No test results will be reported.")
  }
  data <- data[data$cluster %in% cluster_oi,]


  ####calculate means per group x cluster combination
  data_stats <- data %>% dplyr::group_by(group_var, cluster, .drop=T) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)))
  data_stats <- data_stats %>% reshape2::melt(., id.vars=c("cluster", "group_var")) %>%
    tidyr::pivot_wider(names_from = group_var, values_from = value)
  names(data_stats)[names(data_stats) == group1] <- "mean1"
  names(data_stats)[names(data_stats) == group2] <- "mean2"
  names(data_stats)[names(data_stats) == "variable"] <- "marker"
  data_stats$delta_mean <- data_stats$mean1 - data_stats$mean2

  ##add cell counts
  cell_counts <- as.data.frame(cell_counts)
  cell_counts$cluster <- rownames(cell_counts)
  cell_counts <- cell_counts[cell_counts$cluster %in% cluster_oi,]
  names(cell_counts)[names(cell_counts) == group1] <- "n1"
  names(cell_counts)[names(cell_counts) == group2] <- "n2"
  data_stats <- dplyr::left_join(data_stats, cell_counts, by = c("cluster"))

  ####perform wilcoxon test group1 vs group2
  data$group_var <- factor(data$group_var, levels = c(group1, group2))

  wilcox_res <- data.frame()
  for(j in cluster_oi){
    data_cluster <- data[data$cluster == j,]

    for (i in marker_present){
      data_marker <- data_cluster[,c(i,"cluster","group_var")]
      names(data_marker)[1] <- "marker"

      res <- data.frame(cluster = j,
                       marker = i,
                       cluster_marker = paste(j,i))
      #wilcox test
      res$p <- stats::wilcox.test(data_marker$marker ~ data_marker$group_var)$p.value
      #res$estimate <- stats::wilcox.test(data_marker$marker ~ data_marker$group_var,conf.int = TRUE)$estimate

      wilcox_res <- rbind(wilcox_res,res)

      rm(res,data_marker)
    }
    print(paste("Calculations for cluster ",j," are done."))
  }

  ####perform p-value adjsutment
  wilcox_res <- rstatix::adjust_pvalue(data = wilcox_res, p.col = "p", method = p.adjust.method)


  #### add stats
  wilcox_res <- dplyr::left_join(wilcox_res, data_stats, by = c("cluster","marker"))
  #### add groups
  wilcox_res$group1 <- group1
  wilcox_res$group2 <- group2

  #### make results table nice
  wilcox_res <- wilcox_res[,c("cluster","marker","group1","group2","n1","n2","mean1","mean2",
                            "p","p.adj","delta_mean")]

  return(wilcox_res)
}

#' add_diffcyt_statistics
#'
#' @title Add statistical results from the diffcyt package to the fcd.
#' @description
#' Wrapper function to include the statistical results calculated using the diffcyt package into the fcd file.
#' @param fcd flow cytometry dataset, that has been subjected to clustering or cell type label prediction with cyCONDOR before.
#' @param input dataframe containing the statistical results from the diffcyt package. Example: \code{input <- as.data.frame(diffcyt::topTable(res_DA, all = TRUE))}.
#' @param group1 First group of the comparison.
#' @param group2 Second group of the comparison.
#' @import dplyr
#' @import diffcyt
#' @returns the fcd containing the statistical results from the diffcyt package
#'
#' @export
add_diffcyt_statistics <- function(fcd = condor, input, group1 = "ctrl", group2 = "pat") {

  ## check if the statistics slot exists. If not, initiate a list for storage.
  if(is.null(fcd$extras$statistics)){
    fcd$extras$statistics <- list()
  }

  ## add the results with matching colnames to the other statistical tests.
  fcd$extras$statistics$diffcyt <- input %>%
    mutate(
      group1 = group1,
      group2 = group2,

      ## add asteriks indications for significance - similar to the other statistical tests
      p.adj.signif = case_when(
        p_adj > 0.05     ~ "ns",
        p_adj <= 0.0001  ~ "****",
        p_adj <= 0.001   ~ "***",
        p_adj <= 0.01    ~ "**",
        p_adj <= 0.05    ~ "*"
      )
    ) %>%
    select(cluster = cluster_id, group1, group2, p = p_val, p.adj = p_adj, p.adj.signif)
  message("Statistics from diffcyt were saved to fcd$extras$statistics")
  return(fcd)
}
