#' Heatmap of multiple variable selections ordered by importance
#'
#' @param S list of variable selection indices
#' @param labels labels of actual variable names corresponding to selection indices.
#' @param nbcocluster bivariate vector c(number of variable clusters, number of selection clusters).
#' The former number must be specified less than length(labels) and the latter must be less than length(S).
#'
#' @details The list of selected variables is converted to a length(labels)-by-length(S) binary
#' heatmap where each entry is either 1 if variable is selected and 0 otherwise. To help visualize most
#' important variables we perform coclustering of both selections and variables (with the blockcluster package).
#'
#' @return plot of heatmap
#' @export
#'
#' @examples
#' library(seqknockoff)
#'
#' set.seed(1)
#'
#' S <- list(c(sample(1:10,7), sample(11:20,1)),
#' c(sample(1:10,7), sample(11:20,1)),
#' c(sample(1:10,7), sample(11:20,1)))
#'
#' plot_heatmap(S, labels=1:20, nbcocluster=c(5,3))
plot_heatmap <- function(S, labels, nbcocluster=c(5,5)) {

  zeros <- rep(0, length(labels))
  selections <- data.frame(draw=factor(rep(rep(1:length(S)),each=length(labels))),
                           variable=factor(labels),
                           selected=unlist(lapply(S, function(x) {zeros[x] <- 1; return(zeros)})))

  if (!requireNamespace("dplyr", quietly = TRUE) | !requireNamespace("blockcluster", quietly = TRUE)) {
    warning("Packages \"dplyr\" and \"blockcluster\" are needed to order variables of heatmap by importance. Please install them.",
         call. = FALSE)
  } else {

  `%>%` <- dplyr::`%>%`

  biclust <- blockcluster::coclusterBinary(matrix(selections$selected,nrow=length(labels)), nbcocluster=nbcocluster)

  selections$varclass = factor(biclust@rowclass, ordered=TRUE)
  selections$drawclass = factor(rep(biclust@colclass, each=length(labels)), ordered=TRUE)

  # Calculate means per draw cluster (to order heatmap)
  meta.draw <- selections %>%
    dplyr::group_by(drawclass) %>%
    dplyr::summarise(selected = mean(selected)) %>%
    dplyr::arrange(desc(selected))

  # Calculate means per variable block (to order heatmap)
  meta.var <- selections %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(selected = mean(selected))%>%
    dplyr::arrange(selected)

  # Order selections according to block cluster means:
  selections <- selections %>%
    dplyr::mutate(drawclass = factor(drawclass, levels = as.character(meta.draw$drawclass)),
           variable = factor(variable, levels = as.character(meta.var$variable))) %>%
    dplyr::arrange(drawclass) %>%
    dplyr::mutate(draw = factor(draw, levels = unique(draw))) %>%
    dplyr::arrange(variable)
  }

  selections$selected <- factor(selections$selected)

  # ggplot2::ggplot(mapping=ggplot2::aes(x = draw, y = variable)) +
  #   ggplot2::geom_tile(ggplot2::aes(fill = selected), data=selections) +
  #   ggplot2::xlab("Repetition") +
  #   ggplot2::theme(axis.text.x = ggplot2::element_blank(),
  #         axis.ticks.x = ggplot2::element_blank()) +
  # ggplot2::scale_fill_manual(values=c("0"="#132B43","1"="#56B1F7"))
  
  draw_select <- selections %>% 
    dplyr::filter(variable %in% meta.var$variable[meta.var$selected!=0])# varclass 4 is those who never been selected
  draw_freq <- meta.var %>% 
    dplyr::filter(variable %in% draw_select$variable) %>% 
    gridExtra::tableGrob()
  
  ggplot2::ggplot(mapping=ggplot2::aes(x = draw, y = variable)) +
    ggplot2::geom_tile(ggplot2::aes(fill = selected), data=draw_select) +
    ggplot2::xlab("Repetition") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank()) +
    ggplot2::scale_fill_manual(values=c("0"="#132B43","1"="#56B1F7"))#+
    #ggplot2::annotation_custom(draw_freq, xmin = length(S)+2)

}
