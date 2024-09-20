#' Cell-free RNA-seq dataset with elevated expression in the placenta
#'
#' An example dataset containing cell-free RNA counts from 71 samples with 81 genes. The dataset is a list with 2 elements: 'counts' and 'metadata'. 'counts' is a matrix of cell-free RNA counts with 71 rows (samples) and 81 columns (genes). 'metadata' is a data frame with information about the samples, which were retrieved from GSE192902. The genes is a subset of genes with elevated expression in the placenta, according to The Human Protein Atlas database (v19): https://v19.proteinatlas.org/humanproteome/tissue/placenta.
#'
#' @format A list with 2 elements: 'counts' and 'metadata'. 'counts' is a matrix of cell-free RNA counts with 71 rows (samples) and 81 columns (genes). 'metadata' is a data frame with information about the samples.
#' @source The Human Protein Atlas database (v19): https://v19.proteinatlas.org/humanproteome/tissue/placenta
#' @source GSE192902
#'
#' @references Moufarrej MN, Vorperian SK, Wong RJ, Campos AA et al. Early prediction of preeclampsia in pregnancy with cell-free RNA. Nature 2022 Feb;602(7898):689-694. PMID: 35140405
"cfRNA_placenta"
