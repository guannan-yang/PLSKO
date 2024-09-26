tissue_express <- read_tsv("./data/rna_tissue_consensus.tsv")

tissue_express_w <- pivot_wider(tissue_express[,c("Gene", "Tissue", "nTPM")], names_from = "Tissue", values_from = "nTPM")

placenta.ratio <- apply(tissue_express_w[2:41], 1, function(x){x[24]/mean(x[-24])})

tissue_express_w$placenta_ratio <- placenta.ratio

tissue_express_w <- arrange(tissue_express_w, desc(placenta_ratio))

placenta_elevated2 <- filter(tissue_express_w, placenta_ratio >= 2)
length(intersect(placenta_elevated2$Gene, rownames(discv1_counts))) #257 for consensus; 289 for HPA

placenta_elevated3 <- filter(tissue_express_w, placenta_ratio >= 3)
length(intersect(placenta_elevated3$Gene, rownames(discv1_counts))) #96 for consensus; 105 for HPA

saveRDS(placenta_elevated3, "./data/placenta_3fold.rds")

placenta_elevated4 <- filter(tissue_express_w, placenta_ratio >= 4)
length(intersect(placenta_elevated4$Gene, rownames(discv1_counts))) # 61 for HPA
