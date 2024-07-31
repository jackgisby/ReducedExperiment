library(ReducedExperiment)

airway <- ReducedExperiment:::.get_airway_data()

dataset = "hsapiens_gene_ensembl"
ids_to_get = c("hgnc_symbol", "entrezgene_id")
gene_ids <- rownames(airway)
gene_id_type = "ensembl_gene_id"

mart <- biomaRt::useEnsembl(biomart = "genes", dataset = dataset)

biomart_out <- biomaRt::getBM(
    filters = gene_id_type,
    attributes = c(gene_id_type, ids_to_get),
    values = gene_ids, mart = mart
)

head(biomart_out)
write.csv(biomart_out, "inst/extdata/biomart_out.csv", row.names = FALSE)
