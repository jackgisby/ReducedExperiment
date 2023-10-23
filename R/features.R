#' Enrichment analysis for factors and modules
factor_enrich <- function(reduced_experiment, ensembl_id_col="ensembl_id", scale_loadings=TRUE) {

    enrich_res <- data.frame()

    for (comp in componentNames(reduced_experiment)) {



        enrich_res <- rbind(enrich_res, comp_enrich)
    }

}

module_enrich
