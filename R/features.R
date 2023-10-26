#' Overrepresentation analysis
reduced_oa <- function(component_features, database="msigdb_c2_cp", TERM2GENE=NULL,
                                       p_cutoff=1, adj_method="BH",
                                       min_genes=3, universe=NULL, ...) {

    if (is.null(database)) {
        TERM2GENE <- TERM2GENE
    } else if (is.data.frame(database)) {
        TERM2GENE <- database
    } else if (database == "msigdb_c2_cp") {
        TERM2GENE <- get_msigdb_t2g()
    } else {
        stop("Database ", database, " not recognised")
    }

    enrich_res <- data.frame()

    for (comp in names(component_features)) {

        enrich_res_single <- clusterProfiler::enricher(
            component_features[[comp]],
            pvalueCutoff=1, qvalueCutoff=1,
            pAdjustMethod=adj_method,
            TERM2GENE=TERM2GENE,
            universe=universe,
            ...
        )

        if (is.null(enrich_res_single)) next

        enrich_res_single <- enrich_res_single@result
        enrich_res_single$adj_pvalue <- p.adjust(enrich_res_single$pvalue, method=adj_method)

        enrich_res_single <- enrich_res_single[which(enrich_res_single$Count >= min_genes) ,]
        enrich_res_single <- enrich_res_single[which(enrich_res_single$adj_pvalue < p_cutoff) ,]

        if (nrow(enrich_res_single) >= 1) {
            enrich_res_single$component <- comp
            enrich_res_single$method <- "overrepresentation"
            enrich_res_single$adj_method <- adj_method

            keep_cols <- c("ID", "Description", "component", "GeneRatio",
                "BgRatio", "pvalue", "method", "adj_pvalue", "adj_method",
                "Count", "geneID")

            enrich_res_single <- subset(enrich_res_single, select=keep_cols)

            colnames(enrich_res_single) <- c("ID", "description", "component",
                "gene_ratio", "bg_ratio", "pvalue", "method", "adj_pvalue",
                "adj_method", "count", "gene_id")

            enrich_res <- rbind(enrich_res, enrich_res_single)
        }
    }

    return(enrich_res)
}

#' Gene set enrichment analysis
reduced_gsea <- function() {
    # TODO: Implement
    stop("Not implemented yet")
}

#' Get TERM2GENE dataframe from MSigDB
get_msigdb_t2g <- function(species="Homo sapiens", category="C2", subcategory=NULL, subcategory_to_remove="CGP", gene_id="ensembl_gene") {

    t2g <- data.frame(msigdbr::msigdbr(species=species, category=category, subcategory=subcategory))
    
    if (!is.null(subcategory_to_remove)) {
        t2g <- t2g[which(t2g$gs_subcat != subcategory_to_remove),]
    }
    
    t2g <- t2g[, which(colnames(t2g) %in% c("gs_name", gene_id))]

    return(t2g)
}
