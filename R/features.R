#' Overrepresentation analysis
#'
#' Performs overrepresentation analysis given a set of components (factors or
#' modules) that are each linked to a list of genes.
#'
#' @noRd
#' @keywords internal
reduced_oa <- function(component_features, database = "msigdb_c2_cp",
    TERM2GENE = NULL, p_cutoff = 1, adj_method = "BH",
    min_genes = 3, universe = NULL, ...) {
    TERM2GENE <- .get_t2g(database, TERM2GENE)
    enrich_res <- list()

    for (comp in names(component_features)) {
        enrich_res_single <- clusterProfiler::enricher(
            component_features[[comp]],
            pvalueCutoff = 1, qvalueCutoff = 1,
            pAdjustMethod = adj_method,
            TERM2GENE = TERM2GENE,
            universe = universe,
            ...
        )

        enrich_res_single <- .format_enrich_res(
            enrich_res_single,
            adj_method = adj_method,
            min_genes = min_genes,
            p_cutoff = p_cutoff
        )

        if (!is.null(enrich_res_single)) {
            if (nrow(enrich_res_single@result) >= 1) {
                enrich_res_single@result$method <- "overrepresentation"
                enrich_res_single@result$component <- comp
            }
        }

        enrich_res[[comp]] <- enrich_res_single
    }

    return(enrich_res)
}

#' Gets a TERM2GENE data object for clusterProfiler
#'
#' Gets a correctly formatted TERM2GENE object for use by the `clusterProfiler`
#' package for enrichment analysis.
#'
#' @noRd
#' @keywords internal
.get_t2g <- function(database, TERM2GENE) {
    if (is.null(database)) {
        TERM2GENE <- TERM2GENE
    } else if (is.data.frame(database)) {
        TERM2GENE <- database
    } else if (database == "msigdb_c2_cp") {
        TERM2GENE <- get_msigdb_t2g()
    } else {
        stop("Database ", database, " not recognised")
    }

    return(TERM2GENE)
}

#' Reformats enrichment result
#'
#' Reformats `clusterProfiler` enrichment results into a consistent format.
#'
#' @noRd
#' @keywords internal
.format_enrich_res <- function(
    enrich_res_single,
    adj_method,
    p_cutoff,
    min_genes = NULL
) {
    if (is.null(enrich_res_single)) {
        return(NULL)
    }

    if (!is.null(min_genes)) {
        enrich_res_single@result <-
            enrich_res_single@result[which(enrich_res_single@result$Count
                                           >= min_genes), ]
    }

    enrich_res_single@result$p.adjust <-
        stats::p.adjust(enrich_res_single@result$pvalue, method = adj_method)

    enrich_res_single@result <-
        enrich_res_single@result[which(enrich_res_single@result$p.adjust
                                       <= p_cutoff), ]

    if (nrow(enrich_res_single@result) >= 1) {
        enrich_res_single@result$adj_method <- adj_method
        enrich_res_single@result$qvalue <- NA

        if ("pvalueCutoff" %in% slotNames(enrich_res_single)) {
            enrich_res_single@pvalueCutoff <- p_cutoff
        }
    }

    return(enrich_res_single)
}

#' Gene set enrichment analysis
#'
#' Performs \link[clusterProfiler]{GSEA} using the `clusterProfiler` package.
#' Expects a loadings matrix where features (genes) are rows and factors are
#' columns.
#'
#' @noRd
#' @keywords internal
reduced_gsea <- function(S, database = "msigdb_c2_cp", TERM2GENE = NULL,
    p_cutoff = 1, adj_method = "BH", nPermSimple = 1000, eps = 1e-10, ...) {
    TERM2GENE <- .get_t2g(database, TERM2GENE)
    enrich_res <- list()

    for (comp in colnames(S)) {
        S_order <- order(S[, comp], decreasing = TRUE)
        comp_genes <- S[, comp][order(S[, comp], decreasing = TRUE)]
        names(comp_genes) <- rownames(S)[order(S[, comp], decreasing = TRUE)]

        enrich_res_single <- clusterProfiler::GSEA(
            comp_genes,
            pvalueCutoff = 1,
            pAdjustMethod = adj_method,
            TERM2GENE = TERM2GENE,
            nPermSimple = nPermSimple,
            eps = eps,
            ...
        )

        enrich_res_single <- .format_enrich_res(
            enrich_res_single,
            adj_method = adj_method,
            p_cutoff = p_cutoff
        )

        if (!is.null(enrich_res_single)) {
            if (nrow(enrich_res_single@result) >= 1) {
                enrich_res_single@result$method <- "gsea"
                enrich_res_single@result$component <- comp
            }
        }

        enrich_res[[comp]] <- enrich_res_single
    }

    return(enrich_res)
}

#' Get TERM2GENE dataframe from MSigDB
#'
#' Gets C2 canonical pathways from the MSigDB database in the format required
#' by `clusterProfiler` enrichment functions.
#'
#' @param species The species for which to obtain MSigDB pathways. See
#' \link[msigdbr]{msigdbr} for more details.
#'
#' @param category The MSigDB category to retrieve pathways for. See
#' \link[msigdbr]{msigdbr} for more details.
#'
#' @param subcategory The MSigDB subcategory to retrieve pathways for. See
#' \link[msigdbr]{msigdbr} for more details.
#'
#' @param subcategory_to_remove If not NULL, this is a character string
#' indicating a subcategory to be removed from the results of
#' \link[msigdbr]{msigdbr}.
#'
#' @param gene_id The name to be given to the gene_id column of the resulting
#' data.frame.
#'
#' @returns Returns a data.frame, where the "gs_name" column indicates the name
#' of a pathway, and the `gene_id` column indicates genes that belong to
#' said pathway.
#'
#' @export
get_msigdb_t2g <- function(
    species = "Homo sapiens",
    category = "C2",
    subcategory = NULL,
    subcategory_to_remove = "CGP",
    gene_id = "ensembl_gene"
) {
    t2g <- data.frame(msigdbr::msigdbr(
        species = species,
        category = category,
        subcategory = subcategory
    ))

    if (!is.null(subcategory_to_remove)) {
        t2g <- t2g[which(t2g$gs_subcat != subcategory_to_remove), ]
    }

    t2g <- t2g[, which(colnames(t2g) %in% c("gs_name", gene_id))]

    return(t2g)
}

#' Get common factor features
#'
#' Function to count how many genes are aligned with multiple factors.
#'
#' @param factor_features The result of
#' \link[ReducedExperiment]{getAlignedFeatures}.
#'
#' @returns A data.frame for each factor pair with the numbers and proportions
#' of the genes in the input that overlap.
#'
#' @export
get_common_features <- function(factor_features) {
    common_features <- data.frame()

    for (c_1 in unique(factor_features$component)) {
        for (c_2 in unique(factor_features$component)) {
            total_feat_1 <- length(factor_features$feature[
                factor_features$component == c_1
            ])
            total_feat_2 <- length(factor_features$feature[
                factor_features$component == c_2
            ])
            smaller_total <- min(total_feat_1, total_feat_2)

            if (c_1 == c_2) {
                feat_intersect <- NA
            } else {
                feat_intersect <- length(intersect(
                    factor_features$feature[factor_features$component == c_1],
                    factor_features$feature[factor_features$component == c_2]
                ))
            }

            common_features_single <- data.frame(
                c_1 = c_1,
                c_2 = c_2,
                intersect = feat_intersect,
                total_feat_1 = length(factor_features$feature[
                    factor_features$component == c_1
                ]),
                total_feat_2 = length(factor_features$feature[
                    factor_features$component == c_2
                ])
            )

            common_features_single$smaller_total <-
                min(common_features_single$total_feat_1,
                    common_features_single$total_feat_2)

            common_features_single$intersect_prop <-
                common_features_single$intersect /
                common_features_single$smaller_total

            common_features <- rbind(common_features, common_features_single)
        }
    }

    return(common_features)
}

#' Heatmap comparing commonality across factors
#'
#' @param common_features The output of
#' \link[ReducedExperiment]{get_common_features}.
#'
#' @param filename The path at which to save the plot.
#'
#' @param color The colour palette to be used in the heatmap.
#'
#' @returns An object generated by \link[pheatmap]{pheatmap}.
#'
#' @export
plot_common_features <- function(
    common_features,
    filename = NA,
    color = grDevices::colorRampPalette(RColorBrewer::brewer.pal(
        n = 7,
        name = "YlOrRd")
    )(100)
) {
    common_features <- subset(common_features,
                              select = c("c_1", "c_2", "intersect_prop"))

    prop_mat <- stats::reshape(
        common_features,
        idvar = "c_1",
        v.names = "intersect_prop",
        timevar = "c_2",
        direction = "wide",
        sep = "_"
    )

    rownames(prop_mat) <- prop_mat$c_1
    prop_mat <- subset(prop_mat, select = -"c_1")
    colnames(prop_mat) <- gsub("intersect_prop_", "", colnames(prop_mat))

    max_abs <- max(abs(prop_mat), na.rm = TRUE)

    common_hmap <- pheatmap::pheatmap(
        prop_mat,
        na_col = "grey",
        filename = filename,
        color = color
    )

    return(common_hmap)
}


#' Get module preservation statistics
#'
#' Tests whether a set of modules defined in the reference dataset are
#' preserved in the test dataset.
#'
#' @param reference_dataset The dataset that was used to define the modules.
#' Must be a `data.frame` or `matrix` with features as rows and samples as
#' columns, or a \link[ReducedExperiment]{ModularExperiment} or
#' \link[SummarizedExperiment]{SummarizedExperiment} object.
#'
#' @param test_dataset The dataset that will be used to test for module
#' preservation. Must be a `data.frame` or `matrix` with features as rows and
#' samples as columns, or a \link[SummarizedExperiment]{SummarizedExperiment}
#' object. The features of `test_dataset` should be the same as
#' `reference_dataset`.
#'
#' @param reference_assay_name If the reference dataset is a
#' \link[ReducedExperiment]{ModularExperiment} or
#' \link[SummarizedExperiment]{SummarizedExperiment} object, this argument
#' specifies which assay slot was used to define the modules.
#'
#' @param test_assay_name If the reference dataset is a
#' \link[ReducedExperiment]{ModularExperiment} or
#' \link[SummarizedExperiment]{SummarizedExperiment} object, this argument
#' specifies which assay slot is to be used in preservation tests.
#'
#' @param module_assignments If the reference dataset is not a
#' \link[ReducedExperiment]{ModularExperiment} object, this argument is
#' necessary to specify the module assignments.
#'
#' @param greyName The name of the "module" of unassigned genes. Usually
#' "module_0" (ReducedExperiment) or "grey" (WGCNA).
#'
#' @param goldName The name to be used for the "gold" module (which is made up
#' of a random sample of all network genes). Passed to
#' \link[WGCNA]{modulePreservation}.
#'
#' @param networkType A string referring to the type of WGCNA network used for
#' the reference and test datasets. One of"unsigned", "signed" or
#' "signed hybrid". See \link[WGCNA]{adjacency}. Passed to
#' \link[WGCNA]{modulePreservation}.
#'
#' @param corFnc A string referring to the function to be used to calculate
#' correlation. One of "cor" or "bicor". Passed to
#' \link[WGCNA]{modulePreservation}.
#'
#' @param ... Additional arguments to be passed to
#' \link[WGCNA]{modulePreservation}.
#'
#' @returns A `data.frame` containing preservation statistics, as described
#' by \link[WGCNA]{modulePreservation}.
#'
#' @export
module_preservation <- function(reference_dataset, test_dataset,
    reference_assay_name = "normal",
    test_assay_name = "normal",
    module_assignments = NULL,
    greyName = "module_0",
    goldName = "random",
    networkType = "signed",
    corFnc = "cor",
    ...) {
    if (inherits(reference_dataset, "ModularExperiment")) {
        module_assignments <- assignments(reference_dataset)
    } else if (is.null(module_assignments)) {
        stop("If reference_dataset is not a ModularExperiment,
             module_assignments must not be NULL")
    }

    if (inherits(reference_dataset, "SummarizedExperiment")) {
        reference_dataset <- assay(reference_dataset, reference_assay_name)
    }

    if (inherits(test_dataset, "SummarizedExperiment")) {
        test_dataset <- assay(test_dataset, test_assay_name)
    }

    if (!identical(rownames(reference_dataset), rownames(test_dataset))) {
        stop("Rownames of reference_dataset do not match those of test_dataset")
    }

    multi_data <- list(
        "reference" = list("data" = t(reference_dataset)),
        "test" = list("data" = t(test_dataset))
    )

    return(WGCNA::modulePreservation(
        multi_data,
        list("reference" = stats::setNames(names(module_assignments),
                                           module_assignments)),
        dataIsExpr = TRUE,
        networkType = networkType,
        corFnc = corFnc,
        goldName = goldName,
        greyName = greyName,
        ...
    ))
}

#' Plot module preservation statistics
#'
#' @param module_preservation_results The output of
#' \link[ReducedExperiment]{module_preservation}
#'
#' @param show_random If TRUE, shows the random module in the plots.
#'
#' @param remove_module The name of a module to be hidden from the plots.
#'
#' @import ggplot2
#' @import patchwork
#' @export
plot_module_preservation <- function(
    module_preservation_results,
    show_random = TRUE,
    remove_module = NULL
) {
    mr_df <- module_preservation_results$preservation$observed$ref.reference$
        inColumnsAlsoPresentIn.test
    zs_df <- module_preservation_results$preservation$Z$ref.reference$
        inColumnsAlsoPresentIn.test

    mr_df$module <- rownames(mr_df)
    zs_df$module <- rownames(zs_df)

    if (show_random) {
        mr_gold <- mr_df$medianRank.pres[mr_df$module == "random"]
        zs_gold <- zs_df$Zsummary.pres[zs_df$module == "random"]
    } else {
        mr_df <- mr_df[which(mr_df$module != "random"), ]
        zs_df <- zs_df[which(zs_df$module != "random"), ]
    }

    if (!is.null(remove_module)) {
        mr_df <- mr_df[which(mr_df$module != remove_module), ]
        zs_df <- zs_df[which(zs_df$module != remove_module), ]
    }

    max_module_size <- max(mr_df$moduleSize)
    nudge_mr <- 0.14
    nudge_zs <- 0.4

    medianrank_plot <- ggplot(mr_df, aes_string("moduleSize", "medianRank.pres",
                                                col = "module")) +
        geom_point(size = 3) +
        geom_text(aes_string(label = "module"), col = "black",
                  nudge_y = nudge_mr, hjust = 0) +
        theme(legend.position = "none") +
        xlim(c(0, max_module_size * 1.3)) +
        expand_limits(y = 0)

    zsummary_plot <- ggplot(zs_df, aes_string("moduleSize", "Zsummary.pres",
                                              col = "module")) +
        geom_point(size = 3) +
        geom_text(aes_string(label = "module"), col = "black",
                  nudge_y = nudge_zs, hjust = 0) +
        theme(legend.position = "none") +
        xlim(c(0, max_module_size * 1.3)) +
        expand_limits(y = 0) +
        geom_hline(yintercept = 10, col = "green", linetype = "dashed") +
        geom_hline(yintercept = 2, col = "blue", linetype = "dashed")

    if (show_random) {
        medianrank_plot <- medianrank_plot + geom_hline(yintercept = mr_gold,
                                                        col = "gold",
                                                        linetype = "dashed")
        zsummary_plot <- zsummary_plot + geom_hline(yintercept = zs_gold,
                                                    col = "gold",
                                                    linetype = "dashed")
    }

    return(medianrank_plot + zsummary_plot)
}
