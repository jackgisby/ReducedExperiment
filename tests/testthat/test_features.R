context("features")

test_that("FactorisedExperiment enrichment and common features", {
    # Use real data from airway package

    set.seed(2)
    airway <- ReducedExperiment:::.get_airway_data(n_features = 2000)

    set.seed(1)
    airway_fe <- estimate_factors(airway, nc = 2, use_stability = FALSE, method = "imax")

    # Run overrepresentation analysis
    t2g <- read.csv(system.file(
        "extdata",
        "msigdb_t2g_filtered.csv",
        package = "ReducedExperiment"
    ))

    p_cutoff <- 0.9999

    overrep_res <- runEnrich(airway_fe, method = "overrepresentation", feature_id_col = "rownames", as_dataframe = TRUE, p_cutoff = p_cutoff, universe = rownames(airway_fe), TERM2GENE = t2g)
    gsea_res <- runEnrich(airway_fe, method = "gsea", feature_id_col = "rownames", as_dataframe = TRUE, p_cutoff = p_cutoff, TERM2GENE = t2g)

    expect_true("factor_1" %in% overrep_res$component & "factor_2" %in% overrep_res$component)
    expect_true("factor_2" %in% gsea_res$component)
    expect_true(all(overrep_res$p.adjust < p_cutoff))
    expect_true(all(gsea_res$p.adjust < p_cutoff))

    set.seed(1)
    airway_fe <- estimate_factors(airway, nc = 3, use_stability = FALSE, method = "imax")
    cf <- get_common_features(getAlignedFeatures(airway_fe, format = "data.frame"))

    expect_equal(dim(cf), c(9, 7))
    expect_equal(cf$intersect, c(NA, 5, 1, 5, NA, 1, 1, 1, NA))
    expect_equal(cf$total_feat_1, c(20, 20, 20, 20, 20, 20, 16, 16, 16))

    plot_common_features(cf)
})

test_that("ModularExperiment enrichment and preservation", {

    # Use real data from airway package with random modules
    set.seed(2)
    airway <- ReducedExperiment:::.get_airway_data(n_features = 500)
    airway_me <- ReducedExperiment:::.createRandomisedModularExperiment(dim(airway)[1], dim(airway)[2], 4)

    colnames(airway_me) <- colnames(airway)
    rownames(airway_me) <- rownames(airway)
    assay(airway_me, "normal") <- assay(airway, "normal")

    # Run overrepresentation analysis
    t2g <- read.csv(system.file(
        "extdata",
        "msigdb_t2g_filtered.csv",
        package = "ReducedExperiment"
    ))

    p_cutoff <- 0.9999
    enrich_res <- runEnrich(
        airway_me,
        method = "overrepresentation",
        as_dataframe = TRUE,
        p_cutoff = p_cutoff,
        TERM2GENE = t2g
    )

    expect_true(all(paste0("module_", c(1, 3)) %in% enrich_res$component))
    expect_true(all(enrich_res$p.adjust < p_cutoff))

    # Setup for preservation test
    airway_me <- airway_me[names(assignments(airway_me)) %in% c("module_1", "module_4"), ]
    assay(airway_me, "noised") <- assay(airway_me, "normal") + matrix(rnorm(nrow(airway_me) * ncol(airway_me), mean = 0, sd = 0.3), nrow = nrow(airway_me), ncol = ncol(airway_me))

    # Test module preservation
    mp <- module_preservation(airway_me, airway_me, reference_assay_name = "normal", test_assay_name = "noised", verbose = 0, nPermutations = 2)
    plot_module_preservation(mp)

    expect_equal(mp$preservation$Z$ref.reference$inColumnsAlsoPresentIn.test$Zsummary.pres, c(6.17, 13.28, 15.55), tolerance = 0.01)
})

test_that("Get MSGIDB data", {
    t2g <- get_msigdb_t2g()

    expect_equal(ncol(t2g), 2)
    expect_true(nrow(t2g) > 100000)
    expect_true(all(grepl("ENSG", t2g$ensembl_gene) | grepl("ASMPATCH", t2g$ensembl_gene)))
    for (id in c("BIOCARTA", "KEGG", "REACTOME")) expect_true(any(grepl(id, t2g$gs_name)))
})

