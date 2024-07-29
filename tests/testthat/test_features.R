context("features")

test_that("FactorisedExperiment enrichment", {
    # Use real data from airway package

    set.seed(2)
    airway <- .get_airway_data(n_features = 2000)

    set.seed(1)
    airway_fe <- estimate_factors(airway, nc = 2, use_stability = FALSE, method = "imax")

    # Run overrepresentation analysis
    overrep_res <- runEnrich(airway_fe, method = "overrepresentation", feature_id_col = "rownames", as_dataframe = TRUE, p_cutoff = 0.1, universe = rownames(airway_fe))
    gsea_res <- runEnrich(airway_fe, method = "gsea", feature_id_col = "rownames", as_dataframe = TRUE, p_cutoff = 0.1)

    expect_true("factor_1" %in% overrep_res$component & "factor_2" %in% overrep_res$component)
    expect_true("factor_2" %in% gsea_res$component)
    expect_true(all(overrep_res$p.adjust < 0.1))
    expect_true(all(gsea_res$p.adjust < 0.1))
})

test_that("ModularExperiment enrichment", {
    # Use real data from airway package

    set.seed(2)
    airway <- .get_airway_data(n_features = 500)
    airway_me <- identify_modules(airway, verbose = 0, powers = 21)

    # Run overrepresentation analysis
    enrich_res <- runEnrich(airway_me, method = "overrepresentation", as_dataframe = TRUE, p_cutoff = 0.1)

    expect_true(all(paste0("module_", c(1, 2, 4, 5)) %in% enrich_res$component))
    expect_true(all(enrich_res$p.adjust < 0.1))
})

test_that("Get MSGIDB data", {
    t2g <- get_msigdb_t2g()

    expect_equal(ncol(t2g), 2)
    expect_true(nrow(t2g) > 100000)
    expect_true(all(grepl("ENSG", t2g$ensembl_gene) | grepl("ASMPATCH", t2g$ensembl_gene)))
    for (id in c("BIOCARTA", "KEGG", "REACTOME")) expect_true(any(grepl(id, t2g$gs_name)))
})

test_that("Get common features", {
    airway <- .get_airway_data()

    set.seed(1)
    airway_fe <- estimate_factors(airway, nc = 2)

    cf <- get_common_features(getAlignedFeatures(airway_fe, format = "data.frame"))

    expect_equal(dim(cf), c(4, 7))
    expect_equal(cf$intersect, c(NA, 8, 8, NA))
    expect_equal(cf$total_feat_1, c(63, 63, 38, 38))
})

test_that("Module preservation", {

    set.seed(2)
    airway <- .get_airway_data(n_features = 500)
    airway_me <- identify_modules(airway, verbose = 0, powers = 21)

    assay(airway_me, "noised") <- assay(airway_me, "transformed") + matrix(rnorm(nrow(airway_me) * ncol(airway_me), mean = 0, sd = 0.3), nrow = nrow(airway_me), ncol = ncol(airway_me))

    # Test module preservation
    mp <- module_preservation(airway_me, airway_me, reference_assay_name = "transformed", test_assay_name = "noised", verbose = 10, nPermutations = 10)
    plot_module_preservation(mp)

    # Test that module_preservation works with matrices
    mp_matrices <- module_preservation(assay(airway_me, "transformed"), assay(airway_me, "noised"), module_assignments = assignments(airway_me), verbose = 10, nPermutations = 10)
    plot_module_preservation(mp)

    # Ensure results are the same
    expect_equal(
        mp$preservation$Z$ref.reference$inColumnsAlsoPresentIn.test,
        mp_matrices$preservation$Z$ref.reference$inColumnsAlsoPresentIn.test
    )

    # Ensure features are the same
    expect_error(module_preservation(airway_me, airway_me[1:10, ]))
})
