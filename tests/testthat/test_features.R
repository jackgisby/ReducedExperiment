context("features")

test_that("getGeneIDs", {

    data(airway, package="airway")

    rrd <- .makeRandomData(ncol(airway), 10, "sample", "factor")
    rownames(rrd) <- colnames(airway)

    re <- ReducedExperiment(reduced=rrd, assays=assays(airway),
                            rowData=rowData(airway), colData=colData(airway),
                            metadata=metadata(airway))

    re_genes <- getGeneIDs(re)
})



test_that("FactorisedExperiment enrichment", {

    # Use real data from airway package
    airway <- .get_airway_data(n_features=2000)
    airway_fe <- estimate_factors(airway, nc=2, seed=1, use_stability=FALSE, method="imax")

    # Run overrepresentation analysis
    overrep_res <- runEnrich(airway_fe, method="overrepresentation", feature_id_col="gene_id", as_dataframe=TRUE, p_cutoff=0.1, universe=rownames(airway_fe))
    gsea_res <- runEnrich(airway_fe, method="gsea", feature_id_col="rownames", as_dataframe=TRUE, p_cutoff=0.1)

    expect_true("factor_1" %in% overrep_res$component & "factor_2" %in% overrep_res$component)
    expect_true("factor_1" %in% gsea_res$component & "factor_2" %in% gsea_res$component)
    expect_true(all(overrep_res$p.adjust < 0.1))
    expect_true(all(gsea_res$p.adjust < 0.1))
})

test_that("ModularExperiment enrichment", {

    # Use real data from airway package
    airway <- .get_airway_data(n_features=500)
    airway_me <- identify_modules(airway, verbose=0, powers=21)

    # Run overrepresentation analysis
    enrich_res <- runEnrich(airway_me, method="overrepresentation", as_dataframe=TRUE, p_cutoff=0.1)

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
    airway_fe <- estimate_factors(airway, nc=2, seed=1)

    cf <- get_common_features(getAlignedFeatures(airway_fe, z_cutoff=3, n_features=3, format="data.frame"))

    expect_equal(dim(cf), c(4, 7))
    expect_equal(cf$intersect, c(NA, 33, 33, NA))
    expect_equal(cf$total_feat_1, c(107, 107, 356, 356))
})
