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

    # Example expression data formatted as a SummarizedExperiment
    data(airway, package="airway")

    # Make random reduced and loadings data matching airway dimensions
    rrd <- .makeRandomData(ncol(airway), 10, "sample", "factor")
    rownames(rrd) <- colnames(airway)
    rld <- .makeRandomData(nrow(airway), 10, "gene", "factor")
    rownames(rld) <- rownames(airway)

    # Create the FactorExperiment and get entrez gene IDs
    fe <- .se_to_fe(airway, reduced=rrd, loadings=rld, stability=numeric(), scale=FALSE, center=FALSE)

    fe <- getGeneIDs(fe)
    fe <- fe[which(!is.na(rowData(fe)$entrezgene_id)) ,]

    # Run overrepresentation analysis
    enrich_res <- runEnrich(fe, method="overrepresentation", feature_id_col="ensembl_gene_id", as_dataframe=TRUE, p_cutoff=0.01, universe=rownames(fe))
    #TODO: Test enrich result?
})

test_that("ModularExperiment enrichment", {

})

test_that("Get MSGIDB data", {

})

test_that("Get common features", {

})
