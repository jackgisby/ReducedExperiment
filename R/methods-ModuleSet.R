setMethod(
    "initialize",
    signature = "ModuleSet",
    definition = function(.Object, ...)
    {
        callNextMethod(.Object, ...)
    }
)

setMethod(ModuleSet, signature = c(
    assayData = "matrix",
    reducedData = "matrix",
    L = "matrix"
),
function(assayData,
         reducedData,
         L,
         phenoData = annotatedDataFrameFrom(assayData, byrow =
                                                FALSE),
         featureData = annotatedDataFrameFrom(assayData, byrow =
                                                  TRUE),
         experimentData = MIAME(),
         annotation = character(),
         protocolData = annotatedDataFrameFrom(assayData, byrow =
                                                   FALSE),
         ...)
{
    assayData <- assayDataNew(exprs = assayData)
    reducedData <- assayDataNew(reduced = reducedData)
    .ModuleSet(
        assayData = assayData,
        reducedData = reducedData,
        L = L,
        phenoData = phenoData,
        featureData = featureData,
        experimentData = experimentData,
        annotation = annotation,
        protocolData = protocolData,
        ...
    )
})
