setGeneric(
    "reduced",
    function(object, ...) standardGeneric("reduced")
)
setGeneric(
    "reduced<-",
    function(object, value) standardGeneric("reduced<-")
)

setGeneric(
    "loadings",
    function(object, ...) standardGeneric("loadings")
)
setGeneric(
    "loadings<-",
    function(object, value) standardGeneric("loadings<-")
)

setGeneric(
    "stability",
    function(object) standardGeneric("stability")
)
setGeneric(
    "stability<-",
    function(object, value) standardGeneric("stability<-")
)

setGeneric(
    "componentNames",
    function(object) standardGeneric("componentNames")
)
setGeneric(
    "componentNames<-",
    function(object, value) standardGeneric("componentNames<-")
)

setGeneric(
    "moduleNames",
    function(object, ...) standardGeneric("moduleNames")
)
setGeneric(
    "moduleNames<-",
    function(object, value) standardGeneric("moduleNames<-")
)

setGeneric(
    "sampleNames",
    function(object) standardGeneric("sampleNames")
)
setGeneric(
    "sampleNames<-",
    function(object, value) standardGeneric("sampleNames<-")
)

setGeneric(
    "featureNames",
    function(x) standardGeneric("featureNames")
)
setGeneric(
    "featureNames<-",
    function(x, value) standardGeneric("featureNames<-")
)

setGeneric(
    "projectData",
    function(object, newdata, ...) standardGeneric("projectData")
)
setGeneric(
    "calcEigengenes",
    function(object, newdata, ...) standardGeneric("calcEigengenes")
)

setGeneric(
    "assignments",
    function(object, ...) standardGeneric("assignments")
)
setGeneric(
    "assignments<-",
    function(object, value) standardGeneric("assignments<-")
)

setGeneric(
    "nComponents",
    function(object) standardGeneric("nComponents")
)
setGeneric(
    "nModules",
    function(object) standardGeneric("nModules")
)
setGeneric(
    "nSamples",
    function(object) standardGeneric("nSamples")
)
setGeneric(
    "nFeatures",
    function(object) standardGeneric("nFeatures")
)

setGeneric(
    "getGeneIDs",
    function(object, ...) standardGeneric("getGeneIDs")
)

setGeneric(
    "getAlignedFeatures",
    function(object, ...) standardGeneric("getAlignedFeatures")
)
setGeneric(
    "getCentrality",
    function(object, ...) standardGeneric("getCentrality")
)

setGeneric(
    "runEnrich",
    function(object, ...) standardGeneric("runEnrich")
)

setGeneric(
    "plotDendro",
    function(object, ...) standardGeneric("plotDendro")
)

setGeneric(
    "dendrogram",
    function(object) standardGeneric("dendrogram")
)
setGeneric(
    "dendrogram<-",
    function(object, value) standardGeneric("dendrogram<-")
)
