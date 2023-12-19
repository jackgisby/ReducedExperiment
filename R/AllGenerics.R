
setGeneric("reduced", function(x, ...) standardGeneric("reduced"))
setGeneric("reduced<-", function(x, value) standardGeneric("reduced<-"))

setGeneric("loadings", function(x, ...) standardGeneric("loadings"))
setGeneric("loadings<-", function(x, value) standardGeneric("loadings<-"))

setGeneric("stability", function(x) standardGeneric("stability"))
setGeneric("stability<-", function(x, value) standardGeneric("stability<-"))

setGeneric("componentNames", function(x) standardGeneric("componentNames"))
setGeneric("componentNames<-", function(x, value) standardGeneric("componentNames<-"))

setGeneric("moduleNames",       function(x, ...) standardGeneric("moduleNames"))
setGeneric("moduleNames<-",     function(x, value) standardGeneric("moduleNames<-"))

setGeneric("sampleNames", function(x) standardGeneric("sampleNames"))
setGeneric("sampleNames<-", function(x, value) standardGeneric("sampleNames<-"))

setGeneric("featureNames", function(x) standardGeneric("featureNames"))
setGeneric("featureNames<-", function(x, value) standardGeneric("featureNames<-"))

setGeneric("projectData", function(x, newdata, ...) standardGeneric("projectData"))
setGeneric("calcEigengenes", function(x, newdata, ...) standardGeneric("calcEigengenes"))

setGeneric("assignments",       function(x, ...) standardGeneric("assignments"))
setGeneric("assignments<-",     function(x, value) standardGeneric("assignments<-"))

setGeneric("nComponents",     function(x) standardGeneric("nComponents"))
setGeneric("nModules",     function(x) standardGeneric("nModules"))
setGeneric("nSamples",     function(x) standardGeneric("nSamples"))
setGeneric("nFeatures",     function(x) standardGeneric("nFeatures"))

setGeneric("getGeneIDs",     function(x, ...) standardGeneric("getGeneIDs"))

setGeneric("getAlignedFeatures",     function(x, ...) standardGeneric("getAlignedFeatures"))

setGeneric("runEnrich",     function(x, ...) standardGeneric("runEnrich"))

setGeneric("plotDendro",     function(x, ...) standardGeneric("plotDendro"))

setGeneric("dendrogram",     function(x) standardGeneric("dendrogram"))
setGeneric("dendrogram<-",     function(x, value) standardGeneric("dendrogram<-"))
