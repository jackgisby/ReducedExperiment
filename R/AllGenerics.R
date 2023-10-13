setGeneric("reducedData",   function(object) standardGeneric("reducedData"))
setGeneric("reducedData<-", function(object, value) standardGeneric("reducedData<-"))
setGeneric("write.reduced", function(x,...) standardGeneric("write.reduced"))
setGeneric("write.exprs",     function(x,...) standardGeneric("write.exprs"))
setGeneric("exprs",           function(object) standardGeneric("exprs"))
setGeneric("exprs<-",         function(object, value) standardGeneric("exprs<-"))
setGeneric("reduced",       function(object) standardGeneric("reduced"))
setGeneric("reduced<-",     function(object, value) standardGeneric("reduced<-"))
setGeneric("loadings",       function(object) standardGeneric("loadings"))
setGeneric("loadings<-",     function(object, value) standardGeneric("loadings<-"))

setGeneric("FactorSet",
           function(assayData,
                    reducedData,
                    L,
                    phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                    featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                    experimentData=MIAME(),
                    annotation=character(),
                    protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                    ...)
             standardGeneric("FactorSet"),
           signature=c("assayData", "reducedData", "L"))

setGeneric("ModuleSet",
           function(assayData,
                    reducedData,
                    L,
                    phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                    featureData=annotatedDataFrameFrom(assayData, byrow=TRUE),
                    experimentData=MIAME(),
                    annotation=character(),
                    protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE),
                    ...)
             standardGeneric("ModuleSet"),
           signature=c("assayData", "reducedData", "L"))
