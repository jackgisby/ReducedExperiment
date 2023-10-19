#' ReducedSet
ReducedSet <- setClass(
    "ReducedSet",
    contains="eSet",
    representation=representation(reduced="matrix", exprs="matrix"),
    prototype=prototype(
        new("VersionedBiobase", versions=c(classVersion("eSet"), ReducedSet="0.1.0")),
        reduced=matrix(), exprs=matrix())
)

#' A class to represent the results of factor analysis
#' @export
FactorSet <- setClass(
    "FactorSet",
    contains="ReducedSet",
    representation=representation(S="matrix", varexp="numeric"),
    prototype=prototype(new(
        "VersionedBiobase",
        versions=c(classVersion("ReducedSet"), FactorSet="0.1.0")
    ), S=matrix(), varexp=numeric())
)

#' A class to represent modules
#' @export
ModuleSet <- setClass(
    "ModuleSet",
    contains="ReducedSet",
    representation=representation(assignments="numeric"),
    prototype=prototype(new(
        "VersionedBiobase",
        versions=c(classVersion("ReducedSet"), ModuleSet="0.1.0")
    ), assignments=numeric())
)
