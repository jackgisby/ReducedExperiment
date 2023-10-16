.ReducedSet <- setClass(
    "ReducedSet",
    contains = "ExpressionSet",
    representation = representation(reducedData = "matrix",
                                    L = "matrix"),
    prototype = prototype(new(
        "VersionedBiobase",
        versions = c(classVersion("ExpressionSet"),
                     ReducedSet = "0.1.0")
    ))
)

.FactorSet <- setClass("FactorSet",
                       contains = "ReducedSet",
                       prototype = prototype(new(
                           "VersionedBiobase",
                           versions = c(classVersion("ReducedSet"),
                                        FactorSet = "0.1.0")
                       )))

.ModuleSet <- setClass("ModuleSet",
                       contains = "ReducedSet",
                       prototype = prototype(new(
                           "VersionedBiobase",
                           versions = c(classVersion("ReducedSet"),
                                        ModuleSet = "0.1.0")
                       )))
