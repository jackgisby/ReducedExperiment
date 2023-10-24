.single_lmer <- function(data, formula_string, REML = TRUE) {
    out_model <- tryCatch(
        lmerTest::lmer(
            as.formula(formula_string),
            data = data,
            REML = REML,
            control = lme4::lmerControl(check.conv.singular = "ignore")
        ),
        warning = function(w) {
            return(lmerTest::lmer(
                as.formula(formula_string),
                data = data,
                REML = REML,
                control = lme4::lmerControl(optimizer = "Nelder_Mead", check.conv.singular = "ignore")
            ))
        }
    )


    if (class(out_model) == "lmerModLmerTest") {
        return(out_model)
    } else {
        stop("Convergence issue not caught by single_lmer")
    }
}

.run_linear_model <- function(X, pheno, formula, method=c("lm", "lmer"), type="II", ...) {

    pheno$component <- X
    formula <- as.formula(paste0("component", formula))

    if (method == "lmer") {

        linear_model <- .single_lmer(pheno, formula, ...)

        anova_res <- data.frame(anova(linear_model, type=type))
        summary_res <- data.frame(summary(linear_model)$coefficients)

    } else {

        linear_model <- lm(paste0(colnames(sample_by_feature)[i], formula), combined_data)

        anova_res <- data.frame(car::Anova(linear_model, type=type))
        summary_res <- data.frame(summary(linear_model)$coefficients)
    }

    anova_res$term <- rownames(anova_res)
    summary_res$term <- rownames(summary_res)

    return(list("model"=model, "anova"=anova_res, "summary"=summary_res))
}

associate_factors <- function(re, formula, method=c("lm", "lmer"),
                              scale=TRUE, center=TRUE, return_models=FALSE,
                              type="II", ...) {

    models <- list()
    summaries <- anovas <- data.frame()

    for (comp in componentNames(re)) {

        linear_model <- .run_linear_model(
            X=reduced(re[[comp]], scale=scale, center=center),
            pheno=data.frame(colData(re)),
            formula=formula,
            method=method,
            type=type,
            ...
        )

        linear_model$anova$component <- linear_model$summary$component <- comp

        models[[comp]] <- linear_model$model
        anovas <- rbind(anovas, linear_model$anova)
        summaries <- rbind(summaries, linear_model$summary)
    }

    return(list("models"=models, "anovas"=anovas, "summaries"=summaries))
}
