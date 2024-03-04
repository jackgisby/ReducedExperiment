#' @importFrom lmerTest lmer
#' @importFrom lme4 lmerControl
.single_lmer <- function(data, formula_string, REML = TRUE) {
    out_model <- tryCatch(
        lmerTest::lmer(
            stats::as.formula(formula_string),
            data = data,
            REML = REML,
            control = lme4::lmerControl(check.conv.singular = "ignore")
        ),
        warning = function(w) {
            return(lmerTest::lmer(
                stats::as.formula(formula_string),
                data = data,
                REML = REML,
                control = lme4::lmerControl(optimizer = "Nelder_Mead", check.conv.singular = "ignore")
            ))
        }
    )


    if (is(out_model, "lmerModLmerTest")) {
        return(out_model)
    } else {
        stop("Convergence issue not caught by single_lmer")
    }
}

.run_linear_model <- function(X, pheno, formula, method="lm", type="II", ...) {

    pheno$component <- X

    formula <- stats::as.formula(paste0("component", formula))

    if (method == "lmer") {

        linear_model <- .single_lmer(pheno, formula, ...)

        anova_res <- data.frame(stats::anova(linear_model, type=type))
        summary_res <- data.frame(summary(linear_model)$coefficients)

    } else if (method == "lm") {

        linear_model <- stats::lm(formula, pheno)

        anova_res <- data.frame(car::Anova(linear_model, type=type))
        summary_res <- data.frame(summary(linear_model)$coefficients)
    }

    anova_res$term <- rownames(anova_res)
    summary_res$term <- rownames(summary_res)

    return(list("model"=linear_model, "anova"=anova_res, "summary"=summary_res))
}

associate_components <- function(re, formula, method="lm",
                              scale_reduced=TRUE, center_reduced=TRUE,
                              type="II", adj_method="BH", ...) {

    models <- list()
    summaries <- anovas <- data.frame()
    red <- reduced(re, scale_reduced=scale_reduced, center_reduced=center_reduced)

    for (comp in componentNames(re)) {

        linear_model <- .run_linear_model(
            X=red[, comp],
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

    colnames(anovas) <- .rename_results_table(colnames(anovas))
    colnames(summaries) <- .rename_results_table(colnames(summaries))

    rownames(anovas) <- paste(anovas$component, anovas$term, sep="_")
    rownames(summaries) <- paste(summaries$component, summaries$term, sep="_")

    anovas$adj_pvalue <- .adjust_by_term(anovas, method=adj_method)
    summaries$adj_pvalue <- .adjust_by_term(summaries, method=adj_method)

    return(list("models"=models, "anovas"=anovas, "summaries"=summaries))
}

.adjust_by_term <- function(res, method="BH") {

    res$adj_pvalue <- NA

    for (term in unique(res$term)) {
        res$adj_pvalue[which(res$term == term)] <- stats::p.adjust(res$pvalue[which(res$term == term)], method=method)
    }

    return(res$adj_pvalue)
}

.rename_results_table <- function(cnames) {
    cname_conversions <- list(
        "Sum.Sq" = "sum_sq",
        "Mean.Sq" = "mean_sq",
        "NumDF" = "num_df",
        "DenDF" = "den_df",
        "F.value" = "fvalue",
        "Pr..F." = "pvalue",
        "Estimate" = "estimate",
        "Std..Error" = "stderr",
        "t.value" = "tvalue",
        "Pr...t.." = "pvalue"
    )

    for (i in seq_along(cnames))
        if (cnames[i] %in% names(cname_conversions))
            cnames[i] <- cname_conversions[[cnames[i]]]

    return(cnames)
}
