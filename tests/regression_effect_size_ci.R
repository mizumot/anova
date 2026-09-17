args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) {
    stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg))
app_dir <- normalizePath(file.path(dirname(script_path), ".."))

source(file.path(app_dir, "anovakun_489.R"), encoding = "UTF-8")
source(file.path(app_dir, "validation.R"), encoding = "UTF-8")
source(file.path(app_dir, "presentation.R"), encoding = "UTF-8")

example_values <- list()
collect_examples <- function(node) {
    if (is.call(node) && identical(as.character(node[[1]]), "aceEditor")) {
        id <- eval(node[[2]])
        example_values[[id]] <<- eval(node$value)
    }
    if (is.recursive(node)) {
        for (child in as.list(node)[-1]) {
            collect_examples(child)
        }
    }
}

for (expression in parse(file.path(app_dir, "ui.R"))) {
    collect_examples(expression)
}

spec <- make_anova_spec(
    "twoway", NULL, "Factor1Between_Factor2Within", 2, 3
)
checked <- validate_anova_input(example_values$text4, spec)
stopifnot(checked$valid)
labels <- derive_analysis_labels(
    checked$data,
    spec,
    factor1_name = "Method",
    factor2_name = "Time",
    outcome_name = "Score"
)
set.seed(20260917)
result <- do.call(
    anovakun,
    c(
        list(dataset = checked$data, design = spec$design),
        labels$factor_args,
        list(
            mau = TRUE, auto = TRUE, holm = TRUE,
            peta = TRUE, geta = TRUE, besci = TRUE, tech = TRUE
        )
    )
)

ci_table <- result[["EFFECT SIZE INFORMATION"]][[2]]
required_columns <- c("ES", "Source", "Observed", "CI_L", "CI_U")

stopifnot(
    all(required_columns %in% names(ci_table)),
    identical(unique(ci_table$ES), c("p.eta^2", "G.eta^2")),
    all(is.finite(ci_table$Observed)),
    all(is.finite(ci_table$CI_L)),
    all(is.finite(ci_table$CI_U)),
    all(ci_table$CI_L <= ci_table$CI_U),
    all(ci_table$Observed >= 0 & ci_table$Observed <= 1),
    identical(unique(ci_table$Source), c("Method", "Time", "Method x Time"))
)

cat("PASS: bootstrap 95% CIs were calculated for partial and generalized eta squared.\n")
