args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) {
    stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg))
app_dir <- normalizePath(file.path(dirname(script_path), ".."))

source(file.path(app_dir, "validation.R"), encoding = "UTF-8")
source(file.path(app_dir, "presentation.R"), encoding = "UTF-8")
source(file.path(app_dir, "anovakun_489.R"), encoding = "UTF-8")

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

cases <- list(
    text1 = list(
        spec = make_anova_spec("oneway", "Between", NULL, 3),
        names = c("Class", NA, "Score")
    ),
    text2 = list(
        spec = make_anova_spec("oneway", "Within", NULL, 3),
        names = c("Condition", NA, "Score")
    ),
    text3 = list(
        spec = make_anova_spec("twoway", NULL, "Factor1Between_Factor2Between", 2, 3),
        names = c("Classroom", "Class Size", "Score")
    ),
    text4 = list(
        spec = make_anova_spec("twoway", NULL, "Factor1Between_Factor2Within", 2, 3),
        names = c("Method", "Time", "Score")
    ),
    text5 = list(
        spec = make_anova_spec("twoway", NULL, "Factor1Within_Factor2Within", 2, 3),
        names = c("Type", "Time", "Score")
    )
)

results <- lapply(names(cases), function(id) {
    case <- cases[[id]]
    checked <- validate_anova_input(example_values[[id]], case$spec)
    stopifnot(checked$valid)

    labels <- derive_analysis_labels(
        checked$data,
        case$spec,
        factor1_name = case$names[1],
        factor2_name = case$names[2],
        outcome_name = case$names[3]
    )
    long <- to_plot_long(checked$data, case$spec, labels)
    summary1 <- make_plot_summary(long, case$spec, 1L)
    stopifnot(nrow(long) > 0, all(is.finite(summary1$Mean)))

    if (length(labels$factor_names) == 2) {
        summary2 <- make_plot_summary(long, case$spec, 2L)
        stopifnot(nrow(summary2) > 0, all(is.finite(summary2$Mean)))
    }

    call_args <- c(
        list(dataset = checked$data, design = case$spec$design),
        labels$factor_args,
        list(peta = TRUE, geta = TRUE, tech = TRUE)
    )
    result <- do.call(anovakun, call_args)
    friendly <- make_friendly_anova_table(result, labels)
    stopifnot(
        nrow(friendly) > 0,
        all(c("Source", "Partial η²", "Generalized η²") %in% names(friendly))
    )

    list(labels = labels, long = long, friendly = friendly)
})
names(results) <- names(cases)

stopifnot(
    identical(results$text4$labels$factor_names, c("Method", "Time")),
    any(grepl("Participant within Method", results$text4$friendly$Source, fixed = TRUE)),
    identical(results$text5$labels$factor_levels[[1]], c("TypeA", "TypeB")),
    identical(results$text5$labels$factor_levels[[2]], c("Pre", "Post", "Delayed")),
    any(grepl("Participant × Type", results$text5$friendly$Source, fixed = TRUE))
)

cat("PASS: readable labels, ANOVA tables, and plot data work for all five designs.\n")
