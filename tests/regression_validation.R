args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) {
    stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg))
app_dir <- normalizePath(file.path(dirname(script_path), ".."))

source(file.path(app_dir, "validation.R"), encoding = "UTF-8")

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

specs <- list(
    text1 = make_anova_spec("oneway", "Between", NULL, 3),
    text2 = make_anova_spec("oneway", "Within", NULL, 3),
    text3 = make_anova_spec("twoway", NULL, "Factor1Between_Factor2Between", 2, 3),
    text4 = make_anova_spec("twoway", NULL, "Factor1Between_Factor2Within", 2, 3),
    text5 = make_anova_spec("twoway", NULL, "Factor1Within_Factor2Within", 2, 3)
)

valid_results <- Map(
    function(id, spec) validate_anova_input(example_values[[id]], spec),
    names(specs),
    specs
)
stopifnot(all(vapply(valid_results, function(x) x$valid, logical(1))))

wrong_columns <- validate_anova_input(
    "Group\tScore\tExtra\n1\t10\t20\n2\t11\t21",
    specs$text1
)
missing_value <- validate_anova_input(
    "First\tSecond\tThird\n10\t\t12\n11\t13\t14",
    specs$text2
)
nonnumeric_value <- validate_anova_input(
    "First\tSecond\tThird\n10\toops\t12\n11\t13\t14",
    specs$text2
)
nonfinite_value <- validate_anova_input(
    "First\tSecond\tThird\n10\tInf\t12\n11\t13\t14",
    specs$text2
)
wrong_conditions <- validate_anova_input(
    "Group\tScore\n1\t10\n1\t11\n2\t12",
    specs$text1
)

stopifnot(
    !wrong_columns$valid,
    any(grepl("Column count", wrong_columns$errors, fixed = TRUE)),
    !missing_value$valid,
    any(grepl("Missing values", missing_value$errors, fixed = TRUE)),
    !nonnumeric_value$valid,
    any(grepl("Numeric data", nonnumeric_value$errors, fixed = TRUE)),
    !nonfinite_value$valid,
    any(grepl("Finite numeric data", nonfinite_value$errors, fixed = TRUE)),
    !wrong_conditions$valid,
    any(grepl("condition count", wrong_conditions$errors, fixed = TRUE))
)

sample_configs <- list(
    sample_config_for_analysis("oneway", "Between", NULL),
    sample_config_for_analysis("oneway", "Within", NULL),
    sample_config_for_analysis("twoway", NULL, "Factor1Between_Factor2Between"),
    sample_config_for_analysis("twoway", NULL, "Factor1Between_Factor2Within"),
    sample_config_for_analysis("twoway", NULL, "Factor1Within_Factor2Within")
)
stopifnot(identical(
    vapply(sample_configs, `[[`, character(1), "input_id"),
    paste0("text", 1:5)
))

cat("PASS: input validation and automatic sample mappings work for all five designs.\n")
