args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) != 1) {
    stop("Run this file with Rscript.")
}

script_path <- normalizePath(sub("^--file=", "", file_arg))
app_dir <- normalizePath(file.path(dirname(script_path), ".."))

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

specs <- list(
    text1 = list(design = "As", levels = list(3)),
    text2 = list(design = "sA", levels = list(3)),
    text3 = list(design = "ABs", levels = list(2, 3)),
    text4 = list(design = "AsB", levels = list(2, 3)),
    text5 = list(design = "sAB", levels = list(2, 3))
)

results <- lapply(names(specs), function(id) {
    dat <- read.csv(text = example_values[[id]], sep = "\t")
    call_args <- c(
        list(dataset = dat, design = specs[[id]]$design),
        specs[[id]]$levels,
        list(
            mau = TRUE, auto = TRUE, holm = TRUE,
            peta = TRUE, geta = TRUE, tech = TRUE
        )
    )
    result <- do.call(anovakun, call_args)
    table <- result[["ANOVA TABLE"]][[2]]
    stopifnot(nrow(table) > 0, any(is.finite(table$f.col)))
    result
})
names(results) <- names(specs)

within_table <- results$text5[["ANOVA TABLE"]][[2]]
expected_sources <- c("s", "A", "s x A", "B", "s x B", "A x B", "s x A x B", "Total")
expected_ss <- c(
    4045.68333333333,
    2220.41666666667,
    973.416666666667,
    93143.6333333333,
    2387.36666666667,
    1160.23333333333,
    1213.43333333333,
    105144.183333333
)
expected_f <- c(20.5294923379848, 351.136971000125, 8.60541713595033)

stopifnot(
    identical(as.character(within_table$source.col), expected_sources),
    isTRUE(all.equal(within_table$ss.col, expected_ss, tolerance = 1e-9)),
    isTRUE(all.equal(within_table$f.col[c(2, 4, 6)], expected_f, tolerance = 1e-9)),
    all(c("p.eta^2", "G.eta^2") %in% names(within_table)),
    grepl("anovakun 4.8.9", results$text5$INFORMATION[2], fixed = TRUE)
)

cat("PASS: all five bundled examples completed; sAB error terms match the verified reference.\n")
