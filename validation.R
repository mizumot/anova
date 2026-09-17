make_anova_spec <- function(factor_type, one_design, two_design,
                            factor1_level, factor2_level = NA) {
    check_level <- function(value, label) {
        if (length(value) != 1 || is.na(value) || value < 2 || value != floor(value)) {
            stop(label, " must be an integer of 2 or greater.", call. = FALSE)
        }
        as.integer(value)
    }

    factor1_level <- check_level(factor1_level, "The number of Factor 1 levels")

    if (identical(factor_type, "oneway")) {
        design <- switch(
            one_design,
            Between = "As",
            Within = "sA",
            stop("Select a valid one-way design.", call. = FALSE)
        )
        factor2_level <- NA_integer_
    } else if (identical(factor_type, "twoway")) {
        factor2_level <- check_level(factor2_level, "The number of Factor 2 levels")
        design <- switch(
            two_design,
            Factor1Between_Factor2Between = "ABs",
            Factor1Between_Factor2Within = "AsB",
            Factor1Within_Factor2Within = "sAB",
            stop("Select a valid two-way design.", call. = FALSE)
        )
    } else {
        stop("Select one-way or two-way ANOVA.", call. = FALSE)
    }

    expected_columns <- switch(
        design,
        As = 2L,
        sA = factor1_level,
        ABs = 3L,
        AsB = 1L + factor2_level,
        sAB = factor1_level * factor2_level
    )
    numeric_columns <- switch(
        design,
        As = 2L,
        sA = seq_len(expected_columns),
        ABs = 3L,
        AsB = 2L:expected_columns,
        sAB = seq_len(expected_columns)
    )
    factor_columns <- switch(
        design,
        As = 1L,
        sA = integer(0),
        ABs = 1:2,
        AsB = 1L,
        sAB = integer(0)
    )

    list(
        design = design,
        levels = if (design %in% c("As", "sA")) {
            list(factor1_level)
        } else {
            list(factor1_level, factor2_level)
        },
        factor1_level = factor1_level,
        factor2_level = factor2_level,
        expected_columns = expected_columns,
        numeric_columns = numeric_columns,
        factor_columns = factor_columns
    )
}

sample_config_for_analysis <- function(factor_type, one_design, two_design) {
    if (identical(factor_type, "oneway")) {
        switch(
            one_design,
            Between = list(
                input_id = "text1", factor1_level = 3L, factor2_level = NULL,
                factor1_name = "Class", factor2_name = NULL, outcome_name = "Score"
            ),
            Within = list(
                input_id = "text2", factor1_level = 3L, factor2_level = NULL,
                factor1_name = "Condition", factor2_name = NULL, outcome_name = "Score"
            ),
            NULL
        )
    } else if (identical(factor_type, "twoway")) {
        switch(
            two_design,
            Factor1Between_Factor2Between = list(
                input_id = "text3", factor1_level = 2L, factor2_level = 3L,
                factor1_name = "Classroom", factor2_name = "Class Size", outcome_name = "Score"
            ),
            Factor1Between_Factor2Within = list(
                input_id = "text4", factor1_level = 2L, factor2_level = 3L,
                factor1_name = "Method", factor2_name = "Time", outcome_name = "Score"
            ),
            Factor1Within_Factor2Within = list(
                input_id = "text5", factor1_level = 2L, factor2_level = 3L,
                factor1_name = "Type", factor2_name = "Time", outcome_name = "Score"
            ),
            NULL
        )
    } else {
        NULL
    }
}

format_cell_locations <- function(rows, column_name, values = NULL, limit = 5L) {
    shown <- seq_len(min(length(rows), limit))
    locations <- paste0("row ", rows[shown], ", column '", column_name, "'")
    if (!is.null(values)) {
        locations <- paste0(locations, " ('", values[shown], "')")
    }
    if (length(rows) > limit) {
        locations <- c(locations, paste0("and ", length(rows) - limit, " more"))
    }
    paste(locations, collapse = "; ")
}

validate_anova_input <- function(text, spec) {
    errors <- character()

    if (is.null(text) || !nzchar(trimws(text))) {
        return(list(valid = FALSE, errors = "No data were entered.", data = NULL))
    }

    dat <- tryCatch(
        read.csv(
            text = text,
            sep = "\t",
            header = TRUE,
            stringsAsFactors = FALSE,
            check.names = FALSE,
            na.strings = c("", "NA", "NaN")
        ),
        error = function(e) e
    )
    if (inherits(dat, "error")) {
        return(list(
            valid = FALSE,
            errors = paste("The tab-separated data could not be read:", conditionMessage(dat)),
            data = NULL
        ))
    }

    if (nrow(dat) == 0) {
        errors <- c(errors, "The data must contain at least one data row below the header.")
    }
    if (ncol(dat) != spec$expected_columns) {
        errors <- c(
            errors,
            paste0(
                "Column count: the ", spec$design, " design requires ",
                spec$expected_columns, " columns, but ", ncol(dat), " were found."
            )
        )
    }
    if (any(!nzchar(trimws(names(dat))))) {
        errors <- c(errors, "Every column must have a header in the first row.")
    }
    if (anyDuplicated(names(dat))) {
        duplicates <- unique(names(dat)[duplicated(names(dat))])
        errors <- c(
            errors,
            paste("Column headers must be unique. Duplicated:", paste(duplicates, collapse = ", "))
        )
    }

    missing_cells <- which(is.na(dat), arr.ind = TRUE)
    if (nrow(missing_cells) > 0) {
        shown <- seq_len(min(nrow(missing_cells), 5L))
        locations <- mapply(
            function(row, column) paste0("row ", row, ", column '", names(dat)[column], "'"),
            missing_cells[shown, "row"], missing_cells[shown, "col"]
        )
        if (nrow(missing_cells) > 5L) {
            locations <- c(locations, paste0("and ", nrow(missing_cells) - 5L, " more"))
        }
        errors <- c(
            errors,
            paste("Missing values were found at", paste(locations, collapse = "; "))
        )
    }

    numeric_columns <- intersect(spec$numeric_columns, seq_len(ncol(dat)))
    for (column in numeric_columns) {
        original <- dat[[column]]
        converted <- suppressWarnings(as.numeric(trimws(as.character(original))))
        invalid <- which(!is.na(original) & is.na(converted))
        if (length(invalid) > 0) {
            errors <- c(
                errors,
                paste0(
                    "Numeric data are required at ",
                    format_cell_locations(
                        invalid,
                        names(dat)[column],
                        as.character(original[invalid])
                    ),
                    "."
                )
            )
        } else {
            non_finite <- which(!is.na(converted) & !is.finite(converted))
            if (length(non_finite) > 0) {
                errors <- c(
                    errors,
                    paste0(
                        "Finite numeric data are required at ",
                        format_cell_locations(
                            non_finite,
                            names(dat)[column],
                            as.character(original[non_finite])
                        ),
                        "."
                    )
                )
            }
            dat[[column]] <- converted
        }
    }

    factor_values <- lapply(
        intersect(spec$factor_columns, seq_len(ncol(dat))),
        function(column) trimws(as.character(dat[[column]]))
    )

    if (spec$design == "As" && length(factor_values) == 1) {
        found <- length(unique(factor_values[[1]][!is.na(factor_values[[1]])]))
        if (found != spec$factor1_level) {
            errors <- c(errors, paste0(
                "Factor 1 condition count: expected ", spec$factor1_level,
                ", but found ", found, "."
            ))
        }
    } else if (spec$design == "sA") {
        found <- ncol(dat)
        if (found != spec$factor1_level) {
            errors <- c(errors, paste0(
                "Factor 1 condition count: expected ", spec$factor1_level,
                " score columns, but found ", found, "."
            ))
        }
    } else if (spec$design == "ABs" && length(factor_values) == 2) {
        found1 <- length(unique(factor_values[[1]][!is.na(factor_values[[1]])]))
        found2 <- length(unique(factor_values[[2]][!is.na(factor_values[[2]])]))
        if (found1 != spec$factor1_level) {
            errors <- c(errors, paste0(
                "Factor 1 condition count: expected ", spec$factor1_level,
                ", but found ", found1, "."
            ))
        }
        if (found2 != spec$factor2_level) {
            errors <- c(errors, paste0(
                "Factor 2 condition count: expected ", spec$factor2_level,
                ", but found ", found2, "."
            ))
        }
        if (found1 == spec$factor1_level && found2 == spec$factor2_level &&
                !anyNA(factor_values[[1]]) && !anyNA(factor_values[[2]])) {
            factor1_names <- unique(factor_values[[1]])
            factor2_names <- unique(factor_values[[2]])
            expected <- as.vector(outer(factor1_names, factor2_names, paste, sep = " / "))
            observed <- unique(paste(factor_values[[1]], factor_values[[2]], sep = " / "))
            missing_conditions <- setdiff(expected, observed)
            if (length(missing_conditions) > 0) {
                errors <- c(errors, paste(
                    "No observations were found for these Factor 1 / Factor 2 conditions:",
                    paste(missing_conditions, collapse = ", ")
                ))
            }
        }
    } else if (spec$design == "AsB" && length(factor_values) == 1) {
        found1 <- length(unique(factor_values[[1]][!is.na(factor_values[[1]])]))
        found2 <- max(0L, ncol(dat) - 1L)
        if (found1 != spec$factor1_level) {
            errors <- c(errors, paste0(
                "Factor 1 condition count: expected ", spec$factor1_level,
                ", but found ", found1, "."
            ))
        }
        if (found2 != spec$factor2_level) {
            errors <- c(errors, paste0(
                "Factor 2 condition count: expected ", spec$factor2_level,
                " score columns, but found ", found2, "."
            ))
        }
    } else if (spec$design == "sAB") {
        found <- ncol(dat)
        expected <- spec$factor1_level * spec$factor2_level
        if (found != expected) {
            errors <- c(errors, paste0(
                "Within-subject condition count: expected ", expected,
                " score columns (", spec$factor1_level, " x ", spec$factor2_level,
                "), but found ", found, "."
            ))
        }
    }

    list(
        valid = length(errors) == 0,
        errors = unique(errors),
        data = dat,
        summary = paste0(nrow(dat), " rows, ", ncol(dat), " columns; design ", spec$design)
    )
}
