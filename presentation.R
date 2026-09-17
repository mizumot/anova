clean_display_label <- function(value, fallback) {
    if (is.null(value) || length(value) == 0 || is.na(value[1]) || !nzchar(trimws(value[1]))) {
        fallback
    } else {
        trimws(value[1])
    }
}

parse_two_within_levels <- function(headers, factor1_level, factor2_level) {
    tokens <- strsplit(headers, "[._:|]+")
    can_parse <- all(lengths(tokens) >= 2)

    if (can_parse) {
        factor1 <- vapply(tokens, `[[`, character(1), 1)
        factor2 <- vapply(tokens, function(x) x[[length(x)]], character(1))
        factor1_names <- unique(factor1)
        factor2_names <- unique(factor2)
        combinations <- paste(factor1, factor2, sep = "\r")

        can_parse <- length(factor1_names) == factor1_level &&
            length(factor2_names) == factor2_level &&
            length(unique(combinations)) == factor1_level * factor2_level
    }

    if (can_parse) {
        list(
            factor1 = factor1_names,
            factor2 = factor2_names,
            inferred = TRUE
        )
    } else {
        list(
            factor1 = paste("Level", seq_len(factor1_level)),
            factor2 = paste("Level", seq_len(factor2_level)),
            inferred = FALSE
        )
    }
}

derive_analysis_labels <- function(dat, spec, factor1_name = NULL,
                                   factor2_name = NULL, outcome_name = NULL) {
    parsed_within <- NULL

    if (spec$design == "As") {
        factor_names <- names(dat)[1]
        factor_levels <- list(unique(as.character(dat[[1]])))
        outcome <- names(dat)[2]
    } else if (spec$design == "sA") {
        factor_names <- clean_display_label(factor1_name, "Condition")
        factor_levels <- list(names(dat))
        outcome <- clean_display_label(outcome_name, "Score")
    } else if (spec$design == "ABs") {
        factor_names <- names(dat)[1:2]
        factor_levels <- list(
            unique(as.character(dat[[1]])),
            unique(as.character(dat[[2]]))
        )
        outcome <- names(dat)[3]
    } else if (spec$design == "AsB") {
        factor_names <- c(
            names(dat)[1],
            clean_display_label(factor2_name, "Condition")
        )
        factor_levels <- list(
            unique(as.character(dat[[1]])),
            names(dat)[-1]
        )
        outcome <- clean_display_label(outcome_name, "Score")
    } else {
        parsed_within <- parse_two_within_levels(
            names(dat),
            spec$factor1_level,
            spec$factor2_level
        )
        factor_names <- c(
            clean_display_label(factor1_name, "Factor 1"),
            clean_display_label(factor2_name, "Factor 2")
        )
        factor_levels <- list(
            parsed_within$factor1,
            parsed_within$factor2
        )
        outcome <- clean_display_label(outcome_name, "Score")
    }

    factor_names <- trimws(factor_names)
    factor_names[!nzchar(factor_names) | factor_names == "s"] <- paste0(
        "Factor ", which(!nzchar(factor_names) | factor_names == "s")
    )
    factor_names <- make.unique(factor_names, sep = " ")

    subject_position <- match("s", strsplit(spec$design, "", fixed = TRUE)[[1]])
    between_count <- subject_position - 1L

    list(
        factor_names = factor_names,
        factor_levels = factor_levels,
        factor_args = setNames(factor_levels, factor_names),
        outcome = clean_display_label(outcome, "Score"),
        subject = "Participant",
        between_count = between_count,
        within_levels_inferred = if (is.null(parsed_within)) TRUE else parsed_within$inferred
    )
}

to_plot_long <- function(dat, spec, labels) {
    subjects <- paste0("P", seq_len(nrow(dat)))

    if (spec$design == "As") {
        long <- data.frame(
            .subject = subjects,
            .factor1 = as.character(dat[[1]]),
            .value = dat[[2]],
            stringsAsFactors = FALSE
        )
    } else if (spec$design == "sA") {
        long <- data.frame(
            .subject = rep(subjects, times = spec$factor1_level),
            .factor1 = rep(labels$factor_levels[[1]], each = nrow(dat)),
            .value = unlist(dat, use.names = FALSE),
            stringsAsFactors = FALSE
        )
    } else if (spec$design == "ABs") {
        long <- data.frame(
            .subject = subjects,
            .factor1 = as.character(dat[[1]]),
            .factor2 = as.character(dat[[2]]),
            .value = dat[[3]],
            stringsAsFactors = FALSE
        )
    } else if (spec$design == "AsB") {
        long <- data.frame(
            .subject = rep(subjects, times = spec$factor2_level),
            .factor1 = rep(as.character(dat[[1]]), times = spec$factor2_level),
            .factor2 = rep(labels$factor_levels[[2]], each = nrow(dat)),
            .value = unlist(dat[-1], use.names = FALSE),
            stringsAsFactors = FALSE
        )
    } else {
        column_factor1 <- rep(
            labels$factor_levels[[1]],
            each = spec$factor2_level
        )
        column_factor2 <- rep(
            labels$factor_levels[[2]],
            times = spec$factor1_level
        )
        long <- data.frame(
            .subject = rep(subjects, times = ncol(dat)),
            .factor1 = rep(column_factor1, each = nrow(dat)),
            .factor2 = rep(column_factor2, each = nrow(dat)),
            .value = unlist(dat, use.names = FALSE),
            stringsAsFactors = FALSE
        )
    }

    long$.subject <- factor(long$.subject, levels = subjects)
    long$.factor1 <- factor(
        long$.factor1,
        levels = labels$factor_levels[[1]]
    )
    if (length(labels$factor_levels) == 2) {
        long$.factor2 <- factor(
            long$.factor2,
            levels = labels$factor_levels[[2]]
        )
    }
    long
}

summarize_plot_groups <- function(dat, group_columns, value_column) {
    keys <- do.call(
        interaction,
        c(dat[group_columns], list(drop = TRUE, lex.order = TRUE))
    )
    groups <- split(seq_len(nrow(dat)), keys)

    rows <- lapply(groups, function(index) {
        values <- dat[[value_column]][index]
        row <- as.list(dat[index[1], group_columns, drop = FALSE])
        row$n <- length(values)
        row$Mean <- mean(values)
        row$SD <- if (length(values) > 1) stats::sd(values) else NA_real_
        row
    })
    result <- do.call(rbind.data.frame, c(rows, stringsAsFactors = FALSE))
    rownames(result) <- NULL
    result
}

make_plot_summary <- function(long, spec, x_factor) {
    group_columns <- if (".factor2" %in% names(long)) {
        c(".factor1", ".factor2")
    } else {
        ".factor1"
    }
    raw_summary <- summarize_plot_groups(long, group_columns, ".value")

    within_factors <- switch(
        spec$design,
        As = integer(0),
        sA = 1L,
        ABs = integer(0),
        AsB = 2L,
        sAB = 1:2
    )
    use_within_ci <- x_factor %in% within_factors

    if (use_within_ci) {
        subject_mean <- ave(long$.value, long$.subject, FUN = mean)
        if (spec$design == "AsB") {
            group_grand_mean <- ave(long$.value, long$.factor1, FUN = mean)
            repeated_conditions <- spec$factor2_level
        } else {
            group_grand_mean <- rep(mean(long$.value), nrow(long))
            repeated_conditions <- if (spec$design == "sA") {
                spec$factor1_level
            } else {
                spec$factor1_level * spec$factor2_level
            }
        }
        long$.normalized <- long$.value - subject_mean + group_grand_mean
        normalized_summary <- summarize_plot_groups(
            long,
            group_columns,
            ".normalized"
        )
        correction <- sqrt(repeated_conditions / (repeated_conditions - 1))
        half_width <- stats::qt(0.975, normalized_summary$n - 1) *
            normalized_summary$SD / sqrt(normalized_summary$n) * correction
        ci_method <- "within-subject 95% CI"
    } else {
        half_width <- stats::qt(0.975, raw_summary$n - 1) *
            raw_summary$SD / sqrt(raw_summary$n)
        ci_method <- "95% CI"
    }

    raw_summary$CI_L <- raw_summary$Mean - half_width
    raw_summary$CI_U <- raw_summary$Mean + half_width
    raw_summary$ci_method <- ci_method

    raw_summary$.factor1 <- factor(
        raw_summary$.factor1,
        levels = levels(long$.factor1)
    )
    if (".factor2" %in% names(raw_summary)) {
        raw_summary$.factor2 <- factor(
            raw_summary$.factor2,
            levels = levels(long$.factor2)
        )
        raw_summary <- raw_summary[order(raw_summary$.factor1, raw_summary$.factor2), ]
    } else {
        raw_summary <- raw_summary[order(raw_summary$.factor1), ]
    }
    rownames(raw_summary) <- NULL
    raw_summary
}

friendly_source_label <- function(source, labels) {
    if (is.na(source) || source %in% c("Error", "Total")) {
        return(source)
    }
    parts <- strsplit(source, " x ", fixed = TRUE)[[1]]
    if (!"s" %in% parts) {
        return(paste(parts, collapse = " × "))
    }
    if (length(parts) == 1) {
        return(labels$subject)
    }

    terms <- setdiff(parts, "s")
    between_names <- if (labels$between_count > 0) {
        labels$factor_names[seq_len(labels$between_count)]
    } else {
        character()
    }
    between_terms <- terms[terms %in% between_names]
    within_terms <- terms[!terms %in% between_names]

    participant <- if (length(between_terms) > 0) {
        paste(
            labels$subject,
            paste(between_terms, collapse = " × "),
            sep = " within "
        )
    } else {
        labels$subject
    }
    if (length(within_terms) > 0) {
        if (length(between_terms) > 0) {
            participant <- paste0("(", participant, ")")
        }
        participant <- paste(
            participant,
            paste(within_terms, collapse = " × "),
            sep = " × "
        )
    }
    paste0(participant, " (error)")
}

format_table_number <- function(value, digits = 3) {
    ifelse(
        is.na(value),
        "",
        sub("\\.?0+$", "", formatC(value, format = "f", digits = digits))
    )
}

format_table_p <- function(value) {
    vapply(value, function(x) {
        if (is.na(x)) {
            ""
        } else if (x < 0.001) {
            "< .001"
        } else {
            sub("^0", "", sprintf("%.3f", x))
        }
    }, character(1))
}

make_friendly_anova_table <- function(result, labels) {
    table <- result[["ANOVA TABLE"]][[2]]
    data.frame(
        Source = vapply(
            as.character(table$source.col),
            friendly_source_label,
            character(1),
            labels = labels
        ),
        SS = format_table_number(table$ss.col),
        df = format_table_number(table$df.col),
        MS = format_table_number(table$ms.col),
        F = format_table_number(table$f.col),
        p = format_table_p(table$p.col),
        `Partial η²` = format_table_number(table[["p.eta^2"]]),
        `Generalized η²` = format_table_number(table[["G.eta^2"]]),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}
