.plot_font_family <- "sans"
if (identical(Sys.info()[["sysname"]], "Darwin")) {
    .plot_font_family <- "Hiragino Sans"
} else if (identical(Sys.info()[["sysname"]], "Linux")) {
    fontconfig_file <- file.path(getwd(), "fontconfig.conf")
    if (file.exists(fontconfig_file) && dir.exists("/srv/shiny-server/.fonts")) {
        Sys.setenv(FONTCONFIG_FILE = normalizePath(fontconfig_file))
        .plot_font_family <- "Noto Sans CJK JP"
    }
}

library(shiny)
library(shinyAce)
library(ggplot2)

.anovakun_env <- new.env(parent = globalenv())
source("anovakun_489.R", local = .anovakun_env, encoding = "UTF-8")
source("validation.R", local = TRUE, encoding = "UTF-8")
source("presentation.R", local = TRUE, encoding = "UTF-8")

.anovakun_version <- "4.8.9"



shinyServer(function(input, output, session) {

    analysis_spec <- reactive({
        tryCatch(
            make_anova_spec(
                input$factor,
                input$one.design,
                input$two.design,
                input$factor1.level,
                input$factor2.level
            ),
            error = function(e) structure(
                list(message = conditionMessage(e)),
                class = "anova_spec_error"
            )
        )
    })

    input_check <- debounce(reactive({
        spec <- analysis_spec()
        if (inherits(spec, "anova_spec_error")) {
            return(list(valid = FALSE, errors = spec$message, data = NULL, spec = NULL))
        }
        check <- validate_anova_input(input$text, spec)
        check$spec <- spec
        check
    }), millis = 250)

    validated_state <- reactive({
        check <- input_check()
        validate(need(check$valid, paste(check$errors, collapse = "\n")))
        check
    })

    input_data <- reactive(validated_state()$data)

    validated_spec <- reactive(validated_state()$spec)

    analysis_labels <- reactive({
        derive_analysis_labels(
            input_data(),
            validated_spec(),
            factor1_name = input$factor1.name,
            factor2_name = input$factor2.name,
            outcome_name = input$outcome.name
        )
    })

    output$input_validation_status <- renderUI({
        check <- input_check()
        if (isTRUE(check$valid)) {
            div(
                class = "alert alert-success",
                strong("Input checked: "),
                check$summary
            )
        } else {
            div(
                class = "alert alert-danger",
                strong("Input errors"),
                tags$ul(lapply(check$errors, tags$li))
            )
        }
    })

    observeEvent(list(input$factor, input$one.design, input$two.design), {
        config <- sample_config_for_analysis(
            input$factor,
            input$one.design,
            input$two.design
        )
        req(!is.null(config))
        sample_text <- input[[config$input_id]]
        req(!is.null(sample_text))

        updateNumericInput(
            session,
            "factor1.level",
            value = config$factor1_level
        )
        if (!is.null(config$factor2_level)) {
            updateNumericInput(
                session,
                "factor2.level",
                value = config$factor2_level
            )
        }
        if (!is.null(config$factor1_name)) {
            updateTextInput(session, "factor1.name", value = config$factor1_name)
        }
        if (!is.null(config$factor2_name)) {
            updateTextInput(session, "factor2.name", value = config$factor2_name)
        }
        if (!is.null(config$outcome_name)) {
            updateTextInput(session, "outcome.name", value = config$outcome_name)
        }
        updateAceEditor(session, "text", value = sample_text)
    }, ignoreInit = TRUE, priority = 100)

    run_anova <- function(dat, spec, labels, bootstrap_ci = FALSE, technical = FALSE) {
        args <- c(
            list(dataset = dat, design = spec$design),
            labels$factor_args
        )

        do.call(
            .anovakun_env$anovakun,
            c(args, list(
                mau = TRUE, auto = TRUE, holm = TRUE,
                peta = TRUE, geta = TRUE,
                besci = bootstrap_ci, tech = technical
            ))
        )
    }

    anova_result <- reactive({
        spec <- validated_spec()
        validate(need(
            !inherits(spec, "anova_spec_error"),
            if (inherits(spec, "anova_spec_error")) spec$message else "Invalid analysis settings."
        ))
        run_anova(input_data(), spec, analysis_labels(), technical = TRUE)
    })

    anova_output <- reactive({
        spec <- validated_spec()
        paste(
            capture.output(run_anova(
                input_data(),
                spec,
                analysis_labels()
            )),
            collapse = "\n"
        )
    })

    effect_ci <- reactiveValues(
        result = NULL,
        error = NULL,
        stale = FALSE
    )

    effect_ci_inputs <- reactive({
        list(
            text = input$text,
            factor = input$factor,
            one_design = input$one.design,
            two_design = input$two.design,
            factor1_level = input$factor1.level,
            factor2_level = input$factor2.level,
            factor1_name = input$factor1.name,
            factor2_name = input$factor2.name,
            outcome_name = input$outcome.name
        )
    })

    observeEvent(effect_ci_inputs(), {
        if (!is.null(effect_ci$result) || !is.null(effect_ci$error)) {
            effect_ci$result <- NULL
            effect_ci$error <- NULL
            effect_ci$stale <- TRUE
        }
    }, ignoreInit = TRUE)

    observeEvent(input$effect_size_ci, {
        effect_ci$result <- NULL
        effect_ci$error <- NULL
        effect_ci$stale <- FALSE

        result <- tryCatch({
            text <- isolate(input$text)
            spec <- isolate(analysis_spec())
            if (inherits(spec, "anova_spec_error")) {
                stop(spec$message, call. = FALSE)
            }
            check <- validate_anova_input(text, spec)
            if (!check$valid) {
                stop(paste(check$errors, collapse = " "), call. = FALSE)
            }
            labels <- derive_analysis_labels(
                check$data,
                spec,
                factor1_name = isolate(input$factor1.name),
                factor2_name = isolate(input$factor2.name),
                outcome_name = isolate(input$outcome.name)
            )

            withProgress(
                message = "Calculating 95% confidence intervals",
                detail = "Running 2,000 bootstrap samples...",
                value = 0.1,
                {
                    value <- run_anova(
                        check$data,
                        spec,
                        labels,
                        bootstrap_ci = TRUE,
                        technical = TRUE
                    )
                    incProgress(0.9, detail = "Completed")
                    value
                }
            )
        }, error = function(e) {
            effect_ci$error <- conditionMessage(e)
            NULL
        })

        if (!is.null(result)) {
            effect_ci$result <- result
        }
    }, ignoreInit = TRUE)

    output$effect_size_ci_status <- renderUI({
        if (isTRUE(effect_ci$stale)) {
            p(
                class = "text-warning",
                "The data or analysis settings changed. Click the button to recalculate the confidence intervals."
            )
        } else if (!is.null(effect_ci$error)) {
            p(
                class = "text-danger",
                paste("The confidence intervals could not be calculated:", effect_ci$error)
            )
        }
    })

    format_effect_size_ci <- function(result) {
        capture.output({
            info <- result[["EFFECT SIZE INFORMATION"]]
            cat(info[[1]], sep = "\n")
            cat("\n")
            print(info[[2]], row.names = FALSE)

            post_results <- result[["POST ANALYSES"]]
            for (name in names(post_results)) {
                post_result <- post_results[[name]]
                if (is.list(post_result) && !is.null(post_result$bescitab)) {
                    cat("\nSimple effects for", name, "\n")
                    print(post_result$bescitab, row.names = FALSE)
                }
            }

            invisible(NULL)
        })
    }

    output$effect_size_ci_out <- renderText({
        result <- effect_ci$result
        req(!is.null(result))
        paste(format_effect_size_ci(result), collapse = "\n")
    })

    output$download_results <- downloadHandler(
        filename = function() {
            paste0("langtest-anova-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".txt")
        },
        content = function(file) {
            spec <- analysis_spec()
            if (inherits(spec, "anova_spec_error")) {
                stop(spec$message, call. = FALSE)
            }
            check <- validate_anova_input(input$text, spec)
            if (!check$valid) {
                stop(paste(check$errors, collapse = " "), call. = FALSE)
            }
            labels <- derive_analysis_labels(
                check$data,
                spec,
                factor1_name = input$factor1.name,
                factor2_name = input$factor2.name,
                outcome_name = input$outcome.name
            )
            structured_result <- run_anova(
                check$data,
                spec,
                labels,
                technical = TRUE
            )
            friendly_table <- make_friendly_anova_table(
                structured_result,
                labels
            )

            lines <- c(
                "Langtest ANOVA Results",
                paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
                paste("R version:", R.version.string),
                paste("ANOVA-kun version:", .anovakun_version),
                paste("Design:", spec$design),
                paste("Factors:", paste(labels$factor_names, collapse = " x ")),
                paste("Outcome:", labels$outcome),
                "",
                "Readable ANOVA Table",
                paste(capture.output(print(friendly_table, row.names = FALSE)), collapse = "\n"),
                "",
                "Technical Output from ANOVA-kun",
                paste(
                    capture.output(run_anova(check$data, spec, labels)),
                    collapse = "\n"
                )
            )

            if (!is.null(effect_ci$result)) {
                lines <- c(
                    lines,
                    "",
                    format_effect_size_ci(effect_ci$result)
                )
            } else {
                lines <- c(
                    lines,
                    "",
                    "Effect-size 95% confidence intervals were not calculated."
                )
            }

            writeLines(lines, file, useBytes = TRUE)
        }
    )

    output$analysis_context <- renderUI({
        labels <- analysis_labels()
        details <- c(
            paste("Outcome:", labels$outcome),
            paste("Factors:", paste(labels$factor_names, collapse = " × "))
        )
        notes <- list(
            p(strong(details[1]), br(), details[2])
        )
        if (validated_spec()$design %in% c("sA", "AsB", "sAB")) {
            notes <- c(notes, list(
                p(class = "help-block", "Participant identifies repeated-measures error terms.")
            ))
        }
        if (!isTRUE(labels$within_levels_inferred)) {
            notes <- c(notes, list(
                p(
                    class = "text-warning",
                    "Factor levels could not be inferred from the column headers; generic level labels are used."
                )
            ))
        }
        do.call(div, c(list(class = "alert alert-info"), notes))
    })

    output$anova_table <- renderTable({
        make_friendly_anova_table(anova_result(), analysis_labels())
    }, striped = TRUE, bordered = TRUE, spacing = "s", align = "lrrrrrrr")

    plot_long_data <- reactive({
        to_plot_long(input_data(), validated_spec(), analysis_labels())
    })

    build_anova_plot <- function(x_factor = 1L) {
        spec <- validated_spec()
        labels <- analysis_labels()
        long <- plot_long_data()

        if (x_factor == 2L && !".factor2" %in% names(long)) {
            return(NULL)
        }

        summary <- make_plot_summary(long, spec, x_factor)
        within_factors <- switch(
            spec$design,
            As = integer(0),
            sA = 1L,
            ABs = integer(0),
            AsB = 2L,
            sAB = 1:2
        )
        x_is_within <- x_factor %in% within_factors

        if (x_factor == 1L) {
            long$.x <- long$.factor1
            summary$.x <- summary$.factor1
            x_name <- labels$factor_names[1]
            if (".factor2" %in% names(long)) {
                long$.group <- long$.factor2
                summary$.group <- summary$.factor2
                group_name <- labels$factor_names[2]
            }
        } else {
            long$.x <- long$.factor2
            summary$.x <- summary$.factor2
            x_name <- labels$factor_names[2]
            long$.group <- long$.factor1
            summary$.group <- summary$.factor1
            group_name <- labels$factor_names[1]
        }

        plot <- ggplot()
        individual_note <- NULL

        if (isTRUE(input$show_individual)) {
            if (x_is_within) {
                if (".group" %in% names(long)) {
                    plot <- plot +
                        geom_line(
                            data = long,
                            aes(
                                x = .x,
                                y = .value,
                                colour = .group,
                                group = interaction(.subject, .group)
                            ),
                            alpha = 0.13,
                            linewidth = 0.35
                        ) +
                        geom_point(
                            data = long,
                            aes(x = .x, y = .value, colour = .group),
                            alpha = 0.13,
                            size = 0.9
                        )
                } else {
                    plot <- plot +
                        geom_line(
                            data = long,
                            aes(x = .x, y = .value, group = .subject),
                            colour = "grey55",
                            alpha = 0.22,
                            linewidth = 0.35
                        ) +
                        geom_point(
                            data = long,
                            aes(x = .x, y = .value),
                            colour = "grey55",
                            alpha = 0.22,
                            size = 0.9
                        )
                }
                individual_note <- "Thin lines show individual participants."
            } else {
                if (".group" %in% names(long)) {
                    plot <- plot + geom_point(
                        data = long,
                        aes(x = .x, y = .value, colour = .group),
                        position = position_jitterdodge(
                            jitter.width = 0.08,
                            dodge.width = 0.45
                        ),
                        alpha = 0.22,
                        size = 1.2
                    )
                } else {
                    plot <- plot + geom_jitter(
                        data = long,
                        aes(x = .x, y = .value),
                        width = 0.08,
                        height = 0,
                        colour = "grey50",
                        alpha = 0.28,
                        size = 1.2
                    )
                }
                individual_note <- "Faint points show individual observations."
            }
        }

        if (".group" %in% names(summary)) {
            dodge <- position_dodge(width = 0.45)
            plot <- plot +
                geom_errorbar(
                    data = summary,
                    aes(
                        x = .x,
                        ymin = CI_L,
                        ymax = CI_U,
                        colour = .group,
                        group = .group
                    ),
                    width = 0.10,
                    linewidth = 0.65,
                    position = dodge
                ) +
                geom_line(
                    data = summary,
                    aes(
                        x = .x,
                        y = Mean,
                        colour = .group,
                        linetype = .group,
                        group = .group
                    ),
                    linewidth = 0.8,
                    position = dodge
                ) +
                geom_point(
                    data = summary,
                    aes(
                        x = .x,
                        y = Mean,
                        colour = .group,
                        shape = .group,
                        group = .group
                    ),
                    size = 3,
                    position = dodge
                ) +
                scale_colour_discrete(name = group_name) +
                scale_linetype_discrete(name = group_name) +
                scale_shape_discrete(name = group_name)
        } else {
            plot <- plot +
                geom_errorbar(
                    data = summary,
                    aes(x = .x, ymin = CI_L, ymax = CI_U),
                    width = 0.10,
                    linewidth = 0.65,
                    colour = "#0072B2"
                ) +
                geom_point(
                    data = summary,
                    aes(x = .x, y = Mean),
                    size = 3,
                    colour = "#0072B2"
                )
            if (x_is_within) {
                plot <- plot + geom_line(
                    data = summary,
                    aes(x = .x, y = Mean, group = 1),
                    linewidth = 0.8,
                    colour = "#0072B2"
                )
            }
        }

        caption <- paste(
            c(
                individual_note,
                paste0("Error bars show ", unique(summary$ci_method), ".")
            ),
            collapse = " "
        )
        plot <- plot +
            labs(
                title = paste(labels$outcome, "by", x_name),
                x = x_name,
                y = labels$outcome,
                caption = caption
            ) +
            theme_minimal(base_size = 13, base_family = .plot_font_family) +
            theme(
                panel.grid.minor = element_blank(),
                legend.position = "top",
                plot.caption = element_text(hjust = 0)
            )

        if (identical(input$axis, "min.max")) {
            plot <- plot + coord_cartesian(ylim = range(long$.value, na.rm = TRUE))
        } else if (identical(input$axis, "dfn") &&
                is.finite(input$dfn.min) && is.finite(input$dfn.max) &&
                input$dfn.min < input$dfn.max) {
            plot <- plot + coord_cartesian(ylim = c(input$dfn.min, input$dfn.max))
        }

        plot
    }


    output$AnovaPlot1 <- renderPlot({
        print(build_anova_plot(1L))
    })
    
    
    
    

    output$AnovaPlot2 <- renderPlot({
        plot <- build_anova_plot(2L)
        req(!is.null(plot))
        print(plot)
    })




    output$anovakun.out <- renderText({
        anova_output()
    })

})
