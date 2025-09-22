library(shiny)
library(bslib)
library(promises)
library(fastmap)
library(duckdb)
library(DBI)
library(fontawesome)
library(reactable)
library(plotly)
library(ggplot2)
library(ggridges)
library(dplyr)
library(here)
library(ellmer)
library(shinychat)

# Open the duckdb database
conn <- dbConnect(duckdb(), dbdir = here("IC50_expr_cnv.duckdb"), read_only = TRUE)
# Close the database when the app stops
onStop(\() dbDisconnect(conn))

# gpt-4o does much better than gpt-4o-mini, especially at interpreting plots
gemini_model <- "gemini-2.5-flash"

# Dynamically create the system prompt, based on the real data. For an actually
# large database, you wouldn't want to retrieve all the data like this, but
# instead either hand-write the schema or write your own routine that is more
# efficient than system_prompt().
system_prompt_str <- system_prompt(dbGetQuery(conn, "SELECT * FROM IC50_expr_cnv"), "IC50_expr_cnv")

# This is the greeting that should initially appear in the sidebar when the app
# loads.
greeting <- paste(readLines(here("greeting.md")), collapse = "\n")

icon_explain <- tags$img(src = "stars.svg")

ui <- page_sidebar(
  style = "background-color: rgb(248, 248, 248);",
  title = "Compound IC50 and Gene Expression/CNV Explorer",
  includeCSS(here("styles.css")),
  sidebar = sidebar(
    width = 400,
    style = "height: 100%;",
    chat_ui("chat", height = "100%", fill = TRUE)
  ),
  useBusyIndicators(),

  # 🏷️ Header
  textOutput("show_title", container = h3),
  verbatimTextOutput("show_query") |>
    tagAppendAttributes(style = "max-height: 100px; overflow: auto;"),

  # 🎯 Value boxes
  layout_columns(
    fill = FALSE,
    value_box(
      showcase = fa_i("vial-virus"),
      "Total cell lines",
      textOutput("total_cell_lines", inline = TRUE)
    ),
    value_box(
      showcase = fa_i("florin-sign"),
      "Median IC50",
      textOutput("median_IC50", inline = TRUE)
    ),
    value_box(
      showcase = fa_i("arrow-trend-up"),
      "Median AUC",
      textOutput("median_AUC", inline = TRUE)
    ),
  ),
  layout_columns(
    style = "min-height: 800px;",
    col_widths = c(6, 6, 12),

    # 📊 Scatter plot
    card(
      card_header(
        class = "d-flex justify-content-between align-items-center",
        "IC50 vs AUC",
        span(
          actionLink(
            "interpret_scatter",
            icon_explain,
            class = "me-3 text-decoration-none",
            aria_label = "Explain scatter plot"
          ),
          popover(
            title = "Add a color variable", placement = "top",
            fa_i("ellipsis"),
            radioButtons(
              "scatter_color",
              NULL,
              c("none", "SampleCollectionSite", "Mutations.NRAS", "Mutations.BRAF", "Mutations.PIK3CA", 
              "EXPR.NRAS", "EXPR.BRAF", "EXPR.PIK3CA", "CNV.NRAS", "CNV.BRAF", "CNV.PIK3CA"),
              inline = TRUE
            )
          )
        )
      ),
      plotlyOutput("scatterplot")
    ),
  
  # 📊 Box plot for sample collection site
  card(
      card_header(
        class = "d-flex justify-content-between align-items-center",
        "IC50 by Tissue Site",
        span(
          actionLink(
            "interpret_box",
            icon_explain,
            class = "me-3 text-decoration-none",
            aria_label = "Explain box plot"
          ),
          popover(
            title = "Add a color variable", placement = "top",
            fa_i("ellipsis"),
            radioButtons(
              "box_color",
              NULL,
              c("none", "Mutations.NRAS", "Mutations.BRAF", "Mutations.PIK3CA", "Mutations.NRAS.PIK3CA"),
              inline = TRUE
            )
          )
        )
      ),
      plotlyOutput("boxplot")
    ),
  
  # 🔍 Data table
  card(
    style = "height: 400px;",
    card_header("IC50 data"),
    reactableOutput("table", height = "100%")
  ),
  )
)

server <- function(input, output, session) {
  # 🔄 Reactive state/computation --------------------------------------------

  current_title <- reactiveVal(NULL)
  current_query <- reactiveVal("")

  # This object must always be passed as the `.ctx` argument to query(), so that
  # tool functions can access the context they need to do their jobs; in this
  # case, the database connection that query() needs.
  ctx <- list(conn = conn)

  # The reactive data frame. Either returns the entire dataset, or filtered by
  # whatever Sidebot decided.
  IC50_data <- reactive({
    sql <- current_query()
    if (is.null(sql) || sql == "") {
      sql <- "SELECT * FROM IC50_expr_cnv;"
    }
    dbGetQuery(conn, sql)
  })



  # 🏷️ Header outputs --------------------------------------------------------

  output$show_title <- renderText({
    current_title()
  })

  output$show_query <- renderText({
    current_query()
  })



  # 🎯 Value box outputs -----------------------------------------------------

  output$total_cell_lines <- renderText({
    nrow(IC50_data())
  })
  
  output$median_IC50 <- renderText({
    x <- median(10^IC50_data()$IC50, na.rm=T)
    paste0(formatC(x, format = "f", digits = 1, big.mark = ","), " nM")
  })

  output$median_AUC <- renderText({
    x <- median(IC50_data()$AUC, na.rm=T)
    formatC(x, format = "f", digits = 1, big.mark = ",")
  })


  # 📊 Scatter plot ----------------------------------------------------------

  scatterplot <- reactive({
    req(nrow(IC50_data()) > 0)

    color <- input$scatter_color

    data <- IC50_data()
    data <- data[!is.na(data$IC50), ]
    # avoid missing data leading to 
    # non-rendering plot

    p <- plot_ly(data, x = ~IC50, y = ~AUC, type = "scatter", mode = "markers")

    if (color != "none") {
    # get number of levels
      n_col <- length(unique(data[[color]]))
      
      # use a palette with enough colors
      pal <- scales::hue_pal()(n_col)   # generates as many distinct hues as needed
      
      p <- plot_ly(
        data,
        x = ~IC50,
        y = ~AUC,
        color = as.formula(paste0("~", color)),
        colors = pal,                 # <- supply custom palette
        type = "scatter",
        mode = "markers"
      )
    }


    p <- p |> add_lines(
      x = ~IC50, y = fitted(loess(AUC ~ IC50, data = data)),
      line = list(color = "rgba(255, 0, 0, 0.5)"),
      name = "LOESS", inherit = FALSE
    )

    p <- p |> layout(showlegend = FALSE)

    return(p)
  })

  output$scatterplot <- renderPlotly({
    scatterplot()
  })

  observeEvent(input$interpret_scatter, {
    explain_plot(chat, scatterplot(), model = gemini_model, .ctx = ctx)
  })
  
  
  # 📊 Box plot -----------------------------------------------------
  boxplot <- reactive({
    req(nrow(IC50_data()) > 0)

    color_by <- input$box_color

    data <- IC50_data()
    data <- data[!is.na(data$IC50), ] 
    # avoid missing data leading to 
    # non-rendering plot

    # Calculate medians for sorting
    data_medians <- data |>
      group_by(SampleCollectionSite) |>
      summarize(median_IC50 = median(IC50, na.rm = TRUE)) |>
      arrange(median_IC50)

    # Get the top 7 sites and shorten the long label
    top_7_sites <- head(data_medians, 7) |>
      pull(SampleCollectionSite)

    shortened_sites <- c(
      "haematopoietic_and_lymphoid_tissue" = "haem/lymph tissue"
    )

    data <- data |>
      filter(SampleCollectionSite %in% top_7_sites) |>
      mutate(
        SampleCollectionSite = factor(
          SampleCollectionSite,
          levels = top_7_sites,
          labels = sapply(top_7_sites, function(site) {
            if (site %in% names(shortened_sites)) {
              return(shortened_sites[site])
            } else {
              return(site)
            }
          })
        ),
        # Safely create the color_by_factor column
        # If 'none' is selected, assign a constant color
        color_by_factor = if (color_by != "none") {
          factor(get(color_by))
        } else {
          factor("black")
        }
      )
    

    # Create the ggplot with boxplot and jittered points
    p <- ggplot(data, aes(x = SampleCollectionSite, y = IC50,
      text = paste0(
        "Cell line: ", Cell_line, "<br>",
        "IC50: ", round(IC50, 3), "<br>",
        "Tissue: ", SampleCollectionSite
      ))) +
      # Box plot layer
      geom_boxplot(outlier.shape=NA) +
      # Individual data points (jitter) layer
      geom_jitter(
        aes(color = color_by_factor),
        width = 0.2, # Adjust jitter width to prevent overlap
        alpha = 0.8
      ) +
      # Adjust theme and labels
      labs(
        x = "Tissue Site",
        y = "IC50",
        color = color_by
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title  = element_text(size = 7)
      )

    # Conditionally hide the legend if 'none' is selected
    if (color_by == "none") {
      p <- p + guides(color = "none")
    }

    return(p)
  })

  output$boxplot <- renderPlotly({
  gg <- boxplot()
  p  <- ggplotly(gg, tooltip = "text")
  # Apply style only to traces of type "box"
  for(i in seq_along(p$x$data)){
    if(p$x$data[[i]]$type == "box"){
      p$x$data[[i]]$boxpoints <- FALSE
    }
  }
  p
  })

  observeEvent(input$interpret_box, {
    explain_plot(chat, boxplot(), model = gemini_model, .ctx = ctx)
  })

  # 🔍 Data table ------------------------------------------------------------
  
  output$table <- renderReactable({
    reactable(IC50_data(),pagination = FALSE, compact = TRUE)
  })

  
  # ✨ Sidebot ✨ -------------------------------------------------------------

  append_output <- function(...) {
    txt <- paste0(...)
    shinychat::chat_append_message(
      "chat",
      list(role = "assistant", content = txt),
      chunk = TRUE,
      operation = "append",
      session = session
    )
  }

  #' Modifies the data presented in the data dashboard, based on the given SQL
  #' query, and also updates the title.
  #' @param query A DuckDB SQL query; must be a SELECT statement.
  #' @param title A title to display at the top of the data dashboard,
  #'   summarizing the intent of the SQL query.
  update_dashboard <- function(query, title) {
    append_output("\n```sql\n", query, "\n```\n\n")

    tryCatch(
      {
        # Try it to see if it errors; if so, the LLM will see the error
        dbGetQuery(conn, query)
      },
      error = function(err) {
        append_output("> Error: ", conditionMessage(err), "\n\n")
        stop(err)
      }
    )

    if (!is.null(query)) {
      current_query(query)
    }
    if (!is.null(title)) {
      current_title(title)
    }
  }

  #' Reset data dashboard
  #' @param query Empty string
  #' @param title Empty string
  reset_dashboard <- function(query, title) {
    append_output("\n```sql\n", query, "\n```\n\n")
    
    if (!is.null(query)) {
      current_query(query)
    }
    if (!is.null(title)) {
      current_title(title)
    }
  }
  
  #' Perform a SQL query on the data, and return the results as JSON.
  #' @param query A DuckDB SQL query; must be a SELECT statement.
  #' @return The results of the query as a JSON string.
  query <- function(query) {
    # Do this before query, in case it errors
    append_output("\n```sql\n", query, "\n```\n\n")

    tryCatch(
      {
        df <- dbGetQuery(conn, query)
      },
      error = function(e) {
        append_output("> Error: ", conditionMessage(e), "\n\n")
        stop(e)
      }
    )
  

    tbl_html <- df_to_html(df, maxrows = 5)
    append_output(tbl_html, "\n\n")

    df |> jsonlite::toJSON(auto_unbox = TRUE)
  }

  # Preload the conversation with the system prompt. These are instructions for
  # the chat model, and must not be shown to the end user.
  chat <- chat_gemini(model = gemini_model,
                      api_key = Sys.getenv("GEMINI_API_KEY"),
                      system_prompt = system_prompt_str)
  chat$register_tool(tool(
    update_dashboard,
    "Modifies the data presented in the data dashboard, based on the given SQL query, and also updates the title.",
    query = type_string("A DuckDB SQL query; must be a SELECT statement."),
    title = type_string("A title to display at the top of the data dashboard, summarizing the intent of the SQL query.")
  ))
  
  chat$register_tool(tool(
    reset_dashboard,
    "Remove all filters and reset dashboard to original state",
    query = type_string(""),
    title = type_string("")
  ))
  
  chat$register_tool(tool(
    query,
    "Perform a SQL query on the data, and return the results as JSON.",
    query = type_string("A DuckDB SQL query; must be a SELECT statement.")
  ))

  # Prepopulate the chat UI with a welcome message that appears to be from the
  # chat model (but is actually hard-coded). This is just for the user, not for
  # the chat model to see.
  chat_append("chat", greeting)

  # Handle user input
  observeEvent(input$chat_user_input, {
    # Add user message to the chat history
    chat_append("chat", chat$stream_async(input$chat_user_input)) %...>% {
      # print(chat)
    }
  })
}

df_to_html <- function(df, maxrows = 5) {
  df_short <- if (nrow(df) > 10) head(df, maxrows) else df

  tbl_html <- capture.output(
    df_short |>
      xtable::xtable() |>
      print(type = "html", include.rownames = FALSE, html.table.attributes = NULL)
  ) |> paste(collapse = "\n")

  if (nrow(df_short) != nrow(df)) {
    rows_notice <- glue::glue("\n\n(Showing only the first {maxrows} rows out of {nrow(df)}.)\n")
  } else {
    rows_notice <- ""
  }

  paste0(tbl_html, "\n", rows_notice)
}

shinyApp(ui, server)
