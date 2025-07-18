library(tidyverse)
library(here)
library(glue)
library(shiny)
library(bslib)
library(ragg)
library(hrbrthemes)

prefix <- "wfm"
height <- 600
# color <- "firebrick2"
color <- "steelblue"

t_wfm <- read_tsv(
  here("02_results", glue("{prefix}_new_wfm.txt")),
  col_names = T,
  show_col_types = F
) |>
  mutate(nS = S * 1e9) |>
  select(-S)

t_eye_property <- read_tsv(
  here("02_results", glue("{prefix}_eye_property.txt")),
  col_names = T,
  show_col_types = F
)

t_eye_t <- read_tsv(
  here("02_results", glue("{prefix}_eye_t.txt")),
  col_names = T,
  show_col_types = F
) |>
  pull(T)

t_eye_tags <- read_tsv(
  here("02_results", glue("{prefix}_eye_db_tags.txt")),
  col_names = T,
  show_col_types = F
) |>
  pull(tags)

t_eye_db <- as.matrix(read_tsv(
  here("02_results", glue("{prefix}_eye_db.txt")),
  col_names = T,
  show_col_types = F
))

t_eye_i <- as.matrix(read_tsv(
  here("02_results", glue("{prefix}_eye_db_i.txt")),
  col_names = T,
  show_col_types = F
))

t_eye_db_long <- tibble(
  T = rep(t_eye_t, ncol(t_eye_db)),
  UIs = rep(1:ncol(t_eye_db), each = nrow(t_eye_db)),
  tags = rep(t_eye_tags, each = nrow(t_eye_db)),
  value = c(t_eye_db)
)

available_tags <- setdiff(unique(t_eye_db_long$tags), NA)

ui <- page_navbar(
  navbar_options = navbar_options(
    bg = "#1c2951", # "#1c2841",  # "#0062cc",
    underline = T
  ),
  title = "EYE Diagram",

  nav_panel(
    "EYE and WFM",
    fluidRow(
      div(
        style = "margin: auto; width: 95%",
        sliderInput(
          "range",
          "UI:",
          min = 1,
          max = 300,
          step = 1,
          value = c(1, 100),
          width = "100%"
        )
      )
    ),
    fluidRow(
      column(
        8,
        plotOutput("wfm")
      ),
      column(
        4,
        plotOutput("eye")
      )
    )
  ),
  nav_panel(
    "Tags",
    selectInput(
      "tags",
      "Tags",
      choices = available_tags,
      selected = c("1T - 1T", "2T - 1T", "3T - 1T"),
      multiple = T,
      width = "100%"
    ),
    plotOutput("eye_by_tags")
  )
)

server <- function(input, output, session) {
  eye_min_ui <- reactive(input$range[1])
  eye_max_ui <- reactive(input$range[2])

  wfm_min_ui <- reactive(ifelse(
    (eye_max_ui() - eye_min_ui()) < 5,
    eye_min_ui(),
    eye_max_ui() - 5
  ))
  wfm_max_ui <- reactive(eye_max_ui() + 5)

  l_eye_i <- reactive(t_eye_i[, wfm_min_ui():wfm_max_ui()])
  wfm_last <- reactive(t_wfm |> slice(t_eye_i[, eye_max_ui()]))
  wfm_all <- reactive(t_wfm |> slice(min(l_eye_i()):max(l_eye_i())))

  # eye_plot <- reactive({
  #   if (eye_max_ui() == eye_min_ui()) {
  #     matplot(
  #       t_eye_t,
  #       t_eye_db[, eye_max_ui()],
  #       type = "l",
  #       lty = 1,
  #       col = alpha("firebrick2"),
  #       xlab = "T",
  #       ylab = "V"
  #       # ylim = c(-2, 2)
  #     )
  #   } else {
  #     matplot(
  #       t_eye_t,
  #       t_eye_db[, eye_min_ui():(eye_max_ui() - 1)],
  #       type = "l",
  #       lty = 1,
  #       col = alpha("gray", 0.6),
  #       xlab = "T",
  #       ylab = "V"
  #     )
  #     lines(
  #       t_eye_t,
  #       t_eye_db[, eye_max_ui()],
  #       type = "l",
  #       lty = 1,
  #       lwd = 2,
  #       col = alpha("firebrick2")
  #     )
  #   }
  # })

  eye_plot <- reactive({
    if (eye_max_ui() == eye_min_ui()) {
      ggplot(
        t_eye_db_long |>
          filter(UIs == eye_max_ui()) |>
          mutate(UIs = factor(UIs)),
        aes(x = T, y = value, color = UIs)
      ) +
        geom_line(linewidth = 2) +
        scale_color_manual(values = c(color)) +
        labs(y = NULL) +
        theme_ipsum_rc() +
        theme(legend.position = "none")
    } else {
      ggplot(
        t_eye_db_long |>
          filter(UIs >= eye_min_ui() & as.numeric(UIs) < eye_max_ui()) |>
          mutate(UIs = factor(UIs)),
        aes(x = T, y = value, color = UIs)
      ) +
        geom_line(alpha = 0.6) +
        scale_color_manual(values = rep("gray", eye_max_ui() - eye_min_ui())) +
        geom_line(
          data = t_eye_db_long |>
            filter(UIs == eye_max_ui()),
          aes(x = T, y = value),
          linewidth = 2,
          color = color
        ) +
        labs(y = NULL) +
        theme_ipsum_rc() +
        theme(legend.position = "none")
    }
  })

  output$eye <- renderPlot(
    {
      eye_plot()
    },
    height = height
  )

  wfm_plot <- reactive({
    ggplot() +
      geom_line(data = wfm_all(), aes(nS, V), color = "grey") +
      geom_point(data = wfm_all(), aes(nS, V), color = "grey") +
      geom_line(data = wfm_last(), aes(nS, V), color = color) +
      geom_point(data = wfm_last(), aes(nS, V), color = color) +
      # scale_y_continuous(limits = c(-2, 2), breaks = c(-2, -1, 0, 1, 2)) +
      theme_ipsum_rc() +
      theme(legend.title = element_blank(), legend.position = "top")
  })

  output$wfm <- renderPlot(
    {
      wfm_plot()
    },
    height = height
  )

  eye_plot_by_tags <- reactive(
    {
      df <- t_eye_db_long |>
        filter(tags %in% input$tags) |>
        mutate(UIs = factor(UIs))
      ggplot(
        df |> mutate(tags = factor(tags, input$tags)),
        aes(x = T, y = value, color = tags, linetype = UIs)
      ) +
        geom_line(alpha = 0.6) +
        scale_linetype_manual(
          values = rep("solid", length(unique(df$UIs))),
          guide = "none"
        ) +
        labs(y = "V") +
        theme_ipsum_rc() +
        scale_color_ft() +
        scale_fill_ft()
    }
  )

  output$eye_by_tags <- renderPlot(
    {
      eye_plot_by_tags()
    },
    height = height
  )
}

shinyApp(ui, server)
