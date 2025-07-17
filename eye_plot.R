# plot eye diagram
#
# - Wei Du
#

library(tidyverse)
library(gsignal)
library(RColorBrewer)
library(plotly)
library(ggmatplot)
library(grid)
library(gridBase)
library(R6)
library(patchwork)
library(here)
library(glue)

p_fig <- "06_fig"
prefix <- "wfm"

df_wfm <- read_tsv(
  here("02_results", glue("{prefix}_new_wfm.txt")),
  show_col_types = FALSE
)
df_wfm_rise <- read_tsv(
  here("02_results", glue("{prefix}_wfm_rise.txt")),
  show_col_types = FALSE
)
df_wfm_fall <- read_tsv(
  here("02_results", glue("{prefix}_wfm_fall.txt")),
  show_col_types = FALSE
)
df_ref <- read_tsv(
  here("02_results", glue("{prefix}_new_ref.txt")),
  show_col_types = FALSE
)
df_ref_raw_with_tags <- read_tsv(
  here(
    "02_results",
    glue("{prefix}_ref_raw_with_tags.txt")
  ),
  show_col_types = FALSE
)
df_eye_property <- read_tsv(
  here(
    "02_results",
    glue("{prefix}_eye_property.txt")
  ),
  show_col_types = FALSE
)
df_eye_t <- read_tsv(
  here("02_results", glue("{prefix}_eye_t.txt")),
  show_col_types = FALSE
)
df_eye_db <- read_tsv(
  here("02_results", glue("{prefix}_eye_db.txt")),
  show_col_types = FALSE
)
df_eye_db_i <- read_tsv(
  here("02_results", glue("{prefix}_eye_db_i.txt")),
  show_col_types = FALSE
)
df_tie <- read_tsv(
  here("02_results", glue("{prefix}_tie.txt")),
  show_col_types = FALSE
)
df_eye_db_transition_tags <- read_tsv(
  here(
    "02_results",
    glue("{prefix}_eye_db_tags.txt")
  ),
  show_col_types = FALSE
)

f10_view_eye_sliced_wfm <- function(tag = NA, cnt = 4) {
  message("f10: view wfm and sliced wfm for eye")
  if (is.na(tag)) {
    n <- min(cnt, ncol(df_eye_db))
    x_start_min <- min(df_eye_db_i[, 1:n])
    x_start_max <- max(df_eye_db_i[, 1:n])
    x_end_min <- min(df_eye_db_i[,
      ((ncol(df_eye_db_i) - n + 1):ncol(df_eye_db_i))
    ])
    x_end_max <- max(df_eye_db_i[,
      ((ncol(df_eye_db_i) - n + 1):ncol(df_eye_db_i))
    ])
    t_start <- tibble(
      `#` = seq(x_start_min, x_start_max),
      v = df_wfm$V[x_start_min:x_start_max],
      type = "wfm"
    )
    t_start <- rbind(
      t_start,
      tibble(
        `#` = seq(x_start_min, x_start_max),
        v = df_ref$V[x_start_min:x_start_max],
        type = "ref"
      )
    )
    levels <- c("wfm", "ref")
    for (i in 1:n) {
      t_start <- rbind(
        t_start,
        tibble(`#` = df_eye_db_i[[i]], v = df_eye_db[[i]], type = str_c("#", i))
      )
      levels <- c(levels, str_c("#", i))
    }
    t_start$type <- factor(t_start$type, levels = levels)
    t_end <- tibble(
      `#` = seq(x_end_min, x_end_max),
      v = df_wfm$V[x_end_min:x_end_max],
      type = "wfm"
    )
    t_end <- rbind(
      t_end,
      tibble(
        `#` = seq(x_end_min, x_end_max),
        v = df_ref$V[x_end_min:x_end_max],
        type = "ref"
      )
    )
    levels <- c("wfm", "ref")
    for (i in n:1) {
      t_end <- rbind(
        t_end,
        tibble(
          `#` = df_eye_db_i[[(ncol(df_eye_db_i) - i + 1)]],
          v = df_eye_db[[(ncol(df_eye_db_i) - i + 1)]],
          type = str_c("#-", abs(-i + 1))
        )
      )
      levels <- c(levels, str_c("#-", abs(-i + 1)))
    }
    t_end$type <- factor(t_end$type, levels = levels)
  } else {
    tag_index <- which(df_eye_db_transition_tags$tags == tag)

    stopifnot(
      "!!! tag not found, expected format: 3T - 1T !!!" = length(tag_index) > 0
    )

    loc_eye_db <- df_eye_db[, tag_index]
    loc_eye_db_i <- df_eye_db_i[, tag_index]

    n <- min(cnt, ncol(loc_eye_db))

    x_start_min <- min(loc_eye_db_i[, 1:n])
    x_start_max <- max(loc_eye_db_i[, 1:n])
    x_end_min <- min(loc_eye_db_i[,
      ((ncol(loc_eye_db_i) - n + 1):ncol(loc_eye_db_i))
    ])
    x_end_max <- max(loc_eye_db_i[,
      ((ncol(loc_eye_db_i) - n + 1):ncol(loc_eye_db_i))
    ])
    t_start <- tibble(
      `#` = seq(x_start_min, x_start_max),
      v = df_wfm$V[x_start_min:x_start_max],
      type = "wfm"
    )
    t_start <- rbind(
      t_start,
      tibble(
        `#` = seq(x_start_min, x_start_max),
        v = df_ref$V[x_start_min:x_start_max],
        type = "ref"
      )
    )
    levels <- c("wfm", "ref")
    for (i in 1:n) {
      t_start <- rbind(
        t_start,
        tibble(
          `#` = loc_eye_db_i[[i]],
          v = loc_eye_db[[i]],
          type = str_c("#", i)
        )
      )
      levels <- c(levels, str_c("#", i))
    }
    t_start$type <- factor(t_start$type, levels = levels)
    t_end <- tibble(
      `#` = seq(x_end_min, x_end_max),
      v = df_wfm$V[x_end_min:x_end_max],
      type = "wfm"
    )
    t_end <- rbind(
      t_end,
      tibble(
        `#` = seq(x_end_min, x_end_max),
        v = df_ref$V[x_end_min:x_end_max],
        type = "ref"
      )
    )
    levels <- c("wfm", "ref")
    for (i in n:1) {
      t_end <- rbind(
        t_end,
        tibble(
          `#` = loc_eye_db_i[[(ncol(loc_eye_db_i) - i + 1)]],
          v = loc_eye_db[[(ncol(loc_eye_db_i) - i + 1)]],
          type = str_c("#-", abs(-i + 1))
        )
      )
      levels <- c(levels, str_c("#-", abs(-i + 1)))
    }
    t_end$type <- factor(t_end$type, levels = levels)
  }
  p_start <- ggplot(t_start, aes(x = `#`, y = v, color = type, shape = type)) +
    geom_line(data = t_start[!t_start$type == "ref", ]) +
    geom_step(data = t_start[t_start$type == "ref", ]) +
    # geom_point() +
    # facet_grid(type ~ ., scale = "free") +
    facet_grid(type ~ .) +
    labs(title = "start") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    # theme_bw() +
    theme(
      strip.text.y = element_blank(),
      legend.title = element_blank()
    )
  # print(p_start)
  ggp_start <- ggplotly(p_start)
  print(ggp_start)
  wd <- getwd()
  setwd(here(p_fig))
  html_name <- glue("{prefix}_f10_view_sliced_wfm_start.html")
  htmlwidgets::saveWidget(ggp_start, html_name)
  setwd(wd)
  message("save to ", html_name)

  p_end <- ggplot(t_end, aes(x = `#`, y = v, color = type, shape = type)) +
    geom_line(data = t_end[!t_end$type == "ref", ]) +
    geom_step(data = t_end[t_end$type == "ref", ]) +
    # geom_point() +
    # facet_grid(type ~ ., scale = "free") +
    facet_grid(type ~ .) +
    labs(title = "end") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    # theme_bw() +
    theme(
      strip.text.y = element_blank(),
      legend.title = element_blank()
    )
  # print(p_end)
  ggp_end <- ggplotly(p_end)
  print(ggp_end)
  wd <- getwd()
  setwd(here(p_fig))
  html_name <- glue("{prefix}_f10_view_sliced_wfm_end.html")
  htmlwidgets::saveWidget(ggp_end, html_name)
  setwd(wd)
  message("save to ", html_name)
}

f10_ggplot_zoomin_wfm_ref_tie <- function(
  max_lines = NULL,
  color = "steelblue",
  bins = 50
) {
  t_max <- (df_eye_property$OS *
    nrow(df_ref_raw_with_tags) *
    1 /
    df_eye_property$Fs) /
    10 *
    1e9 # nS

  data_wfm <- rbind(
    tibble(x = df_wfm$S, y = df_wfm$V, type = "wfm (zoom in)"),
    tibble(
      x = df_wfm_rise$S,
      y = df_wfm_rise$V,
      type = "zero crossing - rising edge"
    ),
    tibble(
      x = df_wfm_fall$S,
      y = df_wfm_fall$V,
      type = "zero crossing - falling edge"
    )
  ) %>%
    mutate(x = x * 1e9) %>%
    dplyr::filter(x <= t_max)

  vlines <- df_tie$S * 1e9
  vlines <- vlines[vlines <= t_max]

  p_wfm <- ggplot(data_wfm, aes(x, y, color = type, shape = type)) +
    geom_line(data = data_wfm[data_wfm$type == "wfm (zoom in)", ]) +
    geom_point(data = data_wfm[!data_wfm$type == "wfm (zoom in)", ]) +
    geom_vline(xintercept = vlines) +
    scale_x_continuous(limits = c(0, t_max)) +
    scale_shape_manual(
      values = c(
        `wfm (zoom in)` = 4,
        `zero crossing - rising edge` = 4,
        `zero crossing - falling edge` = 3
      )
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      title = str_c(
        "total wfm length: ",
        nrow(df_wfm),
        "points, ",
        max(df_wfm$S) * 1e9,
        "nS"
      ),
      x = "nS",
      y = "V"
    ) +
    theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top")

  # data_ref <- tibble(x = (seq_along(self$ref_v) - 1) * (1 / self$Fs) * 1e12, y = self$ref_v, type = "genie (zoom in)")
  data_ref <- tibble(x = df_ref$S * 1e9, y = df_ref$V, type = "genie (zoom in)")
  p_ref <- ggplot(
    data_ref %>% dplyr::filter(x <= t_max),
    aes(x, y, color = type)
  ) +
    geom_step() +
    scale_x_continuous(limits = c(0, t_max)) +
    scale_y_continuous(breaks = c(-1, 0, 1)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "nS", y = "") +
    theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top")

  data_tie <- tibble(x = df_tie$S * 1e9, y = df_tie$pS, type = "tie (zoom in)")
  p_tie <- ggplot(
    data_tie %>% dplyr::filter(x <= t_max),
    aes(x, y, color = type)
  ) +
    geom_line() +
    geom_point() +
    scale_x_continuous(limits = c(0, t_max)) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(x = "nS", y = "pS") +
    theme_bw() +
    theme(legend.title = element_blank(), legend.position = "top")

  if (is.null(max_lines)) {
    data_tie_for_hist <- tibble(
      x = df_tie$S * 1e9,
      y = df_tie$pS,
      type = "tie (zoom in)"
    )
  } else {
    tie_index <- which(
      df_tie$S <=
        (tail(df_eye_db_i[[max_lines]], 1) - 1) * 1 / df_eye_property$Fs
    )
    data_tie_for_hist <- tibble(
      x = df_tie$S[tie_index] * 1e9,
      y = df_tie$pS[tie_index],
      type = "tie (zoom in)"
    )
  }

  p_tie_hist <- ggplot(data_tie_for_hist, aes(y)) +
    geom_histogram(bins = bins, color = color, fill = color) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(title = "TIE", x = "pS") +
    theme_bw()

  return(list(
    "p_wfm" = p_wfm,
    "p_ref" = p_ref,
    "p_tie" = p_tie,
    "p_tie_hist" = p_tie_hist
  ))
}

f10_ggplot_eye <- function(
  loc_eye_db = NULL,
  index = NULL,
  max_lines = 500,
  color = "steelblue",
  alpha = 0.2,
  bins = 50
) {
  # loc_eye_db is the matrix
  if (is.null(loc_eye_db)) {
    loc_eye_db <- df_eye_db
  }
  if (is.null(index)) {
    if (is.null(df_eye_t)) {
      index <- 1:nrow(loc_eye_db)
    } else {
      index <- df_eye_t$T
    }
  }

  # colnames(loc_eye_db) <- 1:ncol(loc_eye_db)
  # t <- as_tibble(loc_eye_db, .name_repair = "unique") %>% mutate(index = index)
  if (is.null(max_lines)) {
    max_lines <- ncol(loc_eye_db)
  } else {
    max_lines <- min(max_lines, ncol(loc_eye_db))
  }
  loc_eye_db <- loc_eye_db %>% mutate(index = index)

  p <- ggplot(loc_eye_db, aes(x = index))
  for (k in 1:max_lines) {
    p <- p +
      geom_line(aes(y = .data[[names(t)[k]]]), color = color, alpha = alpha)
  }
  p <- p +
    labs(
      title = str_c(
        sprintf("%.3f", df_eye_property$`datarate(Gbps)`),
        "Gbps, ",
        max_lines,
        "UI"
      ),
      x = "T",
      y = "V"
    ) +
    theme_bw()
  return(p)
}

f10_ggplot_eye_and_wfm <- function(
  loc_eye_db = NULL,
  index = NULL,
  max_lines = 500,
  color = "steelblue",
  alpha = 0.2,
  bins = 50
) {
  l <- f10_ggplot_zoomin_wfm_ref_tie(max_lines, color, bins)

  p_eye <- f10_ggplot_eye(loc_eye_db, index, max_lines, color, alpha, bins)

  design <- c(
    area(t = 1, l = 1, b = 1, r = 2),
    area(t = 2, l = 1, b = 2, r = 2),
    area(t = 3, l = 1, b = 3, r = 2),
    area(t = 1, l = 3, b = 2, r = 3),
    area(t = 3, l = 3, b = 3, r = 3)
  )
  p <- l$p_wfm +
    l$p_ref +
    l$p_tie +
    p_eye +
    l$p_tie_hist +
    plot_layout(design = design)
  print(p)
  fig_name <- here(p_fig, "f10_ggplot_eye.png")
  ggsave(fig_name, p, width = 10, height = 6)
  message("save to ", fig_name)
  return(p)
}

f10_ggmatplot_eye <- function(
  loc_eye_db = NULL,
  index = NULL,
  max_lines = 500,
  color = "steelblue",
  alpha = 0.2,
  bins = 50
) {
  # loc_eye_db is the matrix
  if (is.null(loc_eye_db)) {
    loc_eye_db <- df_eye_db
  }
  if (is.null(index)) {
    if (is.null(df_eye_t)) {
      index <- 1:nrow(loc_eye_db)
    } else {
      index <- df_eye_t$T
    }
  }

  if (is.null(max_lines)) {
    max_lines <- ncol(loc_eye_db)
  } else {
    max_lines <- min(max_lines, ncol(loc_eye_db))
  }

  p <- ggmatplot(
    index,
    df_eye_db[, 1:max_lines],
    plot_type = "line",
    linetype = 1,
    color = color,
    alpha = alpha,
    show.legend = F,
    xlab = "T",
    ylab = "V"
  ) +
    labs(
      title = str_c(
        sprintf("%.3f", df_eye_property$`datarate(Gbps)`),
        "Gbps, ",
        max_lines,
        "UI"
      )
      # ylab = "V"
    ) +
    theme_bw()
  return(p)
}

f10_ggmatplot_eye_and_wfm <- function(
  loc_eye_db = NULL,
  index = NULL,
  max_lines = 500,
  color = "steelblue",
  alpha = 0.2,
  bins = 50
) {
  l <- f10_ggplot_zoomin_wfm_ref_tie(max_lines, color, bins)

  p_eye <- f10_ggmatplot_eye(loc_eye_db, index, max_lines, color, alpha, bins)

  design <- c(
    area(t = 1, l = 1, b = 1, r = 2),
    area(t = 2, l = 1, b = 2, r = 2),
    area(t = 3, l = 1, b = 3, r = 2),
    area(t = 1, l = 3, b = 2, r = 3),
    area(t = 3, l = 3, b = 3, r = 3)
  )
  p <- l$p_wfm +
    l$p_ref +
    l$p_tie +
    p_eye +
    l$p_tie_hist +
    plot_layout(design = design)
  # print(p)
  fig_name <- here(p_fig, "f10_ggmatplot_eye.png")
  ggsave(fig_name, p, width = 10, height = 6)
  message("save to ", fig_name)
  return(p)
}

f10_matplot_eye <- function(
  loc_eye_db = NULL,
  index = NULL,
  max_lines = NULL,
  color = "steelblue",
  alpha = 0.2,
  bins = 50
) {
  # loc_eye_db is the matrix
  if (is.null(loc_eye_db)) {
    loc_eye_db <- df_eye_db
  }
  if (is.null(index)) {
    if (is.null(df_eye_t)) {
      index <- 1:nrow(loc_eye_db)
    } else {
      index <- df_eye_t$T
    }
  }
  if (is.null(max_lines)) {
    max_lines <- ncol(loc_eye_db)
  } else {
    max_lines <- min(max_lines, ncol(loc_eye_db))
  }

  matplot(
    index,
    df_eye_db[, 1:max_lines],
    type = "l",
    lty = 1,
    col = alpha(color, alpha),
    xlab = "T",
    ylab = "V",
    main = str_c(
      sprintf("%.3f", df_eye_property$`datarate(Gbps)`),
      "Gbps, ",
      max_lines,
      "UI"
    )
  )
}

f10_matplot_eye_and_wfm <- function(
  loc_eye_db = NULL,
  index = NULL,
  max_lines = NULL,
  color = "steelblue",
  alpha = 0.2,
  bins = 50
) {
  l <- f10_ggplot_zoomin_wfm_ref_tie(max_lines, color, bins)

  # png(str_c(self$p_fig, "f10_matplot_eye.png"), width = 12, height = 6, units = "in", res = 300)

  layout(
    matrix(c(2, 2, 1, 3, 3, 1, 4, 4, 5), nrow = 3, ncol = 3, byrow = T),
    heights = c(2, 2, 2),
    widths = c(2, 2, 2)
  )

  f10_matplot_eye(loc_eye_db, index, max_lines, color, alpha, bins)

  plot.new()
  vps <- baseViewports()
  pushViewport(vps$figure)
  vp <- plotViewport(c(0, 0, 0, 0))
  print(l$p_wfm, vp = vp)
  popViewport()

  plot.new()
  vps <- baseViewports()
  pushViewport(vps$figure)
  vp <- plotViewport(c(0, 0, 0, 0))
  print(l$p_ref, vp = vp)
  popViewport()

  plot.new()
  vps <- baseViewports()
  pushViewport(vps$figure)
  vp <- plotViewport(c(0, 0, 0, 0))
  print(l$p_tie, vp = vp)
  popViewport()

  plot.new()
  vps <- baseViewports()
  pushViewport(vps$figure)
  vp <- plotViewport(c(0, 0, 0, 0))
  print(l$p_tie_hist, vp = vp)
  popViewport()

  par(mfrow = c(1, 1)) # reset

  # dev.off()
}

f10_view_eye_sliced_wfm()
# l <- f10_ggplot_zoomin_wfm_ref_tie()
# p <- f10_ggplot_eye_and_wfm(max_lines = 200)
# p <- f10_ggmatplot_eye_and_wfm(max_lines = 100)
# f10_matplot_eye_and_wfm(max_lines = 100)
f10_matplot_eye_and_wfm(max_lines = NULL)

# expected in rstudio but distortion in positron
dev.print(pdf, file = here(p_fig, glue("{prefix}_f10_matplot_eye.pdf")))
