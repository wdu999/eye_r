# parse wfm file saved from scope
# and reference file
# and genrate eye diagram, jitter
#
# - Wei Du
#

library(tidyverse)
library(plotly)
library(gsignal)
library(here)
library(glue)

p_in <- "01_data"
p_out <- "02_results"
p_fig <- "06_fig"
wfm_file <- "wfm.dat"
ref_file <- "prbs7.txt"

prefix <- tools::file_path_sans_ext(basename(wfm_file))

Fs <- 40e9
F <- 1e9
OS <- Fs / F

prbs_length <- NA
df_ref_raw_with_tags <- NA
eye_window_length_l <- 0.8 # first 0.8T in eye window
eye_window_length_r <- 0.8 # second 0.8T in eye window

df_wfm_raw <- NULL
df_wfm <- NULL
df_wfm_rise <- NULL
df_wfm_fall <- NULL

df_ref_raw <- NULL
df_ref <- NULL
df_ref_rise <- NULL
df_ref_fall <- NULL

df_tie <- NULL

eye_samples_per_T <- NULL
eye_samples_per_Window <- NULL

df_eye_t <- NULL
df_eye_db <- NULL
df_eye_db_i <- NULL
df_eye_db_transition_tags <- NULL
df_eye_jitter <- NULL

df_eye_property <- NULL

format_message <- function(s) {
  message("")
  message(str_flatten(rep("-", str_count(s))))
  message(s)
  message(str_flatten(rep("-", str_count(s))))
  message("")
}

linear_interp <- function(x0, y0, x1, y1, y) {
  return(x0 + (y - y0) * (x1 - x0) / (y1 - y0))
}

interleave <- function(l1, l2) {
  n1 <- length(l1)
  n2 <- length(l2)

  if (n1 != n2) {
    message("length not equal")
  }

  l <- c()
  for (i in 1:max(n1, n2)) {
    if (i <= n1) {
      l <- c(l, l1[i])
    }
    if (i <= n2) l <- c(l, l2[i])
  }
  return(l)
}

zero_crossing <- function(
  x,
  y,
  thr = 0,
  find_rising_edge = TRUE,
  interp = FALSE
) {
  if (find_rising_edge) {
    zc_i <- seq_along(y)[(y[2:length(y)] >= thr) & (y[1:(length(y) - 1)] < thr)]
  } else {
    zc_i <- seq_along(y)[(y[2:length(y)] <= thr) & (y[1:(length(y) - 1)] > thr)]
  }

  if (interp) {
    zc_x <- c()
    for (i in zc_i) {
      zc_x <- c(zc_x, linear_interp(x[i], y[i], x[i + 1], y[i + 1], thr))
    }
    zc_y <- rep(thr, length(zc_i))
  } else {
    zc_x <- x[zc_i]
    zc_y <- y[zc_i]
  }
  return(tibble(i = zc_i, x = zc_x, y = zc_y))
}

zero_crossing_wfm <- function(
  ref_i,
  wfm_t,
  wfm_v,
  find_rising_edge = TRUE,
  interp = TRUE
) {
  d <- 10 # find within a small range
  zc_i <- c()
  zc_x <- c()
  zc_y <- c()

  for (i in ref_i) {
    rng <- seq(max(i - d, 0), min(i + d, length(wfm_v)))
    t <- zero_crossing(wfm_t[rng], wfm_v[rng], 0, find_rising_edge, interp)

    if (nrow(t) == 0) {
      message("\n!!! unable to find zero crossing !!!\n")
      next
    }

    wi <- t$i
    wx <- t$x
    wy <- t$y

    if (length(wi) == 1 && length(wx) == 1 && length(wy) == 1) {
      zc_i <- c(zc_i, rng[wi])
      zc_x <- c(zc_x, wx)
      zc_y <- c(zc_y, wy)
    } else {
      zc_i <- c(zc_i, NA)
      zc_x <- c(zc_x, NA)
      zc_y <- c(zc_y, NA)
    }
  }
  return(tibble(i = zc_i, x = zc_x, y = zc_y))
}

f1_load_raw_wfm <- function() {
  format_message(str_c("f1: load raw wfm ", here(p_in, wfm_file)))
  df_wfm_raw <<- read_delim(
    here(p_in, wfm_file),
    col_names = c("S", "V"),
    col_types = list(col_double(), col_double()),
    delim = " "
  )
}

f1_load_raw_ref <- function() {
  format_message(str_c("f1: load raw ref ", here(p_in, ref_file)))
  df_ref_raw <<- read_delim(
    here(p_in, ref_file),
    col_names = c("nrz"),
    # col_types = list(col_double(), col_integer(), col_integer()),
    delim = ","
  )
}

f1_view_raw_wfm <- function(n = NULL) {
  if (is.null(n)) {
    if (is.null(df_ref_raw)) {
      n <- 1000
    } else {
      n <- nrow(df_ref_raw) * OS * 3
    }
  }

  p_wfm <- ggplot(
    df_wfm_raw[1:n, ],
    aes(S, V)
  ) +
    geom_line() +
    # geom_point() +
    labs(title = "raw wfm")
  print(ggplotly(p_wfm))
}

f1_view_raw_ref <- function() {
  p_ref <- ggplot(
    df_ref_raw %>% mutate(`#` = seq_along(df_ref_raw$nrz)),
    aes(`#`, nrz)
  ) +
    geom_step() +
    labs(title = "raw ref")
  print(ggplotly(p_ref))
}

f2_upsample <- function() {
  format_message("f2: upsample wfm")
  if ((Fs %% F) > 0) {
    message("F = ", F, " Fs = ", Fs, " fractional sampling, do interpolation")

    wfm_v_old <- df_wfm_raw$V
    wfm_t_old <- (seq_along(wfm_v_old) - 1) * (1 / Fs)
    Fs_old <- Fs

    OS <<- ceiling(Fs_old / F / 10) * 10 # new os
    Fs <<- F * OS # new Fs
    up <- Fs / Fs_old

    message("F = ", F, " Fs = ", Fs, "up by ", up, "x, new OS ", OS)

    n <- length(wfm_v_old)
    n_new <- as.integer(n * up)

    x <- pracma::linspace(0, 1, n)
    x_new <- pracma::linspace(0, 1, n_new)

    f <- splinefun(x, wfm_v_old)

    df_wfm <<- tibble(S = (1:n_new - 1) * 1 / Fs, V = f(x_new))
    file_new_wfm <- here(p_out, glue("{prefix}_new_wfm.txt"))
    write_tsv(df_wfm, file_new_wfm)
    message("save to ", file_new_wfm)

    cnt <- 500
    dfp <- rbind(
      tibble(S = wfm_t_old[1:cnt], V = wfm_v_old[1:cnt], type = "original"),
      tibble(
        S = df_wfm$S[1:(cnt * up)],
        V = df_wfm$V[1:(cnt * up)],
        type = "upsampled"
      )
    )

    p <- ggplot(
      dfp,
      aes(S, V, color = type, shape = type)
    ) +
      geom_line() +
      geom_point() +
      scale_shape_manual(values = c(3, 4)) +
      # scale_shape_manual(values = c(20, 3)) +
      # scale_color_brewer(palette = "Dark2") +
      # scale_fill_brewer(palette = "Dark2") +
      scale_color_manual(values = c("#003049", "#d62828", "#f77f00")) +
      scale_fill_manual(values = c("#003049", "#d62828", "#f77f00")) +
      labs(title = "upsample wfm") +
      theme(legend.title = element_blank())
    ggp <- ggplotly(p)
    print(ggp)
    html_name <- glue("{prefix}_f2_upsample.html")
    wd <- getwd()
    setwd(here(p_fig))
    htmlwidgets::saveWidget(ggp, html_name)
    setwd(wd)
    message("save to ", html_name)
  } else {
    message("integer sampling, no interpolation")

    OS <<- Fs / F
    t0 <- df_wfm_raw$S[1]
    df_wfm <<- df_wfm_raw %>% mutate(S = S - t0)

    # save to new_wfm.txt anyway
    file_new_wfm <- here(p_out, glue("{prefix}_new_wfm.txt"))
    write_tsv(df_wfm, file_new_wfm)
    message("still save to ", file_new_wfm)
  }
}

f3_detect_pattern <- function() {
  format_message("f3: detect prbs length and update transition tags")
  ref_one_cycle <- df_ref_raw$nrz

  c <- ccf(
    ref_one_cycle,
    ref_one_cycle,
    plot = TRUE,
    lag.max = (length(ref_one_cycle) - 1)
  )
  c_max <- max(c$acf)
  c_max_lag <- c$lag[which(c$acf == c_max)] # it is 0
  high_corr_lags <- c$lag[which(c$acf > 0.2)]
  prbs_length <<- high_corr_lags[which(high_corr_lags > c_max_lag)][1]
  message("prbs length ", prbs_length)

  ref_patterns <- rep(ref_one_cycle[1:prbs_length], 3)

  pats <- rep(ref_one_cycle * 2 - 1, 3)
  zc <- zerocrossing(seq_along(pats), pats) - 0.5 # the index before zero crossing
  zc_diff <- diff(zc)
  tags <- seq_along(pats) * NA
  for (i in seq_along(zc)) {
    if ((i == 1) | i == length(zc)) {
      tags[zc[i]] <- NA
    } else {
      tags[zc[i]] <- str_c(zc_diff[i - 1], "T - ", zc_diff[i], "T")
    }
  }
  ref_one_prbs_with_tags <- set_names(
    pats[1:prbs_length],
    tags[(prbs_length + 1):(prbs_length * 2)]
  )
  ref_raw_with_tags <- rep(
    ref_one_prbs_with_tags,
    ceiling(length(ref_one_cycle) / prbs_length)
  )[1:length(ref_one_cycle)]
  df_ref_raw_with_tags <<- tibble(
    tags = names(ref_raw_with_tags),
    V = ref_raw_with_tags
  )
  file_ref_raw_with_tags <- here(p_out, glue("{prefix}_ref_raw_with_tags.txt"))
  write_tsv(df_ref_raw_with_tags, file_ref_raw_with_tags)
  message("save to ", file_ref_raw_with_tags)
  message("update raw ref with tags of nT transition")

  dfp <- rbind(
    tibble(x = seq_along(pats), y = pats, type = "ref"),
    tibble(x = zc, y = 0, type = "zc")
  )
  p <- ggplot(dfp, aes(x, y, color = type)) +
    geom_step() +
    geom_point() +
    labs(title = "3x ref raw and zero crossings") +
    scale_color_manual(values = c("#003049", "#d62828", "#f77f00")) +
    scale_fill_manual(values = c("#003049", "#d62828", "#f77f00")) +
    theme(legend.title = element_blank())
  ggp <- ggplotly(p)
  print(ggp)
  html_name <- glue("{prefix}_f3_detect_pattern.html")
  wd <- getwd()
  setwd(here(p_fig))
  htmlwidgets::saveWidget(ggp, html_name)
  setwd(wd)
  message("save to ", html_name)
}

f4_roll_ref <- function() {
  format_message("f4: roll pat to match wfm")
  tags <- df_ref_raw_with_tags$tags
  index_with_tags <- seq_along(df_ref_raw_with_tags$tags)[!is.na(tags)]

  ref_one_cycle <- rep(df_ref_raw_with_tags$V, each = OS)
  names(ref_one_cycle) <- NA
  names(ref_one_cycle)[index_with_tags * OS] <- tags[!is.na(tags)]
  message("update upsampled ref with tags of nT transition")

  wfm_one_cycle <- df_wfm$V[1:length(ref_one_cycle)]
  ref_one_cycle_invert <- -1 * ref_one_cycle

  c <- ccf(
    wfm_one_cycle,
    ref_one_cycle,
    lag.max = (length(wfm_one_cycle) - 1),
    plot = FALSE
  )
  c_max <- max(c$acf)
  c_min <- min(c$acf)
  c_max_lag <- c$lag[which(c$acf == max(c$acf))]
  c_min_lag <- c$lag[which(c$acf == min(c$acf))]

  if (c_max > abs(c_min)) {
    i <- c_max_lag
    ref_one_cycle_rolled <- c(tail(ref_one_cycle, i), head(ref_one_cycle, -i))
    message("ref samples correlate wfm, roll ref ", i, " to align with wfm")
  } else {
    i <- c_min_lag
    ref_one_cycle_rolled <- c(
      tail(rev(ref_one_cycle_invert), i),
      head(rev(ref_one_cycle_invert), -i)
    )
    message(
      "inverted ref samples correlate wfm, roll ref ",
      i,
      " to align with wfm"
    )
  }

  v <- rep(
    ref_one_cycle_rolled,
    ceiling(length(df_wfm$V) / length(ref_one_cycle_rolled))
  )[1:length(df_wfm$V)]
  t <- (seq_along(v) - 1) * (1 / Fs)
  df_ref <<- tibble(
    S = t,
    V = v,
    tags = names(v)
  )
  file_new_ref <- here(p_out, glue("{prefix}_new_ref.txt"))
  write_tsv(df_ref, file_new_ref)
  message("save to ", file_new_ref)

  scale <- 1 / mean(wfm_one_cycle[wfm_one_cycle > 0])

  dfp <- rbind(
    tibble(
      `#` = seq_along(wfm_one_cycle),
      y = wfm_one_cycle * scale,
      type = "wfm (normalized to ref)"
    ),
    tibble(
      `#` = seq_along(wfm_one_cycle),
      y = ref_one_cycle_rolled,
      type = "ref"
    )
  )
  p <- ggplot(
    dfp,
    aes(`#`, y, color = type, shape = type)
  ) +
    geom_line() +
    geom_point() +
    scale_shape_manual(values = c(4, 4)) +
    # scale_color_brewer(palette = "Dark2") +
    # scale_fill_brewer(palette = "Dark2") +
    scale_color_manual(values = c("#003049", "#d62828", "#f77f00")) +
    scale_fill_manual(values = c("#003049", "#d62828", "#f77f00")) +
    labs(title = "one cycle wfm and rolled ref") +
    theme(legend.title = element_blank())
  ggp <- ggplotly(p)
  print(ggp)
  html_name <- glue("{prefix}_f4_roll_ref.html")
  wd <- getwd()
  setwd(here(p_fig))
  htmlwidgets::saveWidget(ggp, html_name)
  setwd(wd)
  message("save to ", html_name)
}
f5_remove_leading_trailing_samples <- function() {
  # remove leading and/or trailing samples
  # to solve the problem that ref/wfm may have more zero crossing than the other one
  format_message("f5: remove leading and/or trailing samples")

  n <- OS / 2

  x <- 1:n
  ref <- df_ref$V[x]
  wfm <- df_wfm$V[x]

  if ((sign(ref[1]) * sign(ref[n])) != (sign(wfm[1]) * sign(wfm[n]))) {
    message(
      "remove leading 0.5T samples due to ref or wfm have more zero crossing than the other"
    )
    matplot(x, cbind(ref, wfm), type = "b", pch = 1, col = 1:2)
    legend("topright", legend = c("ref", "wfm"), pch = 1, col = 1:2)
    df_ref <<- df_ref[(n + 1):nrow(df_ref), ]
    df_wfm <<- df_wfm[(n + 1):nrow(df_wfm), ]
  }

  x <- (nrow(df_ref) - n + 1):nrow(df_ref)
  ref <- df_ref$V[x]
  wfm <- df_wfm$V[x]

  if ((sign(ref[1]) * sign(ref[n])) != (sign(wfm[1]) * sign(wfm[n]))) {
    message(
      "remove trailing 0.5T samples due to ref or wfm have more zero crossing than the other"
    )
    matplot(x, cbind(ref, wfm), type = "b", pch = 1, col = 1:2)
    legend("topright", legend = c("ref", "wfm"), pch = 1, col = 1:2)
    df_ref <<- df_ref[1:(nrow(df_ref) - n), ]
    df_wfm <<- df_wfm[1:(nrow(df_wfm) - n), ]
  }

  # overwrite new_wfm.txt
  file_new_wfm <- here(p_out, glue("{prefix}_new_wfm.txt"))
  write_tsv(df_wfm, file_new_wfm)
  message("overwrite ", file_new_wfm)

  file_new_ref <- here(p_out, glue("{prefix}_new_ref.txt"))
  write_tsv(df_ref, file_new_ref)
  message("overwrite ", file_new_ref)
}
f6_find_zero_crossing <- function() {
  format_message("f6: find zero crossing")

  t <- zero_crossing(
    df_ref$S,
    df_ref$V,
    thr = 0,
    find_rising_edge = TRUE,
    interp = FALSE
  )
  df_ref_rise <<- tibble(
    i = t$i,
    S = t$x,
    V = t$y,
  )

  t <- zero_crossing(
    df_ref$S,
    df_ref$V,
    thr = 0,
    find_rising_edge = FALSE,
    interp = FALSE
  )
  df_ref_fall <<- tibble(
    i = t$i,
    S = t$x,
    V = t$y,
  )

  message("find ref zero crossing, seperate rising and falling edges")

  t <- zero_crossing_wfm(
    df_ref_rise$i,
    df_wfm$S,
    df_wfm$V,
    find_rising_edge = TRUE,
    interp = TRUE
  )
  df_wfm_rise <<- tibble(
    i = t$i,
    S = t$x,
    V = t$y,
  )

  t <- zero_crossing_wfm(
    df_ref_fall$i,
    df_wfm$S,
    df_wfm$V,
    find_rising_edge = FALSE,
    interp = TRUE
  )
  df_wfm_fall <<- tibble(
    i = t$i,
    S = t$x,
    V = t$y,
  )

  file_ref_rise <- here(p_out, glue("{prefix}_ref_rise.txt"))
  write_tsv(df_ref_rise, file_ref_rise)
  message("save to ", file_ref_rise)

  file_ref_fall <- here(p_out, glue("{prefix}_ref_fall.txt"))
  write_tsv(df_ref_fall, file_ref_fall)
  message("save to ", file_ref_fall)

  file_wfm_rise <- here(p_out, glue("{prefix}_wfm_rise.txt"))
  write_tsv(df_wfm_rise, file_wfm_rise)
  message("save to ", file_wfm_rise)

  file_wfm_fall <- here(p_out, glue("{prefix}_wfm_fall.txt"))
  write_tsv(df_wfm_fall, file_wfm_fall)
  message("save to ", file_wfm_fall)

  message(
    "find wfm zero crossing, seperate rising and falling edges, refer ref crossing edges"
  )
}

f7_gen_tie <- function() {
  format_message("f7: generate tie")
  is_rising_edge_leading <- (df_ref_rise$i[1] <= df_ref_fall$i[1])
  first_zero_crossing_i <- min(df_ref_rise$i[1], df_ref_fall$i[1])
  first_zero_cross_wfm_t <- min(df_wfm_rise$S[1], df_wfm_fall$S[1])
  message("determin first wfm zero crossing location and time")

  ref_clk_rising_edges_t <- (df_ref_rise$i - first_zero_crossing_i) *
    1 /
    Fs +
    first_zero_cross_wfm_t
  ref_clk_falling_edges_t <- (df_ref_fall$i - first_zero_crossing_i) *
    1 /
    Fs +
    first_zero_cross_wfm_t
  message("determin tie, seperate rising edges and falling edges")

  tie_rise <- (df_wfm_rise$S - ref_clk_rising_edges_t)
  tie_fall <- (df_wfm_fall$S - ref_clk_falling_edges_t)
  if (is_rising_edge_leading) {
    tie_t <- interleave(ref_clk_rising_edges_t, ref_clk_falling_edges_t)
    tie_in_pS <- interleave(tie_rise, tie_fall) * 1e12
  } else {
    tie_t <- interleave(ref_clk_falling_edges_t, ref_clk_rising_edges_t)
    tie_in_pS <- interleave(tie_fall, tie_rise) * 1e12
  }
  df_tie <<- tibble(
    S = tie_t,
    pS = tie_in_pS
  )

  file_tie <- here(p_out, glue("{prefix}_tie.txt"))
  write_tsv(df_tie, file_tie)
  message("save to ", file_tie)

  message("interleave tie on rise edges and fall edges to one list")

  tie_jitter_pkpk_in_pS <- max(df_tie$pS) - min(df_tie$pS)
  message("TIE jitter PkPk (pS) = ", sprintf("%.3f", tie_jitter_pkpk_in_pS))

  t_max <- (OS * length(df_ref_raw_with_tags$V) * 1 / Fs) / 10
  dfp <- rbind(
    tibble(x = df_wfm$S, y = df_wfm$V, type = "wfm"),
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
    dplyr::filter(x <= t_max)
  p <- ggplot(dfp, aes(x, y, color = type, shape = type)) +
    geom_line(data = dfp[dfp$type == "wfm", ]) +
    geom_point() +
    geom_vline(xintercept = df_tie$S[df_tie$S <= t_max]) +
    scale_shape_manual(
      values = c(
        wfm = 4,
        `zero crossing - rising edge` = 4,
        `zero crossing - falling edge` = 3
      )
    ) +
    # scale_color_brewer(palette = "Dark2") +
    # scale_fill_brewer(palette = "Dark2") +
    scale_color_manual(values = c("#003049", "#d62828", "#f77f00")) +
    scale_fill_manual(values = c("#003049", "#d62828", "#f77f00")) +
    labs(title = "check zero crossing (zoom in)") +
    theme(legend.title = element_blank())
  ggp <- ggplotly(p)
  print(ggp)
  html_name <- glue("{prefix}_f7_gen_tie.html")
  wd <- getwd()
  setwd(here(p_fig))
  htmlwidgets::saveWidget(ggp, html_name)
  setwd(wd)
  message("save to ", html_name)
}

f8_gen_eye_db <- function() {
  format_message("f8: generate eye database")
  eye_samples_per_T <<- OS
  eye_samples_per_Window <<- (eye_window_length_l + eye_window_length_r) *
    eye_samples_per_T

  temp <- min(df_wfm_rise$i[1], df_wfm_fall$i[1]) -
    (eye_window_length_l - 0.5) * eye_samples_per_T
  eye_start_i <- ifelse(
    temp > 0,
    temp,
    (temp + ceiling(-temp / eye_samples_per_T) * eye_samples_per_T)
  )

  max_iter <- floor(
    (length(df_wfm$V) - eye_start_i + 1) / eye_samples_per_Window
  )
  eye_stop_i <- eye_start_i + max_iter * eye_samples_per_Window - 1

  eye_db <- list()
  eye_db_i <- list()
  eye_db_transition_tags <- c()
  for (k in 1:max_iter) {
    slice_index <- eye_start_i +
      (k - 1) * eye_samples_per_T +
      (1:eye_samples_per_Window - 1)
    eye_db_i[[k]] <- slice_index
    eye_db[[k]] <- df_wfm$V[slice_index]

    # only tag the 1st edge of the eye
    n <- floor(length(slice_index) / 2)
    tags <- names(df_ref$V[slice_index[1:n]])
    valid_tags <- tags[!is.na(tags)]
    if (length(valid_tags) > 0) {
      eye_db_transition_tags <- c(eye_db_transition_tags, valid_tags[1])
    } else {
      eye_db_transition_tags <- c(eye_db_transition_tags, NA)
    }
  }

  eye_db <- do.call(cbind, eye_db)
  colnames(eye_db) <- 1:ncol(eye_db)
  df_eye_db <<- as_tibble(eye_db)

  eye_db_i <- do.call(cbind, eye_db_i)
  colnames(eye_db_i) <- 1:ncol(eye_db_i)
  df_eye_db_i <<- as_tibble(eye_db_i)

  df_eye_db_transition_tags <<- tibble(tags = eye_db_transition_tags)

  file_eye_db <- here(p_out, glue("{prefix}_eye_db.txt"))
  write_tsv(df_eye_db, file_eye_db)
  message("save to ", file_eye_db)

  file_eye_db_i <- here(p_out, glue("{prefix}_eye_db_i.txt"))
  write_tsv(df_eye_db_i, file_eye_db_i)
  message("save to ", file_eye_db_i)

  file_eye_db_transition_tags <- here(p_out, glue("{prefix}_eye_db_tags.txt"))
  write_tsv(df_eye_db_transition_tags, file_eye_db_transition_tags)
  message("save to ", file_eye_db_transition_tags)

  message("update eye database")
  message("update eye index (sample index in wfm)")
  message("update eye transition tags")
}

f9_gen_eye_jitter <- function(tag = NA) {
  # use up to 1st crossing + 0.15 of eye database for 1st eye crossing, 2nd crossing - 0.15 to the end for 2nd eye crossing
  # return jitter in pS and datarate
  # and eye timing centered to zero crossing
  # tag: "1T - 1T"
  if (!is.null(df_eye_db_transition_tags) && !is.na(tag)) {
    format_message(str_c("f9: eye jitter for tag ", tag))

    tag_index <- which(df_eye_db_transition_tags$tags == tag)

    stopifnot(
      "!!! tag not found, expected format: 3T - 1T !!!" = length(tag_index) > 0
    )

    loc_eye_db <- df_eye_db[, tag_index]
  } else {
    format_message("f9: eye jitter for all transitions")
    loc_eye_db <- df_eye_db
  }

  n_row <- nrow(loc_eye_db)
  n_col <- ncol(loc_eye_db)
  n0 <- ceiling((eye_window_length_l - 0.15) * eye_samples_per_T)
  n1 <- floor((eye_window_length_l + 0.15) * eye_samples_per_T)
  loc_eye_t <- (0:(n_row - 1)) * 1 / Fs * 1e12 # pS
  x0 <- loc_eye_t[1:n0]
  x1 <- loc_eye_t[n1:length(loc_eye_t)]

  thr <- 0
  message("jitter at level ", thr, "V")

  jitter_t0 <- c()
  jitter_t1 <- c()

  for (j in 1:n_col) {
    # 1st eye crossing
    y0 <- loc_eye_db[[j]][1:n0]
    a <- y0[1]
    # b <- y0[floor(length(y0) / 2)]
    b <- y0[length(y0)]
    if ((a < b) && (a * b < 0)) {
      # rising edge
      t <- zero_crossing(
        x0,
        y0,
        thr = thr,
        find_rising_edge = TRUE,
        interp = TRUE
      )
      wi <- t$i
      wx <- t$x
      wy <- t$y
      jitter_t0 <- c(jitter_t0, wx)
    }

    if ((a > b) && (a * b < 0)) {
      # falling edge
      t <- zero_crossing(
        x0,
        y0,
        thr = thr,
        find_rising_edge = FALSE,
        interp = TRUE
      )
      wi <- t$i
      wx <- t$x
      wy <- t$y
      jitter_t0 <- c(jitter_t0, wx)
    }

    # 2nd eye crossing
    y1 <- loc_eye_db[[j]][n1:n_row]
    # a <- y1[ceiling(length(y1) / 2)]
    a <- y1[1]
    b <- y1[length(y1)]
    if ((a < b) && (a * b < 0)) {
      # rising edge
      t <- zero_crossing(
        x1,
        y1,
        thr = thr,
        find_rising_edge = TRUE,
        interp = TRUE
      )
      wi <- t$i
      wx <- t$x
      wy <- t$y
      jitter_t1 <- c(jitter_t1, wx)
    }

    if ((a > b) && (a * b < 0)) {
      # falling edge
      t <- zero_crossing(
        x1,
        y1,
        thr = thr,
        find_rising_edge = FALSE,
        interp = TRUE
      )
      wi <- t$i
      wx <- t$x
      wy <- t$y
      jitter_t1 <- c(jitter_t1, wx)
    }
  }

  crossing_t0 <- mean(jitter_t0)
  crossing_t1 <- mean(jitter_t1)
  message("first eye crossing at ", sprintf("%.3f", crossing_t0), "pS")
  message("second eye crossing at ", sprintf("%.3f", crossing_t1), "pS")

  df_eye_jitter <<- tibble(
    pS = jitter_t0 - crossing_t0
  )
  file_eye_jitter <- here(p_out, glue("{prefix}_eye_jitter.txt"))
  write_tsv(df_eye_jitter, file_eye_jitter)
  message("save to ", file_eye_jitter)

  eye_jitter_pkpk_in_pS <<- max(df_eye_jitter$pS) - min(df_eye_jitter$pS)
  if (is.na(tag)) {
    message("eye jitter PkPk ", sprintf("%.3f", eye_jitter_pkpk_in_pS), "pS")
  } else {
    message(
      "eye jitter PkPk ",
      sprintf("%.3f", eye_jitter_pkpk_pS_for_tag),
      "pS"
    )
    message("don't update eye width")
    message("don't update eye datarate")
  }
  df_eye_t <<- tibble(T = (loc_eye_t - mean(jitter_t0)) / (1 / F * 1e12))
  file_eye_t <- here(p_out, glue("{prefix}_eye_t.txt"))
  write_tsv(df_eye_t, file_eye_t)
  message("save to ", file_eye_t)

  eye_width_in_pS <- crossing_t1 - crossing_t0
  eye_datarate_Gbps <- 1 / (eye_width_in_pS * 1e-12) * 1e-9
  message("update eye width ", sprintf("%.3f", eye_width_in_pS), "pS")
  message("update eye datarate ", sprintf("%.3f", eye_datarate_Gbps), "Gbps")

  df_eye_property <- tibble(
    Fs = Fs,
    F = F,
    OS = OS,
    prbs_length = prbs_length,
    `cross0(pS)` = crossing_t0,
    `cross1(pS)` = crossing_t1,
    `width(pS)` = eye_width_in_pS,
    `datarate(Gbps)` = eye_datarate_Gbps
  )
  file_eye_property <- here(p_out, glue("{prefix}_eye_property.txt"))
  write_tsv(df_eye_property, file_eye_property)
  message("save to ", file_eye_property)
}

f1_load_raw_wfm()
f1_load_raw_ref()
f1_view_raw_wfm()
f1_view_raw_ref()
f2_upsample()
f3_detect_pattern()
f4_roll_ref()
f5_remove_leading_trailing_samples()
f6_find_zero_crossing()
f7_gen_tie()
f8_gen_eye_db()
f9_gen_eye_jitter()
