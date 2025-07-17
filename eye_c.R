# parse wfm file saved from scope
# and reference file
# and genrate eye diagram, jitter
#
# this file is deprecated
# refer eye_parser.R for updated code
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
# library(ggplotify)
# library(cowplot)

Eye <- R6Class(
  "Eye",
  private = list(
    os = NA,
    prbs_length = NA,
    ref_raw_with_tags = NA,
    linear_interp = function(x0, y0, x1, y1, y) {
      return(x0 + (y - y0) * (x1 - x0) / (y1 - y0))
    },
    format_message = function(s) {
      message("")
      message(str_flatten(rep("-", str_count(s))))
      message(s)
      message(str_flatten(rep("-", str_count(s))))
      message("")
    },
    zero_crossing = function(x, y, thr = 0, find_rising_edge = T, interp = F) {
      if (find_rising_edge) {
        zc_i <- seq_along(y)[
          (y[2:length(y)] >= thr) & (y[1:(length(y) - 1)] < thr)
        ]
      } else {
        zc_i <- seq_along(y)[
          (y[2:length(y)] <= thr) & (y[1:(length(y) - 1)] > thr)
        ]
      }

      if (interp) {
        zc_x <- c()
        for (i in zc_i) {
          zc_x <- c(
            zc_x,
            private$linear_interp(x[i], y[i], x[i + 1], y[i + 1], thr)
          )
        }
        zc_y <- rep(thr, length(zc_i))
      } else {
        zc_x <- x[zc_i]
        zc_y <- y[zc_i]
      }
      return(tibble(i = zc_i, x = zc_x, y = zc_y))
    },
    zero_crossing_wfm = function(
      ref_i,
      wfm_t,
      wfm_v,
      find_rising_edge = T,
      interp = T
    ) {
      d <- 10 # find within a small range
      zc_i <- c()
      zc_x <- c()
      zc_y <- c()

      for (i in ref_i) {
        rng <- seq(max(i - d, 0), min(i + d, length(wfm_v)))
        t <- private$zero_crossing(
          wfm_t[rng],
          wfm_v[rng],
          0,
          find_rising_edge,
          interp
        )
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
    },
    interleave = function(l1, l2) {
      n1 <- length(l1)
      n2 <- length(l2)

      if (n1 != n2) message("length not equal")

      l3 <- c()
      for (i in 1:max(n1, n2)) {
        if (i <= n1) l3 <- c(l3, l1[i])
        if (i <= n2) l3 <- c(l3, l2[i])
      }
      return(l3)
    }
  ),
  public = list(
    p_in = NULL,
    p_out = NULL,
    p_fig = NULL,
    wfm_file = NULL,
    ref_file = NULL,
    Fs = NA,
    F = NA,
    eye_window_length_l = NULL,
    eye_window_length_r = NULL,
    eye_alpha = NULL,
    df_wfm_raw = NULL,
    df_ref_raw = NULL,
    wfm_v = NULL,
    wfm_t = NULL,
    ref_t = NULL,
    ref_v = NULL,
    tie_t = NULL,
    tie_in_pS = NULL,
    tie_jitter_pkpk_in_pS = NULL,
    eye_start_i = NULL,
    eye_stop_i = NULL,
    eye_samples_per_T = NULL,
    eye_samples_per_Window = NULL,
    eye_t = NULL,
    eye_db = NULL,
    eye_db_i = NULL,
    eye_db_transition_tags = NULL,
    eye_datarate_Gbps = NULL,
    eye_jitter_in_pS_for_tag = NULL,
    eye_jitter_pkpk_in_pS_for_tag = NULL,
    eye_jitter_in_pS = NULL, # zero crossing
    eye_jitter_pkpk_in_pS = NULL, # zero crossing
    eye_width_in_pS = NULL,
    wfm_rise_i = NULL,
    wfm_rise_t = NULL,
    wfm_rise_v = NULL,
    wfm_fall_i = NULL,
    wfm_fall_t = NULL,
    wfm_fall_v = NULL,
    ref_rise_i = NULL,
    ref_rise_t = NULL,
    ref_rise_v = NULL,
    ref_fall_i = NULL,
    ref_fall_t = NULL,
    ref_fall_v = NULL,
    ref_clk_rising_edges_t = NULL,
    ref_clk_falling_edges_t = NULL,
    initialize = function(
      p_in = "01_data/",
      p_out = "02_results/",
      p_fig = "06_fig/",
      wfm_file = "3G_OSA25_WFM.dat",
      ref_file = "prbs7.txt",
      Fs = 40e9,
      F = 3e9,
      eye_window_length_l = 0.8, # first 0.8T in eye window
      eye_window_length_r = 0.8, # second 0.8T in eye window
      eye_alpha = 0.3
    ) {
      self$p_in <- p_in
      self$p_out <- p_out
      self$p_fig <- p_fig
      self$wfm_file <- wfm_file
      self$ref_file <- ref_file
      self$Fs <- Fs
      self$F <- F
      private$os <- self$Fs / self$F
      self$eye_window_length_l <- eye_window_length_l
      self$eye_window_length_r <- eye_window_length_r
      self$eye_alpha <- eye_alpha
      self$print()
    },
    print = function(...) {
      cat("Eye: (require Fs and F)\n")
      cat("  p_in:              ", self$p_in, "\n", sep = "")
      cat("  p_out:             ", self$p_out, "\n", sep = "")
      cat("  wfm_file:          ", self$wfm_file, "\n", sep = "")
      cat("  ref_file:          ", self$ref_file, "\n", sep = "")
      cat("  Fs(S/s):           ", self$Fs, "\n", sep = "")
      cat("  F(bps):            ", self$F, "\n", sep = "")
      cat("  os:                ", private$os, "\n", sep = "")
      cat("  eye_window_length_l: ", self$eye_window_length_l, "\n", sep = "")
      cat("  eye_window_length_r: ", self$eye_window_length_r, "\n", sep = "")
      cat("  eye_alpha:         ", self$eye_alpha, "\n", sep = "")
      invisible(self)
    },
    get_os = function() private$os,
    get_prbs_length = function() private$prbs_length,
    get_ref_raw_with_tags = function() private$ref_raw_with_tags,
    f1_read_files = function() {
      private$format_message("f1: load wfm and pat files")
      self$df_wfm_raw <- read_delim(
        str_c(self$p_in, self$wfm_file),
        col_names = c("s", "v"),
        col_types = list(col_double(), col_double()),
        delim = " "
      )
      self$df_ref_raw <- read_delim(
        str_c(self$p_in, self$ref_file),
        col_names = c("ch", "mk1", "mk2"),
        col_types = list(col_double(), col_integer(), col_integer()),
        delim = ","
      )
      message("  read wfm file: ", str_c(self$p_in, self$wfm_file))
      message("  read ref file: ", str_c(self$p_in, self$ref_file))
    },
    f1_view_raw_wfm_and_ref = function() {
      p_wfm <- ggplot(
        self$df_wfm_raw[1:(nrow(self$df_ref_raw) * private$os * 3), ],
        aes(s, v)
      ) +
        geom_line() +
        # geom_point() +
        labs(title = "wfm") +
        theme_classic()
      print(ggplotly(p_wfm))

      p_ref <- ggplot(
        self$df_ref_raw %>% mutate(`#` = seq_along(self$df_ref_raw$mk1)),
        aes(`#`, mk1)
      ) +
        geom_step() +
        # geom_point() +
        labs(title = "ref") +
        theme_classic()
      print(ggplotly(p_ref))
    },
    f2_upsample = function(plot = F) {
      private$format_message("f2: upsample wfm")
      if ((self$Fs %% self$F) > 0) {
        message(
          "  F = ",
          self$F,
          " Fs = ",
          self$Fs,
          " fractional sampling, do interpolation"
        )

        wfm_v_old <- self$df_wfm_raw$v
        wfm_t_old <- (seq_along(wfm_v_old) - 1) * (1 / self$Fs)
        Fs_old <- self$Fs

        private$os <- ceiling(Fs_old / self$F / 10) * 10 # new os
        self$Fs <- self$F * private$os # new Fs
        up <- self$Fs / Fs_old

        message(
          "  F = ",
          self$F,
          " Fs = ",
          self$Fs,
          "up by ",
          up,
          "x, new os ",
          private$os
        )

        n <- length(wfm_v_old)
        n_new <- as.integer(n * up)

        x <- pracma::linspace(0, 1, n)
        x_new <- pracma::linspace(0, 1, n_new)

        f <- splinefun(x, wfm_v_old)
        self$wfm_v <- f(x_new)
        self$wfm_t <- (seq_along(self$wfm_v) - 1) * (1 / self$Fs)

        cnt <- 500
        dfp <- rbind(
          tibble(s = wfm_t_old[1:cnt], v = wfm_v_old[1:cnt], type = "original"),
          tibble(
            s = self$wfm_t[1:(cnt * up)],
            v = self$wfm_v[1:(cnt * up)],
            type = "upsampled"
          )
        )

        if (plot) {
          p <- ggplot(
            dfp,
            aes(s, v, color = type, shape = type)
          ) +
            # geom_line() +
            geom_point() +
            scale_shape_manual(values = c(3, 4)) +
            # scale_shape_manual(values = c(20, 3)) +
            scale_color_brewer(palette = "Dark2") +
            scale_fill_brewer(palette = "Dark2") +
            labs(title = "upsample wfm") +
            theme_classic() +
            theme(legend.title = element_blank())
          ggp <- ggplotly(p)
          print(ggp)
          html_name <- str_c(self$p_fig, "f2_upsample.html")
          htmlwidgets::saveWidget(ggp, html_name)
          message("  save to ", html_name)
        }
      } else {
        message("  integer sampling, no interpolation")

        private$os <- self$Fs / self$F
        self$wfm_v <- self$df_wfm_raw$v
        self$wfm_t <- (seq_along(self$wfm_v) - 1) * (1 / self$Fs)
      }
    },
    f3_detect_pattern = function(plot = F) {
      private$format_message(
        "f3: detect prbs length and update transition tags"
      )
      ref_one_cycle <- self$df_ref_raw$mk1

      c <- ccf(
        ref_one_cycle,
        ref_one_cycle,
        plot = T,
        lag.max = (length(ref_one_cycle) - 1)
      )
      c_max <- max(c$acf)
      c_max_lag <- c$lag[which(c$acf == c_max)] # it is 0
      high_corr_lags <- c$lag[which(c$acf > 0.2)]
      private$prbs_length <- high_corr_lags[which(high_corr_lags > c_max_lag)][
        1
      ]
      message("  prbs length ", private$prbs_length)

      ref_patterns <- rep(ref_one_cycle[1:private$prbs_length], 3)

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
        pats[1:private$prbs_length],
        tags[(private$prbs_length + 1):(private$prbs_length * 2)]
      )
      private$ref_raw_with_tags <- rep(
        ref_one_prbs_with_tags,
        ceiling(length(ref_one_cycle) / private$prbs_length)
      )[1:length(ref_one_cycle)]
      message("  update raw ref with tags of nT transition")

      if (plot) {
        dfp <- rbind(
          tibble(x = seq_along(pats), y = pats, type = "ref"),
          tibble(x = zc, y = 0, type = "zc")
        )
        p <- ggplot(dfp, aes(x, y, color = type)) +
          geom_step() +
          geom_point() +
          labs(title = "ref raw and zero crossings") +
          theme_classic() +
          theme(legend.title = element_blank())
        ggp <- ggplotly(p)
        print(ggp)
        html_name <- str_c(self$p_fig, "f3_detect_pattern.html")
        htmlwidgets::saveWidget(ggp, html_name)
        message("  save to ", html_name)
      }
    },
    f4_roll_ref = function(plot = F) {
      private$format_message("f4: roll pat to match wfm")
      # ref_one_cycle <- rep(self$df_ref$ch1 * 2 - 1, each = private$os)
      tags <- names(private$ref_raw_with_tags)
      index_with_tags <- seq_along(private$ref_raw_with_tags)[!is.na(tags)]

      ref_one_cycle <- rep(private$ref_raw_with_tags, each = private$os)
      names(ref_one_cycle) <- NA
      names(ref_one_cycle)[index_with_tags * private$os] <- tags[!is.na(tags)]
      message("  update upsampled ref with tags of nT transition")

      wfm_one_cycle <- self$wfm_v[1:length(ref_one_cycle)]
      ref_one_cycle_invert <- -1 * ref_one_cycle

      # gsignal::xcorr does the same thing as python np.correlate(wfm_one_cycle, ref_one_cycle, mode="full")
      # c0 <- gsignal::xcorr(wfm_one_cycle, ref_one_cycle)
      # c0_max <- max(c0$R)
      # c0_max_lag <- c0$lags[which(c0$R == max(c0$R))]
      # c1 <- gsignal::xcorr(wfm_one_cycle, ref_one_cycle_invert)
      # c1_max <- max(c1$R)
      # c1_max_lag <- c1$lags[which(c1$R == max(c1$R))]
      c <- ccf(
        wfm_one_cycle,
        ref_one_cycle,
        lag.max = (length(wfm_one_cycle) - 1),
        plot = F
      )
      c_max <- max(c$acf)
      c_min <- min(c$acf)
      c_max_lag <- c$lag[which(c$acf == max(c$acf))]
      c_min_lag <- c$lag[which(c$acf == min(c$acf))]

      if (c_max > abs(c_min)) {
        i <- c_max_lag
        # ref_one_cycle_rolled <- c(ref_one_cycle[(length(ref_one_cycle) - i + 1):length(ref_one_cycle)], ref_one_cycle[1:(length(ref_one_cycle) - i)])
        ref_one_cycle_rolled <- c(
          tail(ref_one_cycle, i),
          head(ref_one_cycle, -i)
        )
        message(
          "  ref samples correlate wfm, roll ref ",
          i,
          " to align with wfm"
        )
      } else {
        i <- c_min_lag
        ref_one_cycle_rolled <- c(
          tail(rev(ref_one_cycle_invert), i),
          head(rev(ref_one_cycle_invert), -i)
        )
        message(
          "  inverted ref samples correlate wfm, roll ref ",
          i,
          " to align with wfm"
        )
      }

      self$ref_v <- rep(
        ref_one_cycle_rolled,
        ceiling(length(self$wfm_v) / length(ref_one_cycle_rolled))
      )[1:length(self$wfm_v)]
      self$ref_t <- (seq_along(self$ref_v) - 1) * (1 / self$Fs)

      if (plot) {
        dfp <- rbind(
          tibble(
            `#` = seq_along(wfm_one_cycle),
            y = wfm_one_cycle,
            type = "wfm"
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
          scale_color_brewer(palette = "Dark2") +
          scale_fill_brewer(palette = "Dark2") +
          labs(title = "one cycle wfm and rolled ref") +
          theme_classic() +
          theme(legend.title = element_blank())
        ggp <- ggplotly(p)
        print(ggp)
        html_name <- str_c(self$p_fig, "f4_roll_ref.html")
        htmlwidgets::saveWidget(ggp, html_name)
        message("  save to ", html_name)
      }
    },
    f5_find_zero_crossing = function(plot = F) {
      private$format_message("f5: find zero crossing")
      t <- private$zero_crossing(
        self$ref_t,
        self$ref_v,
        thr = 0,
        find_rising_edge = T,
        interp = F
      )
      self$ref_rise_i <- t$i
      self$ref_rise_t <- t$x
      self$ref_rise_v <- t$y
      t <- private$zero_crossing(
        self$ref_t,
        self$ref_v,
        thr = 0,
        find_rising_edge = F,
        interp = F
      )
      self$ref_fall_i <- t$i
      self$ref_fall_t <- t$x
      self$ref_fall_v <- t$y
      message("  find ref zero crossing, seperate rising and falling edges")

      t <- private$zero_crossing_wfm(
        self$ref_rise_i,
        self$wfm_t,
        self$wfm_v,
        find_rising_edge = T,
        interp = T
      )
      self$wfm_rise_i <- t$i
      self$wfm_rise_t <- t$x
      self$wfm_rise_v <- t$y
      t <- private$zero_crossing_wfm(
        self$ref_fall_i,
        self$wfm_t,
        self$wfm_v,
        find_rising_edge = F,
        interp = T
      )
      self$wfm_fall_i <- t$i
      self$wfm_fall_t <- t$x
      self$wfm_fall_v <- t$y
      message(
        "  find wfm zero crossing, seperate rising and falling edges, refer ref crossing edges"
      )
    },
    f6_gen_tie = function(plot = F) {
      private$format_message("f6: generate tie")
      is_rising_edge_leading <- (self$ref_rise_i[1] <= self$ref_fall_i[1])
      first_zero_crossing_i <- ifelse(
        is_rising_edge_leading,
        self$ref_rise_i[1],
        self$ref_fall_i[1]
      )
      first_zero_cross_wfm_t <- ifelse(
        is_rising_edge_leading,
        self$wfm_rise_t[1],
        self$wfm_fall_t[1]
      )
      message("  determin first wfm zero crossing location and time")

      self$ref_clk_rising_edges_t <- (self$ref_rise_i - first_zero_crossing_i) *
        1 /
        self$Fs +
        first_zero_cross_wfm_t
      self$ref_clk_falling_edges_t <- (self$ref_fall_i -
        first_zero_crossing_i) *
        1 /
        self$Fs +
        first_zero_cross_wfm_t
      message("  determin tie, seperate rising edges and falling edges")

      tie_rise <- (self$wfm_rise_t - self$ref_clk_rising_edges_t)
      tie_fall <- (self$wfm_fall_t - self$ref_clk_falling_edges_t)
      # if (is_rising_edge_leading) {
      #   self$tie_t <- unlist(mapply(list, self$ref_clk_rising_edges_t, self$ref_clk_falling_edges_t, SIMPLIFY = F))
      #   self$tie_in_pS <- unlist(mapply(list, tie_rise, tie_fall, SIMPLIFY = F)) * 1e12
      # } else {
      #   self$tie_t <- unlist(mapply(list, self$ref_clk_falling_edges_t, self$ref_clk_rising_edges_t, SIMPLIFY = F))
      #   self$tie_in_pS <- unlist(mapply(list, tie_fall, tie_rise, SIMPLIFY = F)) * 1e12
      # }
      if (is_rising_edge_leading) {
        # self$tie_t <- unlist(mapply(list, self$ref_clk_rising_edges_t, self$ref_clk_falling_edges_t, SIMPLIFY = F))
        # self$tie_in_pS <- unlist(mapply(list, tie_rise, tie_fall, SIMPLIFY = F)) * 1e12
        self$tie_t <- private$interleave(
          self$ref_clk_rising_edges_t,
          self$ref_clk_falling_edges_t
        )
        self$tie_in_pS <- private$interleave(tie_rise, tie_fall) * 1e12
      } else {
        # self$tie_t <- unlist(mapply(list, self$ref_clk_falling_edges_t, self$ref_clk_rising_edges_t, SIMPLIFY = F))
        # self$tie_in_pS <- unlist(mapply(list, tie_fall, tie_rise, SIMPLIFY = F)) * 1e12
        self$tie_t <- private$interleave(
          self$ref_clk_falling_edges_t,
          self$ref_clk_rising_edges_t
        )
        self$tie_in_pS <- private$interleave(tie_fall, tie_rise) * 1e12
      }
      message("  interleave tie on rise edges and fall edges to one list")

      self$tie_jitter_pkpk_in_pS <- max(self$tie_in_pS) - min(self$tie_in_pS)
      message(
        "  TIE jitter PkPk (pS) = ",
        sprintf("%.3f", self$tie_jitter_pkpk_in_pS)
      )

      t_max <- (private$os * length(private$ref_raw_with_tags) * 1 / self$Fs) /
        10
      if (plot) {
        dfp <- rbind(
          tibble(x = self$wfm_t, y = self$wfm_v, type = "wfm"),
          tibble(
            x = self$wfm_rise_t,
            y = self$wfm_rise_v,
            type = "zero crossing - rising edge"
          ),
          tibble(
            x = self$wfm_fall_t,
            y = self$wfm_fall_v,
            type = "zero crossing - falling edge"
          )
        ) %>%
          dplyr::filter(x <= t_max)
        p <- ggplot(dfp, aes(x, y, color = type, shape = type)) +
          geom_line(data = dfp[dfp$type == "wfm", ]) +
          geom_point() +
          geom_vline(xintercept = self$tie_t[self$tie_t <= t_max]) +
          scale_shape_manual(
            values = c(
              wfm = 4,
              `zero crossing - rising edge` = 4,
              `zero crossing - falling edge` = 3
            )
          ) +
          scale_color_brewer(palette = "Dark2") +
          scale_fill_brewer(palette = "Dark2") +
          labs(title = "check zero crossing (zoom in)") +
          theme_classic() +
          theme(legend.title = element_blank())
        ggp <- ggplotly(p)
        print(ggp)
        html_name <- str_c(self$p_fig, "f6_gen_tie.html")
        htmlwidgets::saveWidget(ggp, html_name)
        message("  save to ", html_name)
      }
    },
    f7_gen_eye_db = function() {
      private$format_message("f7: generate eye database")
      self$eye_samples_per_T <- private$os
      self$eye_samples_per_Window <- (self$eye_window_length_l +
        self$eye_window_length_r) *
        self$eye_samples_per_T

      temp <- min(self$wfm_rise_i[1], self$wfm_fall_i[1]) -
        (self$eye_window_length_l - 0.5) * self$eye_samples_per_T
      self$eye_start_i <- ifelse(
        temp > 0,
        temp,
        (temp +
          ceiling(-temp / self$eye_samples_per_T) * self$eye_samples_per_T)
      )

      max_iter <- floor(
        (length(self$wfm_v) - self$eye_start_i + 1) /
          self$eye_samples_per_Window
      )
      self$eye_stop_i <- self$eye_start_i +
        max_iter * self$eye_samples_per_Window -
        1

      eye_db <- list()
      eye_db_i <- list()
      self$eye_db_transition_tags <- c()
      for (k in 1:max_iter) {
        slice_index <- self$eye_start_i +
          (k - 1) * self$eye_samples_per_T +
          (1:self$eye_samples_per_Window - 1)
        eye_db_i[[k]] <- slice_index
        eye_db[[k]] <- self$wfm_v[slice_index]

        # only tag the 1st edge of the eye
        n <- floor(length(slice_index) / 2)
        tags <- names(self$ref_v[slice_index[1:n]])
        valid_tags <- tags[!is.na(tags)]
        if (length(valid_tags) > 0) {
          self$eye_db_transition_tags <- c(
            self$eye_db_transition_tags,
            valid_tags[1]
          )
        } else {
          self$eye_db_transition_tags <- c(self$eye_db_transition_tags, NA)
        }
      }
      self$eye_db <- do.call(cbind, eye_db)
      self$eye_db_i <- do.call(cbind, eye_db_i)
      message("  update eye database")
      message("  update eye index (sample index in wfm)")
      message("  update eye transition tags")
    },
    f8_gen_eye_jitter = function(tag = NA) {
      # use up to 1st crossing + 0.15 of eye database for 1st eye crossing, 2nd crossing - 0.15 to the end for 2nd eye crossing
      # return jitter in pS and datarate
      # and eye timing centered to zero crossing
      # tag: "1T - 1T"
      if (!is.null(self$eye_db_transition_tags) && !is.na(tag)) {
        private$format_message(str_c("f8: eye jitter for tag ", tag))

        tag_index <- which(self$eye_db_transition_tags == tag)

        stopifnot(
          "  !!! tag not found, expected format: 3T - 1T !!!" = length(
            tag_index
          ) >
            0
        )

        loc_eye_db <- self$eye_db[, tag_index]
      } else {
        private$format_message("f8: eye jitter for all transitions")
        loc_eye_db <- self$eye_db
      }

      n_row <- nrow(loc_eye_db)
      n_col <- ncol(loc_eye_db)
      n0 <- ceiling((self$eye_window_length_l - 0.15) * self$eye_samples_per_T)
      n1 <- floor((self$eye_window_length_l + 0.15) * self$eye_samples_per_T)
      loc_eye_t <- (0:(n_row - 1)) * 1 / self$Fs * 1e12 # pS
      x0 <- loc_eye_t[1:n0]
      x1 <- loc_eye_t[n1:length(loc_eye_t)]

      thr <- 0
      message("  jitter at level ", thr, "V")

      jitter_t0 <- c()
      jitter_t1 <- c()

      for (j in 1:n_col) {
        # 1st eye crossing
        y0 <- loc_eye_db[1:n0, j]
        a <- y0[1]
        # b <- y0[floor(length(y0) / 2)]
        b <- y0[length(y0)]
        if ((a < b) && (a * b < 0)) {
          # rising edge
          t <- private$zero_crossing(
            x0,
            y0,
            thr = thr,
            find_rising_edge = T,
            interp = T
          )
          wi <- t$i
          wx <- t$x
          wy <- t$y
          jitter_t0 <- c(jitter_t0, wx)
        }

        if ((a > b) && (a * b < 0)) {
          # falling edge
          t <- private$zero_crossing(
            x0,
            y0,
            thr = thr,
            find_rising_edge = F,
            interp = T
          )
          wi <- t$i
          wx <- t$x
          wy <- t$y
          jitter_t0 <- c(jitter_t0, wx)
        }

        # 2nd eye crossing
        y1 <- loc_eye_db[n1:n_row, j]
        # a <- y1[ceiling(length(y1) / 2)]
        a <- y1[1]
        b <- y1[length(y1)]
        if ((a < b) && (a * b < 0)) {
          # rising edge
          t <- private$zero_crossing(
            x1,
            y1,
            thr = thr,
            find_rising_edge = T,
            interp = T
          )
          wi <- t$i
          wx <- t$x
          wy <- t$y
          jitter_t1 <- c(jitter_t1, wx)
        }

        if ((a > b) && (a * b < 0)) {
          # falling edge
          t <- private$zero_crossing(
            x1,
            y1,
            thr = thr,
            find_rising_edge = F,
            interp = T
          )
          wi <- t$i
          wx <- t$x
          wy <- t$y
          jitter_t1 <- c(jitter_t1, wx)
        }
      }

      crossing_t0 <- mean(jitter_t0)
      crossing_t1 <- mean(jitter_t1)
      message("  first eye crossing at ", sprintf("%.3f", crossing_t0), "pS")
      message("  second eye crossing at ", sprintf("%.3f", crossing_t1), "pS")

      if (is.na(tag)) {
        self$eye_jitter_in_pS <- jitter_t0 - crossing_t0
        self$eye_jitter_pkpk_in_pS <- max(jitter_t0) - min(jitter_t0)
        message(
          "  eye jitter PkPk ",
          sprintf("%.3f", self$eye_jitter_pkpk_in_pS),
          "pS"
        )
      } else {
        self$eye_jitter_in_pS_for_tag <- jitter_t0 - crossing_t0
        self$eye_jitter_pkpk_in_pS_for_tag <- max(jitter_t0) - min(jitter_t0)
        message(
          "  eye jitter PkPk ",
          sprintf("%.3f", self$eye_jitter_pkpk_pS_for_tag),
          "pS"
        )
        message("  don't update eye width")
        message("  don't update eye datarate")
      }
      self$eye_t <- (loc_eye_t - mean(jitter_t0)) / (1 / self$F * 1e12)
      self$eye_width_in_pS <- crossing_t1 - crossing_t0
      self$eye_datarate_Gbps <- 1 / (self$eye_width_in_pS * 1e-12) * 1e-9
      message(
        "  update eye width ",
        sprintf("%.3f", self$eye_width_in_pS),
        "pS"
      )
      message(
        "  update eye datarate ",
        sprintf("%.3f", self$eye_datarate_Gbps),
        "Gbps"
      )
    },
    f9_view_sliced_wfm = function(tag = NA, cnt = 4) {
      private$format_message("f9: view wfm and sliced wfm for eye")
      if (is.na(tag)) {
        n <- min(cnt, ncol(self$eye_db))
        x_start_min <- min(self$eye_db_i[, 1:n])
        x_start_max <- max(self$eye_db_i[, 1:n])
        x_end_min <- min(self$eye_db_i[,
          ((ncol(self$eye_db_i) - n + 1):ncol(self$eye_db_i))
        ])
        x_end_max <- max(self$eye_db_i[,
          ((ncol(self$eye_db_i) - n + 1):ncol(self$eye_db_i))
        ])
        t_start <- tibble(
          `#` = seq(x_start_min, x_start_max),
          v = self$wfm_v[x_start_min:x_start_max],
          type = "wfm"
        )
        t_start <- rbind(
          t_start,
          tibble(
            `#` = seq(x_start_min, x_start_max),
            v = self$ref_v[x_start_min:x_start_max],
            type = "ref"
          )
        )
        levels <- c("wfm", "ref")
        for (i in 1:n) {
          t_start <- rbind(
            t_start,
            tibble(
              `#` = self$eye_db_i[, i],
              v = self$eye_db[, i],
              type = str_c("#", i)
            )
          )
          levels <- c(levels, str_c("#", i))
        }
        t_start$type <- factor(t_start$type, levels = levels)
        t_end <- tibble(
          `#` = seq(x_end_min, x_end_max),
          v = self$wfm_v[x_end_min:x_end_max],
          type = "wfm"
        )
        t_end <- rbind(
          t_end,
          tibble(
            `#` = seq(x_end_min, x_end_max),
            v = self$ref_v[x_end_min:x_end_max],
            type = "ref"
          )
        )
        levels <- c("wfm", "ref")
        for (i in n:1) {
          t_end <- rbind(
            t_end,
            tibble(
              `#` = self$eye_db_i[, (ncol(self$eye_db_i) - i + 1)],
              v = self$eye_db[, (ncol(self$eye_db_i) - i + 1)],
              type = str_c("#-", abs(-i + 1))
            )
          )
          levels <- c(levels, str_c("#-", abs(-i + 1)))
        }
        t_end$type <- factor(t_end$type, levels = levels)
      } else {
        tag_index <- which(self$eye_db_transition_tags == tag)

        stopifnot(
          "  !!! tag not found, expected format: 3T - 1T !!!" = length(
            tag_index
          ) >
            0
        )

        loc_eye_db <- self$eye_db[, tag_index]
        loc_eye_db_i <- self$eye_db_i[, tag_index]

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
          v = self$wfm_v[x_start_min:x_start_max],
          type = "wfm"
        )
        t_start <- rbind(
          t_start,
          tibble(
            `#` = seq(x_start_min, x_start_max),
            v = self$ref_v[x_start_min:x_start_max],
            type = "ref"
          )
        )
        levels <- c("wfm", "ref")
        for (i in 1:n) {
          t_start <- rbind(
            t_start,
            tibble(
              `#` = loc_eye_db_i[, i],
              v = loc_eye_db[, i],
              type = str_c("#", i)
            )
          )
          levels <- c(levels, str_c("#", i))
        }
        t_start$type <- factor(t_start$type, levels = levels)
        t_end <- tibble(
          `#` = seq(x_end_min, x_end_max),
          v = self$wfm_v[x_end_min:x_end_max],
          type = "wfm"
        )
        t_end <- rbind(
          t_end,
          tibble(
            `#` = seq(x_end_min, x_end_max),
            v = self$ref_v[x_end_min:x_end_max],
            type = "ref"
          )
        )
        levels <- c("wfm", "ref")
        for (i in n:1) {
          t_end <- rbind(
            t_end,
            tibble(
              `#` = loc_eye_db_i[, (ncol(loc_eye_db_i) - i + 1)],
              v = loc_eye_db[, (ncol(loc_eye_db_i) - i + 1)],
              type = str_c("#-", abs(-i + 1))
            )
          )
          levels <- c(levels, str_c("#-", abs(-i + 1)))
        }
        t_end$type <- factor(t_end$type, levels = levels)
      }
      p_start <- ggplot(
        t_start,
        aes(x = `#`, y = v, color = type, shape = type)
      ) +
        geom_line(data = t_start[!t_start$type == "ref", ]) +
        geom_step(data = t_start[t_start$type == "ref", ]) +
        # geom_point() +
        facet_grid(type ~ ., scale = "free") +
        labs(title = "start") +
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") +
        # theme_classic() +
        theme(
          strip.text.y = element_blank(),
          legend.title = element_blank()
        )
      # print(p_start)
      ggp_start <- ggplotly(p_start)
      print(ggp_start)
      html_name <- str_c(self$p_fig, "f9_view_sliced_wfm_start.html")
      htmlwidgets::saveWidget(ggp_start, html_name)
      message("  save to ", html_name)

      p_end <- ggplot(t_end, aes(x = `#`, y = v, color = type, shape = type)) +
        geom_line(data = t_end[!t_end$type == "ref", ]) +
        geom_step(data = t_end[t_end$type == "ref", ]) +
        # geom_point() +
        facet_grid(type ~ ., scale = "free") +
        labs(title = "end") +
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") +
        # theme_classic() +
        theme(
          strip.text.y = element_blank(),
          legend.title = element_blank()
        )
      # print(p_end)
      ggp_end <- ggplotly(p_end)
      print(ggp_end)
      html_name <- str_c(self$p_fig, "f9_view_sliced_wfm_end.html")
      htmlwidgets::saveWidget(ggp_end, html_name)
      message("  save to ", html_name)
    },
    f9_ggplot_zoomin_wfm_ref_tie = function(
      max_lines = NULL,
      color = "steelblue",
      bins = 50
    ) {
      t_max <- (private$os * length(private$ref_raw_with_tags) * 1 / self$Fs) /
        10 *
        1e9 # nS

      data_wfm <- rbind(
        tibble(x = self$wfm_t, y = self$wfm_v, type = "wfm (zoom in)"),
        tibble(
          x = self$wfm_rise_t,
          y = self$wfm_rise_v,
          type = "zero crossing - rising edge"
        ),
        tibble(
          x = self$wfm_fall_t,
          y = self$wfm_fall_v,
          type = "zero crossing - falling edge"
        )
      ) %>%
        mutate(x = x * 1e9) %>%
        dplyr::filter(x <= t_max)

      vlines <- self$tie_t * 1e9
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
            length(self$wfm_v),
            "points, ",
            max(self$wfm_t) * 1e9,
            "nS"
          ),
          x = "nS",
          y = "V"
        ) +
        theme_bw() +
        theme(legend.title = element_blank(), legend.position = "top")

      # data_ref <- tibble(x = (seq_along(self$ref_v) - 1) * (1 / self$Fs) * 1e12, y = self$ref_v, type = "genie (zoom in)")
      data_ref <- tibble(
        x = self$ref_t * 1e9,
        y = self$ref_v,
        type = "genie (zoom in)"
      )
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

      data_tie <- tibble(
        x = self$tie_t * 1e9,
        y = self$tie_in_pS,
        type = "tie (zoom in)"
      )
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
          x = self$tie_t * 1e9,
          y = self$tie_in_pS,
          type = "tie (zoom in)"
        )
      } else {
        tie_index <- which(
          self$tie_t <= (tail(self$eye_db_i[, max_lines], 1) - 1) * 1 / self$Fs
        )
        data_tie_for_hist <- tibble(
          x = self$tie_t[tie_index] * 1e9,
          y = self$tie_in_pS[tie_index],
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
    },
    f9_ggplot_eye = function(
      loc_eye_db = NULL,
      index = NULL,
      max_lines = 500,
      color = "steelblue",
      alpha = 0.2,
      bins = 50
    ) {
      # loc_eye_db is the matrix
      if (is.null(loc_eye_db)) loc_eye_db <- self$eye_db
      if (is.null(index)) {
        if (is.null(self$eye_t)) {
          index <- 1:nrow(loc_eye_db)
        } else {
          index <- self$eye_t
        }
      }

      colnames(loc_eye_db) <- 1:ncol(loc_eye_db)
      t <- as_tibble(loc_eye_db, .name_repair = "unique") %>%
        mutate(index = index)
      if (is.null(max_lines)) {
        max_lines <- ncol(loc_eye_db)
      } else {
        max_lines <- min(max_lines, ncol(loc_eye_db))
      }

      p <- ggplot(t, aes(x = index))
      for (k in 1:max_lines) {
        p <- p +
          geom_line(aes(y = .data[[names(t)[k]]]), color = color, alpha = alpha)
      }
      p <- p +
        labs(
          title = str_c(
            sprintf("%.3f", self$eye_datarate_Gbps),
            "Gbps, ",
            max_lines,
            "UI"
          ),
          x = "T",
          y = "V"
        ) +
        theme_bw()
      return(p)
    },
    f9_ggplot_eye_and_wfm = function(
      loc_eye_db = NULL,
      index = NULL,
      max_lines = 500,
      color = "steelblue",
      alpha = 0.2,
      bins = 50
    ) {
      l <- self$f9_ggplot_zoomin_wfm_ref_tie(max_lines, color, bins)

      p_eye <- self$f9_ggplot_eye(
        loc_eye_db,
        index,
        max_lines,
        color,
        alpha,
        bins
      )

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
      fig_name <- str_c(self$p_fig, "f9_ggplot_eye.png")
      ggsave(fig_name, p, width = 10, height = 6)
      message("  save to ", fig_name)
      return(p)
    },
    f9_ggmatplot_eye = function(
      loc_eye_db = NULL,
      index = NULL,
      max_lines = 500,
      color = "steelblue",
      alpha = 0.2,
      bins = 50
    ) {
      # loc_eye_db is the matrix
      if (is.null(loc_eye_db)) loc_eye_db <- self$eye_db
      if (is.null(index)) {
        if (is.null(self$eye_t)) {
          index <- 1:nrow(loc_eye_db)
        } else {
          index <- self$eye_t
        }
      }
      if (is.null(max_lines)) {
        max_lines <- ncol(loc_eye_db)
      } else {
        max_lines <- min(max_lines, ncol(loc_eye_db))
      }

      p <- ggmatplot(
        index,
        eye$eye_db[, 1:max_lines],
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
            sprintf("%.3f", self$eye_datarate_Gbps),
            "Gbps, ",
            max_lines,
            "UI"
          )
          # ylab = "V"
        ) +
        theme_bw()
      return(p)
    },
    f9_ggmatplot_eye_and_wfm = function(
      loc_eye_db = NULL,
      index = NULL,
      max_lines = 500,
      color = "steelblue",
      alpha = 0.2,
      bins = 50
    ) {
      l <- self$f9_ggplot_zoomin_wfm_ref_tie(max_lines, color, bins)

      p_eye <- self$f9_ggmatplot_eye(
        loc_eye_db,
        index,
        max_lines,
        color,
        alpha,
        bins
      )

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
      fig_name <- str_c(self$p_fig, "f9_ggmatplot_eye.png")
      ggsave(fig_name, p, width = 10, height = 6)
      message("  save to ", fig_name)
      return(p)
    },
    f9_matplot_eye = function(
      loc_eye_db = NULL,
      index = NULL,
      max_lines = NULL,
      color = "steelblue",
      alpha = 0.2,
      bins = 50
    ) {
      # loc_eye_db is the matrix
      if (is.null(loc_eye_db)) loc_eye_db <- self$eye_db
      if (is.null(index)) {
        if (is.null(self$eye_t)) {
          index <- 1:nrow(loc_eye_db)
        } else {
          index <- self$eye_t
        }
      }
      if (is.null(max_lines)) {
        max_lines <- ncol(loc_eye_db)
      } else {
        max_lines <- min(max_lines, ncol(loc_eye_db))
      }
      matplot(
        index,
        self$eye_db[, 1:max_lines],
        self$eye_db[, 1:max_lines],
        type = "l",
        lty = 1,
        col = alpha(color, alpha),
        xlab = "T",
        ylab = "V",
        main = str_c(
          sprintf("%.3f", self$eye_datarate_Gbps),
          "Gbps, ",
          max_lines,
          "UI"
        )
      )
    },
    f9_matplot_eye_and_wfm = function(
      loc_eye_db = NULL,
      index = NULL,
      max_lines = NULL,
      color = "steelblue",
      alpha = 0.2,
      bins = 50
    ) {
      l <- self$f9_ggplot_zoomin_wfm_ref_tie(max_lines, color, bins)

      # png(str_c(self$p_fig, "f9_matplot_eye.png"), width = 12, height = 6, units = "in", res = 300)

      layout(
        matrix(c(2, 2, 1, 3, 3, 1, 4, 4, 5), nrow = 3, ncol = 3, byrow = T),
        heights = c(2, 2, 2),
        widths = c(2, 2, 2)
      )

      self$f9_matplot_eye(loc_eye_db, index, max_lines, color, alpha, bins)

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
  )
)

eye <- Eye$new()
eye$f1_read_files()
# eye$f1_view_raw_wfm_and_ref()
eye$f2_upsample(plot = F)
eye$f3_detect_pattern(plot = F)
eye$f4_roll_ref(plot = F)
eye$f5_find_zero_crossing(plot = F)
eye$f6_gen_tie(plot = F)
eye$f7_gen_eye_db()
eye$f8_gen_eye_jitter()
eye$f9_view_sliced_wfm(tag = NA)
# p <- eye$f9_ggplot_eye_and_wfm(max_lines = 200)
# p <- eye$f9_ggmatplot_eye_and_wfm(max_lines = 100)
# eye$f9_matplot_eye_and_wfm(max_lines = 100)
eye$f9_matplot_eye_and_wfm(max_lines = NULL)
dev.print(pdf, file = str_c(eye$p_fig, "f9_matplot_eye.pdf"))
