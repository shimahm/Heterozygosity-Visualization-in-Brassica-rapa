#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(forcats)
  library(ggplot2)
  library(ggideogram)
  library(ggnewscale)  # allow multiple fill scales
  library(grid)
  library(scales)
})

# ----------------------------- USER SETTINGS -----------------------------
infile         <- "filter_uniq_paired_turnip15.txt"  # your file
out_prefix     <- "heterozygosity_turnip15_ggideo"

# Filters & windowing
win_size_bp    <- 250000
min_QUAL       <- 30
min_DP         <- 10

# Define "het-like" band (expected ~0.5 if truly heterozygous)
het_low        <- 0.30
het_high       <- 0.70

# Alarm when % het-like sites in a window exceeds this (tune 0.01–0.10)
alarm_pct_het  <- 0.05   # 5%

# Right-edge color bar max (cap for contrast). Increase if needed.
pct_het_max    <- 0.25   # 25%

# Per-chromosome point cap (0 = plot all het-like points)
max_points_per_chr <- 20000

# Output sizes
pdf_w <- 7.5; pdf_h <- 10
png_w <- 2100; png_h <- 3000; png_dpi <- 180

# ------------------------ CHROM SIZES & CENTROMERES ----------------------
chrom_lengths <- tibble::tribble(
  ~Chrom, ~Len,
  "A01", 43052898, "A02", 37474953, "A03", 42291893, "A04", 29842065, "A05", 46744448,
  "A06", 52009500, "A07", 38256961, "A08", 31123574, "A09", 73365914, "A10", 30432017
)

centromeres <- tibble::tribble(
  ~Chrom, ~CenStart, ~CenEnd,
  "A01", 19110000, 26170000, "A02", 19100000, 26200000, "A03", 33800000, 38240000,
  "A04",  6710000, 13800000, "A05", 27280000, 30800000, "A06", 16100000, 20100000,
  "A07",  5950000, 16700000, "A08",  6500000, 14680000, "A09", 44250000, 49350000,
  "A10",  5550000, 17300000
)

chrom_order <- paste0("A", sprintf("%02d", 1:10))

# ------------------------------- LOAD DATA -------------------------------
message("Reading: ", infile)
dt <- fread(infile, header = TRUE, data.table = FALSE)

# normalize expected column names
names(dt) <- gsub("^chrom$", "CHROM", names(dt), ignore.case = TRUE)
names(dt) <- gsub("^pos$",   "POS",   names(dt), ignore.case = TRUE)
names(dt) <- gsub("^refdp$", "refDP", names(dt), ignore.case = TRUE)
names(dt) <- gsub("^altdp$", "altDP", names(dt), ignore.case = TRUE)
names(dt) <- gsub("^qual$",  "QUAL",  names(dt), ignore.case = TRUE)
names(dt) <- gsub("^dp$|^totaldp$|^total_depth$|^depth$|^total$", "total", names(dt), ignore.case = TRUE)

req <- c("CHROM","POS","refDP","altDP","QUAL")
miss <- setdiff(req, names(dt)); if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse=", "))

dt <- dt %>%
  filter(CHROM %in% chrom_order) %>%
  mutate(
    POS   = as.numeric(POS),
    QUAL  = as.numeric(QUAL),
    refDP = as.numeric(refDP),
    altDP = as.numeric(altDP),
    total = if ("total" %in% names(dt)) as.numeric(total) else (refDP + altDP),
    AB    = ifelse(refDP + altDP > 0, altDP / (refDP + altDP), NA_real_)
  ) %>%
  filter(!is.na(AB), QUAL >= min_QUAL, total >= min_DP)

if (nrow(dt) == 0) stop("No rows after filtering. Relax thresholds or check columns.")

# ----------------------- IDEOGRAM SCAFFOLD (p/q arms) --------------------
k <- chrom_lengths %>%
  left_join(centromeres, by = "Chrom") %>%
  mutate(
    Chrom = fct_relevel(Chrom, chrom_order),
    Mid   = floor((CenStart + CenEnd)/2L)
  ) %>%
  rowwise() %>%
  mutate(bands = list(tibble::tibble(
    Chrom = Chrom,
    Start = c(0L,      CenStart, Mid,      CenEnd),
    End   = c(CenStart, Mid,      CenEnd,  Len),
    Stain = c("gneg",  "acen",    "acen",  "gneg"),
    Arm   = c("p",     "p",       "q",     "q")
  ))) %>%
  ungroup() %>%
  tidyr::unnest(bands, names_sep = ".") %>%
  dplyr::select(
    Chrom, Len, Mid,
    Start = bands.Start,
    End   = bands.End,
    Stain = bands.Stain,
    Arm   = bands.Arm
  ) %>%
  mutate(
    Chrom = fct_relevel(Chrom, chrom_order),
    y0 = Start - Mid, y1 = End - Mid
  ) %>%
  filter(y0 < y1)

mid_lookup <- k %>% distinct(Chrom, Mid)
len_lookup <- chrom_lengths %>% rename(ChrLen = Len)
x_lookup   <- tibble::tibble(Chrom = levels(k$Chrom)) %>% mutate(x = row_number())

# --------------------- MAP VARIANTS TO IDEOGRAM COORDS -------------------
dt <- dt %>%
  mutate(Chrom = factor(CHROM, levels = chrom_order)) %>%
  left_join(mid_lookup, by = "Chrom") %>%
  left_join(len_lookup, by = "Chrom") %>%
  left_join(x_lookup,  by = "Chrom") %>%
  mutate(y = POS - Mid)

# Window assignment (clamped to chr end)
dt <- dt %>%
  mutate(
    win_start = ((POS - 1L) %/% win_size_bp) * win_size_bp + 1L,
    win_end   = pmin(win_start + win_size_bp - 1L, ChrLen)
  )

# ------------------- HET-FOCUSED STATS & POINTS --------------------------
# Flag het-like sites
dt <- dt %>%
  mutate(is_het = AB >= het_low & AB <= het_high)

# Window stats: % het-like
wstats <- dt %>%
  group_by(Chrom, x, Mid, win_start, win_end) %>%
  summarise(
    n_sites = n(),
    pct_het = mean(is_het),
    .groups = "drop"
  ) %>%
  mutate(
    y0 = win_start - Mid,
    y1 = win_end   - Mid,
    alarm = (pct_het >= alarm_pct_het)
  )

# Points: only het-like sites to reduce clutter
dt_het <- dt %>% filter(is_het)

# Optional per-chrom downsampling of het points (robust; no n() pitfalls)
if (max_points_per_chr > 0) {
  dt_het <- dt_het %>%
    group_by(Chrom) %>%
    mutate(.row_in_grp = row_number()) %>%
    ungroup()

  keep_idx <- dt_het %>%
    group_by(Chrom) %>%
    summarise(
      keep = list({
        n_g <- n()
        if (n_g <= max_points_per_chr) seq_len(n_g) else sample.int(n_g, max_points_per_chr)
      }),
      .groups = "drop"
    ) %>%
    tidyr::unnest_longer(keep, values_to = ".row_in_grp")

  dt_het <- dt_het %>%
    inner_join(keep_idx, by = c("Chrom", ".row_in_grp")) %>%
    select(-.row_in_grp)
}

# ------------------------------ PLOTTING ---------------------------------
cyto_cols <- c(
  gneg="#FFFFFF", gpos25="#C8C8C8", gpos50="#A0A0A0",
  gpos75="#787878", gpos100="#000000", gvar="#E0E0E0",
  stalk="#708090", acen="#2ca02c"
)

chr_width <- 0.80
off_left  <- -chr_width/2 + 0.03    # het point ticks hug left edge
off_right <-  chr_width/2 - 0.03    # %het track hugs right edge
track_half_width <- 0.06

# Gradient for % het-like (0 → pct_het_max)
het_pal <- seq_gradient_pal("#edf8fb", "#2171b5", "Lab")  # light -> dark blue

p <- ggplot(k) +
  # ideogram with discrete fill (stains)
  geom_ideogram(
    aes(x = Chrom, ymin = y0, ymax = y1, chrom = Chrom, fill = Stain, arm = Arm),
    radius = grid::unit(7, "pt"), width = chr_width, linewidth = 0.25, colour = "black",
    show.legend = FALSE
  ) +
  scale_fill_manual(values = cyto_cols) +
  ggnewscale::new_scale_fill() +  # start a new fill scale for %het track

  # Alarm windows (light grey across full ideogram width)
  geom_rect(
    data = wstats %>% filter(alarm),
    aes(xmin = as.numeric(x) - chr_width/2, xmax = as.numeric(x) + chr_width/2,
        ymin = y0, ymax = y1),
    inherit.aes = FALSE, fill = "#D9D9D980", colour = NA
  ) +

  # Right-edge track colored by % het-like in window
  geom_rect(
    data = wstats,
    aes(xmin = as.numeric(x) + off_right - track_half_width,
        xmax = as.numeric(x) + off_right + track_half_width,
        ymin = y0, ymax = y1, fill = pmin(pct_het, pct_het_max)),
    inherit.aes = FALSE, colour = NA, alpha = 0.95, show.legend = TRUE
  ) +
  scale_fill_gradientn(
    colours = het_pal(seq(0,1,length.out=7)),
    limits = c(0, pct_het_max),
    name = "% het-like",
    labels = percent_format(accuracy = 1)
  ) +

  # Left-edge points: het-like sites only (solid, small)
  geom_point(
    data = dt_het,
    aes(x = as.numeric(x) + off_left, y = y),
    inherit.aes = FALSE, size = 0.5, stroke = 0.2, alpha = 0.8, shape = 16
  ) +

  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  scale_y_continuous(
    name = "Position (Mb, centromere aligned at 0)",
    labels = function(z) z/1e6,
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    x = "Chromosome",
    title = "A-genome ideogram highlighting heterozygous-like regions",
    subtitle = paste0(
      "Windows: ", format(win_size_bp, big.mark=","), " bp; Filters: QUAL ≥ ", min_QUAL,
      ", DP ≥ ", min_DP, "; Alarm: % het-like ≥ ", percent(alarm_pct_het)
    ),
    caption = "Green = centromere (acen); grey = alarm windows; right bar = % het-like per window"
  ) +
  theme(
    panel.background = element_blank(),
    axis.title.x     = element_blank(),
    axis.text.x      = element_text(),
    axis.ticks.x     = element_blank(),
    axis.title.y     = element_text(),
    legend.position  = "right"
  )

# ------------------------------ SAVE -------------------------------------
ggsave(paste0(out_prefix, ".pdf"), p, width = pdf_w, height = pdf_h, units = "in")
ggsave(paste0(out_prefix, ".png"), p, width = pdf_w, height = pdf_h, units = "in", dpi = png_dpi)

message("Done: ", out_prefix, ".pdf and ", out_prefix, ".png")

