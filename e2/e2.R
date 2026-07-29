## clear workspace
rm(list = ls())

options(download.file.method = "libcurl")

if (!require(pacman)) { install.packages("pacman") }
pacman::p_load('ggplot2',
               'effsize',
               'tidyr',
               'dplyr',
               'lmerTest',
               'MatchIt',
               'cobalt',
               'interactions',
               'exact2x2',
               'reshape2',
               'ltm',
               'lavaan',
               'pwr',
               'BayesFactor',
               'patchwork')

# Helper: extract BF10 from a BayesFactor object and classify it against the
# predefined thresholds (Jeffreys): BF10 > 3 = moderate+ evidence for H1,
# BF10 < 1/3 = moderate+ evidence for H0, otherwise inconclusive.
bf_report <- function(x, y, paired = FALSE, label = "") {
    bf <- ttestBF(x, y, paired = paired)
    bf10 <- as.numeric(extractBF(bf)$bf)
    verdict <- if (bf10 > 3) "moderate+ evidence for H1" else
               if (bf10 < 1 / 3) "moderate+ evidence for H0" else "inconclusive"
    cat(paste0("    BF10", if (label != "") paste0(" [", label, "]") else "",
               " = ", signif(bf10, 3), " (", verdict, ")\n"))
    invisible(bf10)
}

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory
}

# ---- Nature Research Journals figure settings (per NRJ "Guide to preparing final artwork") ----
# - Sans-serif typeface (Helvetica/Arial); all text between 5 pt and 7 pt.
# - Single-column = 88 mm, double-column = 180 mm; figures saved as vector PDF in RGB.
# The longitudinal plots below overlay individual participant observations on the condition means.
NATURE_FONT    <- "Helvetica"
NATURE_1COL_MM <- 88
NATURE_2COL_MM <- 180

DAYS <- 1:10

# ---- Modal (majority) condition per participant ----------------------------------
# Condition is recorded fresh in each day's record. A couple of participants closed
# and re-entered the Day-1 survey, producing an inconsistent Day-1 condition tag
# while every subsequent day agreed. We treat each participant's modal (most
# frequent) recorded condition as authoritative and assign it consistently across
# all of their days, so no participant is split across conditions.
.cond_all <- do.call(rbind, lapply(DAYS, function(DAY) {
    x <- read.csv(paste0("d", DAY, ".csv"))
    x <- x[x$Finished == "1" & !is.na(x$worker_id) & x$worker_id != "" &
           !is.na(x$condition) & x$condition != "", c("worker_id", "condition")]
    x
}))
modal_condition <- tapply(.cond_all$condition, .cond_all$worker_id,
                          function(z) names(which.max(table(z))))
apply_modal_condition <- function(df) {
    hit <- df$worker_id %in% names(modal_condition)
    df$condition[hit] <- modal_condition[df$worker_id[hit]]
    df
}

############################ READ DAY 1 & APPLY EXCLUSIONS ############################

d <- read.csv("d1.csv")
d <- d[d$Finished == "1", ]
d <- d[!is.na(d$worker_id) & d$worker_id != "", ]
d$day <- 1

size_before <- dim(d)[1]
d <- d[(d$att1 == "2" & d$att2 == "2"), ]
print(paste0("Number of participants hired (before loneliness screening): ", dim(d)[1]))

d_screened_out <- d[d$screenedOutLowLoneliness == "true", ]
d <- d[d$screenedOutLowLoneliness != "true", ]
print(paste0("Number of participants hired (after loneliness screening): ", dim(d)[1]))

# Exclude participants based on comprehension check (all 3 conditions)
size_before <- dim(d)[1]
d <- d[d$comp_1 == "1", ]
print(paste0("Exclusions from comprehension check: ", size_before - dim(d)[1]))
print(paste0("Number of participants after exclusions: ", dim(d)[1]))

# Assign the modal condition (see top of script) so enrollment counts and all
# downstream analyses use one consistent condition per participant.
d <- apply_modal_condition(d)

sample_before <- table(d$condition)
print(sample_before)

############################ AGGREGATE RATINGS FROM ALL DAYS ############################
GET_ONLY_FINISHED <- TRUE  # If TRUE, only include participants who finished the study (all 10 days)

for (DAY in DAYS) {
    if (DAY == 1) { next }  # We already read day 1

    d_day <- read.csv(paste0("d", DAY, ".csv"))
    d_day <- d_day[d_day$Finished == "1", ]
    d_day <- d_day[!is.na(d_day$worker_id) & d_day$worker_id != "", ]
    d_day$day <- DAY

    # Remove duplicate worker_ids within the same day. Some participants have a
    # blank "shell" submission (Finished == 1, Progress == 100, but no condition
    # and no answers) alongside their real response. Keep-first would sometimes
    # retain the shell and discard the real submission, so we instead keep the
    # most complete row per worker: prefer a non-blank condition, then more
    # non-missing loneliness items.
    if (anyDuplicated(d_day$worker_id)) {
        cond_ok <- as.integer(!(is.na(d_day$condition) | d_day$condition == ""))
        lon_cols <- intersect(c("lonely_before_1", "lonely_before_2", "lonely_before_3",
                                 "lonely_after_1",  "lonely_after_2",  "lonely_after_3"),
                               names(d_day))
        lon_ok <- rowSums(!is.na(d_day[, lon_cols, drop = FALSE]) &
                          d_day[, lon_cols, drop = FALSE] != "")
        score <- cond_ok * 100 + lon_ok
        d_day <- d_day[order(d_day$worker_id, -score), ]
        d_day <- d_day[!duplicated(d_day$worker_id), ]
    }

    # Drop any remaining blank "shell" rows (no condition assigned): these are
    # not genuine completions of the day's content.
    d_day <- d_day[!(is.na(d_day$condition) | d_day$condition == ""), ]

    print(paste0("Number of participants on Day ", DAY, ": ", dim(d_day)[1]))

    # Bind with running dataframe, ignoring missing columns
    d <- dplyr::bind_rows(d, d_day)

    if (DAY == DAYS[length(DAYS)]) {  # Last day
        if (GET_ONLY_FINISHED) {
            ids_per_day <- lapply(DAYS, function(x) d[d$day == x, 'worker_id'])
            common_ids <- Reduce(intersect, ids_per_day)
            d <- d[d$worker_id %in% common_ids, ]
        }
        d <- d[!is.na(d$worker_id), ]
    }
}

# Apply the modal condition across all retained rows so every participant has a
# single, consistent condition in the analyses below (shell detection above used
# the raw per-day condition and is unaffected).
d <- apply_modal_condition(d)

print(paste0("Number of participants who finished the study: ", length(unique(d$worker_id))))
print(table(d[, c('condition', 'day')]))

# Attrition rates by condition
print("Attrition rates:")
for (cnd in c("both_text", "both_voice", "control")) {
    attr_cnd <- sample_before[cnd] - table(d[d$condition == cnd, 'day'])[1]
    print(paste0(cnd, ": ", attr_cnd, " of ", sample_before[cnd]))
}

# Test whether attrition differs across conditions (condition x completed/dropped).
conds3 <- c("both_text", "both_voice", "control")
completed_ct <- sapply(conds3, function(cnd) sum(d$condition == cnd & d$day == 1))
dropped_ct   <- as.numeric(sample_before[conds3]) - completed_ct
attr_tab <- rbind(completed = completed_ct, dropped = dropped_ct)
colnames(attr_tab) <- conds3
cat("\n--- Attrition: condition x completed/dropped ---\n"); print(attr_tab)
print(chisq.test(attr_tab))
attr_pair <- function(a, b) {
    pt <- prop.test(c(dropped_ct[a], dropped_ct[b]),
                    c(as.numeric(sample_before[a]), as.numeric(sample_before[b])))
    cat(sprintf("Attrition pairwise %s vs %s: chi2(1) = %.2f, p = %.3f\n",
                a, b, pt$statistic, pt$p.value))
}
cat("--- Attrition pairwise comparisons ---\n")
attr_pair("both_voice", "control")
attr_pair("both_text",  "control")
attr_pair("both_voice", "both_text")

print(paste0("Age: ", mean(as.numeric(d[d$day == 1, ]$age), trim = 0, na.rm = TRUE)))
print(paste0("Percentage of females: ", sum(d[d$day == 1, 'gender'] == "2")[1] / dim(d[d$day == 1, ])[1]))

table(d$condition, d$day)

############## PREPARING THE DATASET ##############

# Convert relevant columns to numeric
cols <- colnames(d)
cols <- cols[grepl("lone", cols) | grepl("comp", cols) | grepl("preference", cols) |
             grepl("ai_capability", cols) | grepl("ai_companion_exp", cols)]
d[cols] <- lapply(d[cols], as.numeric)

# Aggregate the 3 loneliness items (mean) for before/after measurements
d$lonely_before <- rowMeans(d[, c("lonely_before_1", "lonely_before_2", "lonely_before_3")], na.rm = TRUE)
d$lonely_after  <- rowMeans(d[, c("lonely_after_1",  "lonely_after_2",  "lonely_after_3")],  na.rm = TRUE)

# Cronbach's alpha for the loneliness items, per condition
for (cnd in c("both_text", "both_voice", "control")) {
    l_b <- d[d$condition == cnd, c('lonely_before_1', 'lonely_before_2', 'lonely_before_3')]
    l_b <- mutate_all(l_b, function(x) as.numeric(as.character(x)))
    cat("Cronbach alpha (lonely_before),", cnd, ":\n")
    print(cronbach.alpha(l_b, na.rm = TRUE))

    if (cnd != "control") {
        l_a <- d[d$condition == cnd, c('lonely_after_1', 'lonely_after_2', 'lonely_after_3')]
        l_a <- mutate_all(l_a, function(x) as.numeric(as.character(x)))
        cat("Cronbach alpha (lonely_after),", cnd, ":\n")
        print(cronbach.alpha(l_a, na.rm = TRUE))
    }
}

################## FLATTEN THE DATASET ##################

# Keep relevant columns and flatten to long format
d_plot <- d[, c('worker_id', 'condition', 'day',
                'lonely_before', 'lonely_after',
                'preference_before_1', 'preference_after_1')]

d_long <- d_plot %>%
    pivot_longer(cols = c('lonely_before', 'lonely_after'),
                 names_to = "timepoint",
                 values_to = "loneliness") %>%
    mutate(timepoint = ifelse(grepl("before", timepoint), "before", "after"))

# Drop NA loneliness rows. Control has no 'after' measurement, so those drop out here.
d_long <- d_long[!is.na(d_long$loneliness) & !is.nan(d_long$loneliness), ]

# Preference: use post-interaction for text/voice; for control use the single measurement (preference_before_1).
d_long$preference <- ifelse(d_long$condition == "control",
                            d_long$preference_before_1,
                            ifelse(d_long$timepoint == "after",
                                   d_long$preference_after_1,
                                   d_long$preference_before_1))

############ PLOT LONELINESS SCORES OVER 10 DAYS ############

# Human-readable labels for legends
condition_labels <- c("control" = "Control", "both_text" = "Text", "both_voice" = "Voice")
timepoint_labels <- c("before" = "Before interaction", "after"  = "After interaction")

d_long$condition_label <- factor(condition_labels[d_long$condition],
                                 levels = c("Control", "Text", "Voice"))
d_long$timepoint_label <- factor(timepoint_labels[d_long$timepoint],
                                 levels = c("Before interaction", "After interaction"))

dir.create("plots", showWarnings = FALSE)

# Axis-break glyph: two short parallel diagonal strokes drawn on the y-axis near its
# base to signal that the y-axis is truncated (does not start at zero). `xpos` is the
# x position of the y-axis line (the panel's expanded x-range minimum, obtained from
# ggplot_build); the strokes start on the axis and rise to the right, into the panel.
axis_break <- function(xpos, ybot, dx, h, gap, lwd = 0.3) {
    list(
        # white mask between the two strokes to "break" (erase) the axis line there
        annotate("polygon",
                 x = c(xpos - dx, xpos + dx, xpos + dx, xpos - dx),
                 y = c(ybot,      ybot + h,  ybot + gap + h, ybot + gap),
                 fill = "white", colour = NA),
        annotate("segment", x = xpos - dx, xend = xpos + dx, y = ybot,       yend = ybot + h,
                 linewidth = lwd, lineend = "round", colour = "black"),
        annotate("segment", x = xpos - dx, xend = xpos + dx, y = ybot + gap, yend = ybot + gap + h,
                 linewidth = lwd, lineend = "round", colour = "black")
    )
}

# Combined factor used to drive shape so each condition × timepoint gets a distinct
# pch (open for "before", filled for "after"). The shape legend itself is hidden;
# Condition and Timepoint each have their own legend with override.aes.
d_long$cond_time <- interaction(d_long$condition_label, d_long$timepoint_label,
                                sep = " - ", drop = TRUE)

p_lonely <- ggplot(data = d_long, aes(x = day, y = loneliness,
                          color = condition_label,
                          fill = condition_label,
                          linetype = timepoint_label,
                          shape = cond_time,
                          group = interaction(condition_label, timepoint_label))) +
    # Note: individual points are not overlaid here — with 3 conditions x 2 timepoints
    # the raw points cannot cleanly distinguish before/after; the 95% CI ribbons show spread.
    stat_summary(fun = mean, geom = "line", linewidth = 0.5) +
    stat_summary(fun = mean, geom = "point", size = 1.6, stroke = 0.6) +
    stat_summary(fun.data = mean_cl_normal, fun.args = list(conf.int = 0.95),
                 geom = "ribbon", alpha = 0.10, colour = NA, show.legend = FALSE) +
    # Okabe-Ito colorblind-safe palette: gray, blue, vermillion
    scale_color_manual(name = NULL,
                       values = c("Control" = "#555555",
                                  "Text"    = "#0072B2",
                                  "Voice"   = "#D55E00")) +
    scale_fill_manual(name = NULL,
                      values = c("Control" = "#555555",
                                 "Text"    = "#0072B2",
                                 "Voice"   = "#D55E00")) +
    # Shape encodes condition × timepoint:
    # Text = circle, Voice = square, Control = triangle
    # Before = open (pch 1/0/2), After = filled (pch 16/15/17)
    scale_shape_manual(values = c(
        "Control - Before interaction" = 2,  "Control - After interaction" = 17,
        "Text - Before interaction"    = 1,  "Text - After interaction"    = 16,
        "Voice - Before interaction"   = 0,  "Voice - After interaction"   = 15
    )) +
    scale_linetype_manual(name = NULL,
                          values = c("Before interaction" = "23",
                                     "After interaction"  = "solid")) +
    scale_x_continuous(breaks = 1:10) +
    scale_y_continuous(breaks = seq(30, 80, by = 10)) +
    # Pin the x display range so the axis position is fixed and the break glyph
    # (drawn left of the axis, below) cannot shift it.
    coord_cartesian(xlim = c(1, 10), ylim = c(20, 80), clip = "off") +
    theme_classic() +
    labs(x = "Day", y = "Loneliness Score") +
    guides(
        # Hide raw shape legend (covered by Condition + Timepoint legends below).
        shape = "none",
        # Condition legend: filled shape per condition for clarity.
        color = guide_legend(order = 1, title.position = "top", title.hjust = 0.5,
                             override.aes = list(linetype = 0,
                                                 shape = c(17, 16, 15),
                                                 size = 1.6)),
        fill  = guide_legend(order = 1, title.position = "top", title.hjust = 0.5,
                             override.aes = list(linetype = 0,
                                                 shape = c(17, 16, 15),
                                                 size = 1.6)),
        # Timepoint legend: open vs filled circle to mirror what's in the plot.
        linetype = guide_legend(order = 2,
                                override.aes = list(shape = NA,
                                                    color = "black",
                                                    linewidth = 0.6))
    ) +
    theme(text = element_text(size = 7, family = NATURE_FONT),
          axis.text = element_text(size = 6, colour = "black"),
          axis.line = element_line(linewidth = 0.25),
          axis.ticks = element_line(linewidth = 0.25),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "top",
          legend.box = "horizontal",
          legend.box.spacing = unit(0.1, "cm"),
          legend.box.margin = margin(0, 0, 2, 0),
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.key = element_blank(),
          legend.key.height = unit(0.8, "lines"),
          legend.key.width = unit(1.8, "lines"),
          legend.spacing.x = unit(0.3, "cm"),
          legend.margin = margin(2, 6, 2, 6))

# (saved below as panel a of the combined two-panel Figure 2)

############ PLOT PREFERENCE (AI vs. HUMANS) OVER 10 DAYS ############

# For text/voice: preference_before_1 (before interaction) and preference_after_1
# (after interaction). For control: a single rating (preference_before_1) per day,
# carried on the "before" row of d_long.
d_pref_long <- d_long[!is.na(d_long$preference), ]

# Zoom: data spans ~25-50 on the 0-100 slider; tighten y-range so lines are readable.
# Endpoints of the slider: 0 = prefer Humans, 100 = prefer AI.
pref_min <- 20
pref_max <- 55

p_pref <- ggplot(data = d_pref_long, aes(x = day, y = preference,
                               color = condition_label,
                               fill = condition_label,
                               linetype = timepoint_label,
                               shape = cond_time,
                               group = interaction(condition_label, timepoint_label))) +
    # Note: individual points are not overlaid here — with 3 conditions x 2 timepoints
    # the raw points cannot cleanly distinguish before/after; the 95% CI ribbons show spread.
    stat_summary(fun = mean, geom = "line", linewidth = 0.5) +
    stat_summary(fun = mean, geom = "point", size = 1.6, stroke = 0.6) +
    stat_summary(fun.data = mean_cl_normal, fun.args = list(conf.int = 0.95),
                 geom = "ribbon", alpha = 0.10, colour = NA, show.legend = FALSE) +
    scale_color_manual(name = NULL,
                       values = c("Control" = "#555555",
                                  "Text"    = "#0072B2",
                                  "Voice"   = "#D55E00")) +
    scale_fill_manual(name = NULL,
                      values = c("Control" = "#555555",
                                 "Text"    = "#0072B2",
                                 "Voice"   = "#D55E00")) +
    scale_shape_manual(values = c(
        "Control - Before interaction" = 2,  "Control - After interaction" = 17,
        "Text - Before interaction"    = 1,  "Text - After interaction"    = 16,
        "Voice - Before interaction"   = 0,  "Voice - After interaction"   = 15
    )) +
    scale_linetype_manual(name = NULL,
                          values = c("Before interaction" = "23",
                                     "After interaction"  = "solid")) +
    scale_x_continuous(breaks = 1:10, expand = expansion(mult = c(0.04, 0.02))) +
    scale_y_continuous(breaks = seq(pref_min + 5, pref_max, by = 5)) +
    coord_cartesian(xlim = c(1, 10), ylim = c(pref_min, pref_max), clip = "off") +
    theme_classic() +
    labs(x = "Day", y = "Preference for AI vs. Humans (0 = Humans, 100 = AI)") +
    guides(
        shape = "none",
        color = guide_legend(order = 1, title.position = "top", title.hjust = 0.5,
                             override.aes = list(linetype = 0,
                                                 shape = c(17, 16, 15),
                                                 size = 1.6)),
        fill  = guide_legend(order = 1, title.position = "top", title.hjust = 0.5,
                             override.aes = list(linetype = 0,
                                                 shape = c(17, 16, 15),
                                                 size = 1.6)),
        linetype = guide_legend(order = 2,
                                override.aes = list(shape = NA,
                                                    color = "black",
                                                    linewidth = 0.6))
    ) +
    theme(text = element_text(size = 7, family = NATURE_FONT),
          axis.text = element_text(size = 6, colour = "black"),
          axis.line = element_line(linewidth = 0.25),
          axis.ticks = element_line(linewidth = 0.25),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          legend.position = "top",
          legend.box = "horizontal",
          legend.box.spacing = unit(0.1, "cm"),
          legend.box.margin = margin(0, 0, 2, 0),
          legend.title = element_blank(),
          legend.text = element_text(size = 6),
          legend.key = element_blank(),
          legend.key.height = unit(0.8, "lines"),
          legend.key.width = unit(1.8, "lines"),
          legend.spacing.x = unit(0.3, "cm"),
          legend.margin = margin(2, 6, 2, 6))

# (saved below as panel b of the combined two-panel Figure 2)

############ COMBINED FIGURE 2: (a) loneliness and (b) preference over 10 days ############
# Two-panel figure following the Nature Research Journals "Guide to preparing final
# artwork": double-column width (180 mm), Helvetica, 5-7 pt text, vector PDF in RGB,
# panels labelled with bold lower-case letters. Panel a (loneliness) on top, panel b
# (preference) on the bottom; the shared Condition/Timepoint legend is collected once.
# Add the y-axis break exactly on each axis line. The axis line sits at the panel's
# expanded x-range minimum, which we read back from the built plot so the glyph lands
# on the axis regardless of the scale expansion.
xr_l <- ggplot_build(p_lonely)$layout$panel_params[[1]]$x.range
p_lonely <- p_lonely + axis_break(xpos = xr_l[1], ybot = 22.5, dx = 0.11, h = 1.7, gap = 1.1)
xr_p <- ggplot_build(p_pref)$layout$panel_params[[1]]$x.range
p_pref  <- p_pref  + axis_break(xpos = xr_p[1], ybot = pref_min + 2.0, dx = 0.11, h = 1.1, gap = 0.7)

fig2 <- (p_lonely / p_pref) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "top",
          plot.tag = element_text(size = 8, face = "bold", family = NATURE_FONT))

ggsave("plots/figure2_loneliness_and_preference.pdf", fig2,
       width = NATURE_2COL_MM, height = 190, units = "mm")

# ---- Slides version: identical figure with larger text all around ----
# Reuses the same built panels (p_lonely / p_pref already carry the axis-break
# glyphs) and only enlarges text via a theme applied across both panels with `&`.
# Saved to a separate file; the Nature-spec figure above is left unchanged.
slide_theme <- theme(
    text         = element_text(size = 16, family = NATURE_FONT),
    axis.title   = element_text(size = 16, colour = "black"),
    axis.text    = element_text(size = 14, colour = "black"),
    legend.text  = element_text(size = 14),
    legend.title = element_blank(),
    plot.tag     = element_text(size = 20, face = "bold", family = NATURE_FONT),
    # Tight top/bottom margins so the panels fill the height; extra left margin so
    # the enlarged y-axis titles are not clipped.
    plot.margin       = margin(1, 4, 1, 6, unit = "mm"),
    legend.box.spacing = unit(1, "mm"),
    legend.margin      = margin(1, 6, 1, 6),
    # Larger key boxes so the enlarged legend glyphs have room.
    legend.key.height  = unit(1.4, "lines"),
    legend.key.width   = unit(2.6, "lines")
)

# Bigger legend glyphs for slides (re-declares the same guides as the base panels
# but with larger override.aes sizes; shapes/colours are preserved).
slide_guides <- guides(
    shape = "none",
    color = guide_legend(order = 1, title.position = "top", title.hjust = 0.5,
                         override.aes = list(linetype = 0, shape = c(17, 16, 15), size = 3.6)),
    fill  = guide_legend(order = 1, title.position = "top", title.hjust = 0.5,
                         override.aes = list(linetype = 0, shape = c(17, 16, 15), size = 3.6)),
    linetype = guide_legend(order = 2,
                            override.aes = list(shape = NA, color = "black", linewidth = 1.1))
)

# The panel-b y-axis title is long; at the larger slide font it would clip, so use
# a shorter wording for the slide figure only (the Nature figure keeps the full label).
fig2_slides <- (p_lonely +
                    labs(y = "Loneliness Score") + slide_guides) /
               (p_pref +
                    labs(y = "Preference for AI vs. Humans\n(0 = Humans, 100 = AI)") + slide_guides) +
    plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    (theme(legend.position = "top") + slide_theme)

ggsave("plots/figure2_loneliness_and_preference_slides.pdf", fig2_slides,
       width = 250, height = 230, units = "mm")
ggsave("plots/figure2_loneliness_and_preference_slides.png", fig2_slides,
       width = 250, height = 230, units = "mm", dpi = 300, bg = "white")

# Means summary
d_summary <- d_long %>%
    group_by(condition, timepoint, day) %>%
    summarise(mean_loneliness = mean(loneliness, na.rm = TRUE),
              n = n(), .groups = "drop")
print(d_summary, n = 300)

# With a t-test compare day 1 preference between before chatbot interaction vs. control
pref_control <- d_long$preference[d_long$condition == "control" & d_long$day == 1]
pref_text <- d_long$preference[d_long$condition == "both_text" & d_long$day == 1 & d_long$timepoint == "before"]
pref_voice <- d_long$preference[d_long$condition == "both_voice" & d_long$day == 1 & d_long$timepoint == "before"]
tt_text <- t.test(pref_control, pref_text, paired = FALSE)
cd_text <- cohen.d(pref_control, pref_text, paired = FALSE)
tt_voice <- t.test(pref_control, pref_voice, paired = FALSE)
cd_voice <- cohen.d(pref_control, pref_voice, paired = FALSE)
cat("\n--- Day 1 preference: control vs. before interaction ---\n")
print(paste0("Control: M = ", round(mean(pref_control), 2), " (", round(sd(pref_control), 2), ")"))
print(paste0("Text before: M = ", round(mean(pref_text), 2), " (", round(sd(pref_text), 2), ")"))
print(paste0("Voice before: M = ", round(mean(pref_voice), 2), " (", round(sd(pref_voice), 2), ")"))
print(tt_text); print(cd_text)
print(tt_voice); print(cd_voice)

###################### RESULTS ######################

###################### H1: Within-experience condition: timing * condition * day ######################

# Loneliness ~ Timing * Condition * Day + (1 | Participant ID); text vs voice only
d_long_exp <- d_long %>%
    filter(condition %in% c("both_text", "both_voice"))

d_long_exp$timing_f <- factor(d_long_exp$timepoint, levels = c("before", "after"))
d_long_exp$condition_f <- factor(d_long_exp$condition, levels = c("both_text", "both_voice"))
d_long_exp$day_n <- as.numeric(d_long_exp$day)

cat("\n========== H1: Loneliness ~ Timing * Condition * Day + (1|worker_id) ==========\n")
model_h1 <- lmer(loneliness ~ timing_f * condition_f * day_n + (1 | worker_id), data = d_long_exp)
print(summary(model_h1))
cat("\n--- H1 fixed-effect 95% CIs (Wald) ---\n"); print(round(confint(model_h1, method = "Wald"), 3))

# Cumulative paired t-tests: before vs. after, aggregated across all 10 days, per condition
cat("\n--- H1 cumulative paired t-tests (aggregated across all 10 days) ---\n")
for (cnd in c("both_text", "both_voice")) {
    cat("\n------- Condition:", cnd, "-------\n")
    sub <- d_long_exp %>% filter(condition == cnd)
    # Aggregate by worker_id within timepoint (mean across all 10 days)
    agg <- sub %>%
        group_by(worker_id, timepoint) %>%
        summarise(mean_lonely = mean(loneliness, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = timepoint, values_from = mean_lonely) %>%
        filter(!is.na(before) & !is.na(after))

    ttest <- t.test(agg$before, agg$after, paired = TRUE)
    cd <- cohen.d(agg$before, agg$after, paired = TRUE)
    print(ttest); print(cd)
    bf_report(agg$before, agg$after, paired = TRUE, label = cnd)
    print(paste0("Mean before: ", round(mean(agg$before), 2),
                 " | Mean after: ", round(mean(agg$after), 2)))
}

# Daily paired t-tests: before vs. after on each day, within each condition
cat("\n--- H1 daily paired t-tests ---\n")
print("day; condition; mbefore (sd) | mafter (sd); t(df); p; d")
for (cnd in c("both_text", "both_voice")) {
    for (i in DAYS) {
        sub <- d_long_exp %>%
            filter(condition == cnd & day == i) %>%
            dplyr::select(worker_id, timepoint, loneliness) %>%
            pivot_wider(names_from = timepoint, values_from = loneliness) %>%
            filter(!is.na(before) & !is.na(after))

        if (nrow(sub) < 3) { next }
        tt <- t.test(sub$before, sub$after, paired = TRUE)
        cd <- cohen.d(sub$before, sub$after, paired = TRUE)
        print(paste0(i, "; ", cnd, "; MBefore = ", round(mean(sub$before), 2),
                     " (", round(sd(sub$before), 2), ") vs. MAfter = ",
                     round(mean(sub$after), 2), " (", round(sd(sub$after), 2),
                     "), t(", round(tt$parameter, 2), ") = ", round(tt$statistic, 2),
                     ", p = ", round(tt$p.value, 4), ", d = ", round(cd$estimate, 2),
                     ", 95% CI [", round(cd$conf.int[1], 2), ", ", round(cd$conf.int[2], 2), "]"))
    }
}

# Robustness check (SI Section B, "Robustness of the within-session reduction excluding Day 1"):
# because the within-session reduction was largest on Day 1, re-estimate the H1 model on days
# 2-10 only, to verify the effect is not driven solely by the first session.
cat("\n--- H1 robustness: Timing * Condition * Day model excluding Day 1 (days 2-10) ---\n")
model_h1_noday1 <- lmer(loneliness ~ timing_f * condition_f * day_n + (1 | worker_id),
                        data = d_long_exp[d_long_exp$day_n >= 2, ])
print(summary(model_h1_noday1))
print(round(confint(model_h1_noday1, method = "Wald"), 2))

###################### H2: Post-interaction vs. control across days ######################

# Loneliness ~ Condition (control, text, voice) * Day + (1 | Participant ID)
# Use only post-interaction measurements (after) for text/voice, and the control measurement.
d_long_h2 <- d_long %>%
    filter(condition == "control" |
           (condition %in% c("both_text", "both_voice") & timepoint == "after"))

d_long_h2$condition_f <- factor(d_long_h2$condition, levels = c("control", "both_text", "both_voice"))
d_long_h2$day_n <- as.numeric(d_long_h2$day)

cat("\n========== H2: Loneliness ~ Condition * Day + (1|worker_id) (post-interaction vs control) ==========\n")
model_h2 <- lmer(loneliness ~ condition_f * day_n + (1 | worker_id), data = d_long_h2)
print(summary(model_h2))
cat("\n--- H2 fixed-effect 95% CIs (Wald) ---\n"); print(round(confint(model_h2, method = "Wald"), 3))

# Cumulative t-tests: aggregated loneliness over all 10 days, post-interaction vs control
cat("\n--- H2 cumulative t-tests (aggregated over all 10 days) ---\n")
agg_h2 <- d_long_h2 %>%
    group_by(worker_id, condition) %>%
    summarise(mean_lonely = mean(loneliness, na.rm = TRUE), .groups = "drop")

for (cnd in c("both_text", "both_voice")) {
    x <- agg_h2$mean_lonely[agg_h2$condition == "control"]
    y <- agg_h2$mean_lonely[agg_h2$condition == cnd]
    tt <- t.test(x, y, paired = FALSE)
    cd <- cohen.d(x, y, paired = FALSE)
    cat("\n------- Cumulative: Control vs. ", cnd, " -------\n")
    print(tt); print(cd)
    bf_report(x, y, paired = FALSE, label = paste0("control vs ", cnd))
    print(paste0("Mean control: ", round(mean(x), 2),
                 " | Mean ", cnd, ": ", round(mean(y), 2)))
}

# Text vs. voice on aggregated post-interaction loneliness LEVEL (a null comparison
# referenced in the manuscript). Reported with a Bayes factor to interpret the null.
cat("\n--- H2: post-interaction loneliness level, text vs. voice ---\n")
xt <- agg_h2$mean_lonely[agg_h2$condition == "both_text"]
yv <- agg_h2$mean_lonely[agg_h2$condition == "both_voice"]
tt <- t.test(xt, yv, paired = FALSE); cd <- cohen.d(xt, yv, paired = FALSE)
print(tt); print(cd)
bf_report(xt, yv, paired = FALSE, label = "text vs voice post-interaction level")
print(paste0("Mean text: ", round(mean(xt), 2), " | Mean voice: ", round(mean(yv), 2)))


# Two-sided t-tests comparing post-interaction loneliness vs. control per day, per condition
cat("\n--- H2 daily two-sided t-tests (post-interaction vs. control) ---\n")
for (cnd in c("both_text", "both_voice")) {
    cat("\n------- Condition:", cnd, "vs control -------\n")
    for (i in DAYS) {
        x <- d_long_h2$loneliness[d_long_h2$condition == "control" & d_long_h2$day == i]
        y <- d_long_h2$loneliness[d_long_h2$condition == cnd & d_long_h2$day == i]
        if (length(x) < 3 | length(y) < 3) { next }
        tt <- t.test(x, y, paired = FALSE)
        cd <- cohen.d(x, y, paired = FALSE)
        print(paste0(i, "; MControl = ", round(mean(x), 2), " (", round(sd(x), 2),
                     ") vs. MAfter[", cnd, "] = ", round(mean(y), 2), " (", round(sd(y), 2),
                     "), t(", round(tt$parameter, 2), ") = ", round(tt$statistic, 2),
                     ", p = ", round(tt$p.value, 4), ", d = ", round(cd$estimate, 2),
                     ", 95% CI [", round(cd$conf.int[1], 2), ", ", round(cd$conf.int[2], 2), "]"))
    }
}

###################### H3: Loneliness Difference ~ Condition (voice vs. text) * Day ######################

# Compute per-worker, per-day before-after difference (positive = reduction)
d_long_diff <- d_long %>%
    filter(condition %in% c("both_text", "both_voice")) %>%
    dplyr::select(worker_id, condition, day, timepoint, loneliness) %>%
    pivot_wider(names_from = timepoint, values_from = loneliness) %>%
    mutate(diff = before - after) %>%
    filter(!is.na(diff))

d_long_diff$condition_f <- factor(d_long_diff$condition, levels = c("both_text", "both_voice"))
d_long_diff$day_n <- as.numeric(d_long_diff$day)

cat("\n========== H3: Loneliness Difference ~ Condition (voice vs. text) * Day + (1|worker_id) ==========\n")
model_h3 <- lmer(diff ~ condition_f * day_n + (1 | worker_id), data = d_long_diff)
print(summary(model_h3))
cat("\n--- H3 fixed-effect 95% CIs (Wald) ---\n"); print(round(confint(model_h3, method = "Wald"), 3))

# Cumulative t-test on the aggregated reduction over all 10 days
cat("\n--- H3 cumulative t-test (aggregated reduction over all 10 days, voice vs. text) ---\n")
agg_h3 <- d_long_diff %>%
    group_by(worker_id, condition) %>%
    summarise(mean_diff = mean(diff, na.rm = TRUE), .groups = "drop")
x <- agg_h3$mean_diff[agg_h3$condition == "both_text"]
y <- agg_h3$mean_diff[agg_h3$condition == "both_voice"]
tt <- t.test(x, y, paired = FALSE); print(tt)
cd <- cohen.d(x, y, paired = FALSE); print(cd)
bf_report(x, y, paired = FALSE, label = "text vs voice reduction")
print(paste0("Mean diff text: ", round(mean(x), 2),
             " | Mean diff voice: ", round(mean(y), 2)))


# Daily two-sided t-tests on the difference, voice vs. text
cat("\n--- H3 daily two-sided t-tests (voice vs. text on loneliness reduction) ---\n")
for (i in DAYS) {
    x <- d_long_diff$diff[d_long_diff$condition == "both_text" & d_long_diff$day == i]
    y <- d_long_diff$diff[d_long_diff$condition == "both_voice" & d_long_diff$day == i]
    if (length(x) < 3 | length(y) < 3) { next }
    tt <- t.test(x, y, paired = FALSE)
    cd <- cohen.d(x, y, paired = FALSE)
    bf10 <- as.numeric(extractBF(ttestBF(x, y, paired = FALSE))$bf)
    print(paste0(i, "; MText = ", round(mean(x), 2), " (", round(sd(x), 2),
                 ") vs. MVoice = ", round(mean(y), 2), " (", round(sd(y), 2),
                 "), t(", round(tt$parameter, 2), ") = ", round(tt$statistic, 2),
                 ", p = ", round(tt$p.value, 4), ", d = ", round(cd$estimate, 2),
                 ", 95% CI [", round(cd$conf.int[1], 2), ", ", round(cd$conf.int[2], 2), "]",
                 ", BF10 = ", signif(bf10, 3)))
}

###################### H4: Preference for AI vs. humans, interaction vs. control ######################

cat("\n========== H4: Preference for AI vs humans, interaction vs control ==========\n")

# Cumulative preference tests (Table 3g in the main text): per-participant mean preference
# across the ten days (post-interaction rating for text/voice; the single daily rating for
# control), compared between conditions with Welch t-tests, Cohen's d, and BF10.
cat("\n--- H4 cumulative preference t-tests (aggregated across all 10 days; Table 3g) ---\n")
d_pref_cum <- d
d_pref_cum$pref_cum <- ifelse(d_pref_cum$condition == "control",
                              d_pref_cum$preference_before_1, d_pref_cum$preference_after_1)
agg_pref_cum <- aggregate(pref_cum ~ worker_id + condition, data = d_pref_cum, FUN = mean, na.rm = TRUE)
for (pair in list(c("control", "both_text"), c("control", "both_voice"), c("both_text", "both_voice"))) {
    x <- agg_pref_cum$pref_cum[agg_pref_cum$condition == pair[1]]
    y <- agg_pref_cum$pref_cum[agg_pref_cum$condition == pair[2]]
    tt <- t.test(x, y); cd <- cohen.d(x, y); bf <- extractBF(ttestBF(x, y))$bf
    cat(sprintf("%s (M=%.2f) vs %s (M=%.2f): t(%.1f) = %.2f, p = %.4g, d = %.2f [%.2f, %.2f], BF10 = %.3g\n",
                pair[1], mean(x), pair[2], mean(y), tt$parameter, tt$statistic, tt$p.value,
                cd$estimate, cd$conf.int[1], cd$conf.int[2], bf))
}


# For control: preference is the single measurement (preference_before_1).
# For text/voice: use post-interaction preference (after).
d_pref <- d_long %>%
    filter(condition == "control" |
           (condition %in% c("both_text", "both_voice") & timepoint == "after")) %>%
    filter(!is.na(preference))

cat("\n--- H4 daily two-sided t-tests (preference: interaction conditions vs. control) ---\n")
for (cnd in c("both_text", "both_voice")) {
    cat("\n------- Condition:", cnd, "vs control -------\n")
    for (i in DAYS) {
        x <- d_pref$preference[d_pref$condition == "control" & d_pref$day == i]
        y <- d_pref$preference[d_pref$condition == cnd & d_pref$day == i]
        if (length(x) < 3 | length(y) < 3) { next }
        tt <- t.test(x, y, paired = FALSE)
        cd <- cohen.d(x, y, paired = FALSE)
        print(paste0(i, "; MControl = ", round(mean(x), 2), " (", round(sd(x), 2),
                     ") vs. M", cnd, " = ", round(mean(y), 2), " (", round(sd(y), 2),
                     "), t(", round(tt$parameter, 2), ") = ", round(tt$statistic, 2),
                     ", p = ", round(tt$p.value, 4), ", d = ", round(cd$estimate, 2),
                     ", 95% CI [", round(cd$conf.int[1], 2), ", ", round(cd$conf.int[2], 2), "]"))
    }
}

# Exploratory (not pre-registered): does the preference gap between the interaction
# conditions and control GROW over the ten days? Test the Condition x Day interaction.
cat("\n--- H4 exploratory: preference ~ Condition * Day (does the gap grow?) ---\n")
d_pref_trend <- d_pref
d_pref_trend$condition_f <- factor(d_pref_trend$condition, levels = c("control", "both_text", "both_voice"))
d_pref_trend$day_n <- as.numeric(d_pref_trend$day)
model_pref <- lmer(preference ~ condition_f * day_n + (1 | worker_id), data = d_pref_trend)
print(summary(model_pref))
cat("\n--- preference trend 95% CIs (Wald) ---\n")
print(round(confint(model_pref, method = "Wald")[c("condition_fboth_text:day_n",
                                                   "condition_fboth_voice:day_n", "day_n"), ], 3))


# Exploratory (not pre-registered): the three mixed-effects models reported in Table B6
# (SI Section B), decomposing the preference trends above:
# (i)  the within-session change in preference (after - before) in the text and voice
#      conditions, as a function of condition, day, and their interaction — the transient
#      daily ("substitution") shift toward AI;
# (ii) control participants' daily preference as a function of day — the control decline; and
# (iii) the standing (before-session) preference in the text and voice conditions — whether
#      the daily shift accumulates across days (it does not).
d_b6 <- d
d_b6$day_n <- as.numeric(d_b6$day)
d_b6$condition_f <- factor(d_b6$condition, levels = c("both_text", "both_voice"))
d_b6_ia <- d_b6[d_b6$condition %in% c("both_text", "both_voice"), ]
d_b6_ia$pref_delta <- d_b6_ia$preference_after_1 - d_b6_ia$preference_before_1

cat("\n--- Mean within-session preference change (after - before), by condition ---\n")
print(round(tapply(d_b6_ia$pref_delta, d_b6_ia$condition, mean, na.rm = TRUE), 2))

cat("\n--- Table B6 (i): within-session change ~ Condition * Day + (1 | participant) ---\n")
m_b6_within <- lmer(pref_delta ~ condition_f * day_n + (1 | worker_id), data = d_b6_ia)
print(summary(m_b6_within))
print(round(confint(m_b6_within, method = "Wald"), 2))

cat("\n--- Table B6 (ii): control daily preference ~ Day + (1 | participant) ---\n")
m_b6_control <- lmer(preference_before_1 ~ day_n + (1 | worker_id),
                     data = d_b6[d_b6$condition == "control", ])
print(summary(m_b6_control))
print(round(confint(m_b6_control, method = "Wald"), 2))

cat("\n--- Table B6 (iii): standing (before-session) preference ~ Condition * Day + (1 | participant) ---\n")
m_b6_standing <- lmer(preference_before_1 ~ condition_f * day_n + (1 | worker_id), data = d_b6_ia)
print(summary(m_b6_standing))
print(round(confint(m_b6_standing, method = "Wald"), 2))

cat("\n--- Table B6 (iv): post-interaction preference ~ Condition * Day + (1 | participant) ---\n")
m_b6_post <- lmer(preference_after_1 ~ condition_f * day_n + (1 | worker_id), data = d_b6_ia)
print(summary(m_b6_post))
print(round(confint(m_b6_post, method = "Wald"), 2))
# Simple day slope for the voice condition (releveled), quoted alongside the text
# slope in the manuscript and SI Section B.
d_b6_ia$condition_rev <- relevel(d_b6_ia$condition_f, ref = "both_voice")
m_b6_post_v <- lmer(preference_after_1 ~ condition_rev * day_n + (1 | worker_id), data = d_b6_ia)
cat("\n--- (iv) voice-condition simple day slope ---\n")
print(round(summary(m_b6_post_v)$coefficients["day_n", , drop = FALSE], 4))
print(round(confint(m_b6_post_v, method = "Wald")["day_n", ], 3))

###################### CHATBOT MESSAGE METRICS (words/message & messages sent) ######################

# Merge the chatbot conversation logs (voice_longitudinal DB export) with the
# analysis sample to get, per condition: average words per participant-sent
# message and number of messages sent per participant.
cat("\n========== Chatbot message metrics (words per message & messages sent) ==========\n")

# NOTE: reads the raw chatbot conversation logs (voice_longitudinal/ DB export), which
# are not committed to the public repo (participant identifiers / large files). Guarded
# so the script runs from committed data; the message metrics only build when present.
# The downstream LGM, Bayes-factor tables, and PSM (below) run from the committed
# d1-d10 data and are unaffected by this guard.
if (file.exists("voice_longitudinal/mturk_app_message.csv") &&
    file.exists("voice_longitudinal/mturk_app_participant.csv")) {

msg <- read.csv("voice_longitudinal/mturk_app_message.csv", stringsAsFactors = FALSE)
ppt <- read.csv("voice_longitudinal/mturk_app_participant.csv", stringsAsFactors = FALSE)

# Attach worker_id to every message via room_id. Chatbot turns have an empty
# user_id, so room_id is the common key; rooms are unique per participant and the
# room_id join reproduces the user_id join exactly for human messages.
room_lookup <- ppt[ppt$room_id != "" & !is.na(ppt$room_id), c("room_id", "worker_id")]
room_lookup <- room_lookup[!duplicated(room_lookup$room_id), ]
msg_m <- merge(msg, room_lookup, by = "room_id", all.x = FALSE)

# Word count per message (whitespace-delimited; empty content counts as 0).
msg_m$content <- trimws(msg_m$content)
msg_m$n_words <- lengths(strsplit(msg_m$content, "\\s+"))
msg_m$n_words[msg_m$content == ""] <- 0

# Restrict to the analysis sample and take condition from the Qualtrics data (d).
# Use only non-empty condition labels and keep one row per worker (condition is
# fixed at assignment; blank entries on some day-rows would otherwise duplicate workers).
cond_lookup <- unique(d[d$condition %in% c("both_text", "both_voice", "control"),
                        c("worker_id", "condition")])
cond_lookup <- cond_lookup[!duplicated(cond_lookup$worker_id), ]
msg_m <- merge(msg_m, cond_lookup, by = "worker_id", all.x = FALSE)
print(paste0("Participants with message data in the analysis sample: ",
             length(unique(msg_m$worker_id)), " (of ", nrow(cond_lookup), " in d)"))

# Summarize per condition x sender_type (human = participant-sent, chatbot = AI).
# Words per message is pooled (total words / total messages); messages-per-participant
# is the per-worker message count averaged across workers in the condition.
per_ppt_msg <- msg_m %>%
    group_by(worker_id, condition, sender_type) %>%
    summarise(n_messages = n(), .groups = "drop") %>%
    group_by(condition, sender_type) %>%
    summarise(n_participants        = n(),
              messages_per_ppt_mean = mean(n_messages),
              messages_per_ppt_sd   = sd(n_messages), .groups = "drop")

words_msg <- msg_m %>%
    group_by(condition, sender_type) %>%
    summarise(total_messages     = n(),
              words_per_msg_mean  = mean(n_words),
              words_per_msg_sd    = sd(n_words), .groups = "drop")

msg_metrics <- merge(per_ppt_msg, words_msg, by = c("condition", "sender_type"))
print(as.data.frame(msg_metrics))

} else {
  message("[reproducibility] Skipping chatbot message metrics: voice_longitudinal/ exports not present (excluded from the public repo).")
}

######################## EXPLORATORY: LATENT GROWTH CURVE MODELING ########################

cat("\n========== Exploratory: Latent Growth Curve Modeling ==========\n")

# Control provides a single daily rating WITHOUT an interaction, which is conceptually a
# PRE-interaction measure. We therefore (a) compare the three groups on pre-interaction
# (before-session) loneliness -- a like-for-like test of CARRYOVER to the next day's
# baseline -- and (b) compare the two interaction conditions on their POST-interaction
# loneliness (MODALITY). The pre-interaction trajectory is non-linear (sharp Day 1 -> Day 2
# drop), so the carryover models use a latent-basis specification (freely estimated shape).

mk_wide <- function(dat) {
    dat %>% dplyr::select(worker_id, condition, day, loneliness) %>%
        pivot_wider(names_from = day, values_from = loneliness, names_prefix = "L") %>%
        filter(if_all(starts_with("L"), ~ !is.na(.)))
}
lin_i <- "i =~ 1*L1 + 1*L2 + 1*L3 + 1*L4 + 1*L5 + 1*L6 + 1*L7 + 1*L8 + 1*L9 + 1*L10"
lin_s <- "s =~ 0*L1 + 1*L2 + 2*L3 + 3*L4 + 4*L5 + 5*L6 + 6*L7 + 7*L8 + 8*L9 + 9*L10"
lb_s  <- "s =~ 0*L1 + L2 + L3 + L4 + L5 + L6 + L7 + L8 + L9 + 1*L10"  # latent-basis (free shape)

# (a) CARRYOVER: pre-interaction loneliness, all three conditions (control's daily rating is a 'before')
pre_wide <- mk_wide(d_long %>% filter(timepoint == "before"))
pre_wide$cond_text  <- as.numeric(pre_wide$condition == "both_text")
pre_wide$cond_voice <- as.numeric(pre_wide$condition == "both_voice")
print(paste0("Carryover LGM N: ", nrow(pre_wide))); print(table(pre_wide$condition))

cat("\n--- Carryover: linear fit (expected poor, motivates latent basis) ---\n")
print(fitMeasures(growth(paste(lin_i, lin_s, sep = "\n"), data = pre_wide), c("cfi", "rmsea", "srmr")))

cat("\n--- Carryover: unconditional latent-basis (individual differences / variance components) ---\n")
print(summary(growth(paste(lin_i, lb_s, sep = "\n"), data = pre_wide),
              fit.measures = TRUE, standardized = TRUE, ci = TRUE))

cat("\n--- Carryover: conditional latent-basis (control reference) ---\n")
print(summary(growth(paste(lin_i, lb_s, "i ~ cond_text + cond_voice", "s ~ cond_text + cond_voice", sep = "\n"),
                     data = pre_wide), fit.measures = TRUE, ci = TRUE))

cat("\n--- Carryover: direct voice vs. text contrast (latent-basis, text reference) ---\n")
pre_iv <- pre_wide %>% filter(condition %in% c("both_text", "both_voice"))
pre_iv$voice <- as.numeric(pre_iv$condition == "both_voice")
print(summary(growth(paste(lin_i, lb_s, "i ~ voice", "s ~ voice", sep = "\n"), data = pre_iv),
              fit.measures = TRUE, ci = TRUE))

# (b) MODALITY: post-interaction loneliness, voice vs. text (linear; post trajectory is approximately linear)
post_wide <- mk_wide(d_long %>% filter(timepoint == "after", condition %in% c("both_text", "both_voice")))
post_wide$voice <- as.numeric(post_wide$condition == "both_voice")
print(paste0("Modality LGM N: ", nrow(post_wide))); print(table(post_wide$condition))
cat("\n--- Modality: post-interaction voice vs. text (linear) ---\n")
print(summary(growth(paste(lin_i, lin_s, "i ~ voice", "s ~ voice", sep = "\n"), data = post_wide),
              fit.measures = TRUE, ci = TRUE))

# Pre-interaction (carryover) t-tests corroborating the LGM (Table B7): control's single
# daily rating serves as its pre-interaction value.
cat("\n--- Carryover t-tests on pre-interaction loneliness (Table B7) ---\n")
pre_long <- d_long %>% filter(timepoint == "before")
carry_tt <- function(dat, lab) {
    g <- function(cn) dat$loneliness[dat$condition == cn]
    for (pair in list(c("both_text", "control"), c("both_voice", "control"), c("both_voice", "both_text"))) {
        x <- g(pair[1]); y <- g(pair[2]); tt <- t.test(x, y); cd <- cohen.d(x, y)
        bf <- BayesFactor::extractBF(BayesFactor::ttestBF(x, y))$bf
        cat(sprintf("%-9s %-12s (M=%.1f) vs %-10s (M=%.1f): t(%.0f)=%.2f, p=%.4g, d=%.2f, BF10=%.3g\n",
                    lab, pair[1], mean(x, na.rm = TRUE), pair[2], mean(y, na.rm = TRUE),
                    tt$parameter, tt$statistic, tt$p.value, cd$estimate, bf))
    }
}

carry_tt(pre_long %>% filter(day == 1), "Day 1")
agg_pre <- pre_long %>% filter(day >= 2) %>% group_by(worker_id, condition) %>%
    summarise(loneliness = mean(loneliness, na.rm = TRUE), .groups = "drop")
carry_tt(agg_pre, "Days 2-10")
carry_tt(pre_long %>% filter(day == 10), "Day 10")

# Day-by-day pre-interaction loneliness vs. control, for text and voice (Table B8)
cat("\n--- Day-by-day pre-interaction loneliness vs. control (Table B8) ---\n")
for (cn in c("both_text", "both_voice")) {
    for (dy in DAYS) {
        ctl <- pre_long$loneliness[pre_long$condition == "control" & pre_long$day == dy]
        pr  <- pre_long$loneliness[pre_long$condition == cn      & pre_long$day == dy]
        tt <- t.test(ctl, pr); cd <- cohen.d(ctl, pr)
        cat(sprintf("Day %2d %-10s: control %.2f (%.2f) vs pre %.2f (%.2f), t(%.1f)=%.2f, p=%.4g, d=%.2f [%.2f, %.2f]\n",
                    dy, cn, mean(ctl), sd(ctl), mean(pr), sd(pr),
                    tt$parameter, tt$statistic, tt$p.value, cd$estimate, cd$conf.int[1], cd$conf.int[2]))
    }
}

# Daily Bayes factors for the daily tables (B1 within-condition, B2 post-vs-control,
# B4 preference-vs-control, B8 pre-vs-control; paired BF for B1, independent otherwise) and
# the Day-1 baseline-preference comparison (B5).
cat("\n--- Daily Bayes factors (Tables B1, B2, B4, B5, B8) ---\n")
d$lb <- rowMeans(sapply(d[, c("lonely_before_1", "lonely_before_2", "lonely_before_3")], as.numeric), na.rm = TRUE)
d$la <- rowMeans(sapply(d[, c("lonely_after_1", "lonely_after_2", "lonely_after_3")], as.numeric), na.rm = TRUE)
d$pb <- as.numeric(d$preference_before_1); d$pa <- as.numeric(d$preference_after_1)
bfi <- function(x, y) BayesFactor::extractBF(BayesFactor::ttestBF(x, y))$bf
bfp <- function(x, y) BayesFactor::extractBF(BayesFactor::ttestBF(x, y, paired = TRUE))$bf
for (cn in c("both_text", "both_voice")) {
    for (dy in DAYS) {
        sub <- d[d$condition == cn & d$day == dy, ]
        ctl_lb <- d$lb[d$condition == "control" & d$day == dy]; ctl_lb <- ctl_lb[!is.na(ctl_lb)]
        ctl_pb <- d$pb[d$condition == "control" & d$day == dy]; ctl_pb <- ctl_pb[!is.na(ctl_pb)]
        mb <- !is.na(sub$lb) & !is.na(sub$la)
        # %.5g (not %.3g): enough significant digits to read the 2-dp Bayes factors
        # reported in Tables B2/B4/B8, while keeping very large BFs in compact
        # scientific notation.
        cat(sprintf("Day %2d %-10s | B1 within=%.5g | B2 post-vs-ctrl=%.5g | B4 pref-vs-ctrl=%.5g | B8 pre-vs-ctrl=%.5g\n",
            dy, cn, bfp(sub$lb[mb], sub$la[mb]),
            bfi(ctl_lb, sub$la[!is.na(sub$la)]), bfi(ctl_pb, sub$pa[!is.na(sub$pa)]),
            bfi(ctl_lb, sub$lb[!is.na(sub$lb)])))
    }
}
for (cn in c("both_text", "both_voice")) {
    x <- d$pb[d$condition == "control" & d$day == 1]; y <- d$pb[d$condition == cn & d$day == 1]
    cat(sprintf("B5 Day-1 baseline preference %-10s BF10=%.3g\n", cn, bfi(x[!is.na(x)], y[!is.na(y)])))
}

###################### APPENDIX: PROPENSITY SCORE MATCHING (SELECTION ROBUSTNESS) ######################

# Robustness check for the selection concern (chatbot conditions show higher baseline
# AI preference than control). We match each interaction condition to control on day-1
# baseline covariates (age, gender, AI experience, baseline loneliness, baseline AI
# preference) via 1:1 nearest-neighbor propensity-score matching, then re-run the key
# control comparisons (H2 post-interaction loneliness, H4 AI-vs-humans preference) on the
# matched samples. H1 (within-condition before/after) and H3 (text vs. voice) do not
# involve the control group and are not subject to this chatbot-vs-control selection.
cat("\n========== Appendix: Propensity Score Matching (selection robustness) ==========\n")

base <- d[d$day == 1, c("worker_id", "condition", "age", "gender", "edu", "ai_companion_exp",
                        "lonely_before", "preference_before_1")]
base$age         <- as.numeric(base$age)
base$female      <- as.numeric(base$gender == "2")
base$education   <- as.numeric(base$edu)
base$ai_exp      <- as.numeric(base$ai_companion_exp == 1)
base$pref_before <- as.numeric(base$preference_before_1)
base <- base[!is.na(base$age) & !is.na(base$female) & !is.na(base$education) & !is.na(base$ai_exp) &
             !is.na(base$lonely_before) & !is.na(base$pref_before), ]

# Per-worker outcomes (aggregated across the 10 days)
agg_lon <- d_long_h2 %>% group_by(worker_id, condition) %>%
    summarise(lonely = mean(loneliness, na.rm = TRUE), .groups = "drop")
agg_prf <- d_pref %>% group_by(worker_id, condition) %>%
    summarise(pref = mean(preference, na.rm = TRUE), .groups = "drop")

run_psm <- function(cond) {
    sub <- base[base$condition %in% c("control", cond), ]
    sub$treat <- as.numeric(sub$condition == cond)
    # link = "linear.logit" puts the propensity score on the logit scale, so the
    # caliper is 0.2 SD of the propensity-score LOGIT — the width recommended by
    # Austin (2011). Without it, distance = "glm" leaves the score on the probability
    # scale and caliper = 0.2 would be 0.2 SD of the probability, a different width.
    mm <- MatchIt::matchit(treat ~ age + female + education + ai_exp + lonely_before + pref_before,
                           data = sub, method = "nearest", distance = "glm",
                           link = "linear.logit", ratio = 1, caliper = 0.2)
    cat("\n----- ", cond, " vs control -----\n", sep = "")
    # binary = "std" so binary covariates (female, ai_exp) are also reported as standardized
    # mean differences — these are the SMD values in Table B9 of the SI.
    print(cobalt::bal.tab(mm, un = TRUE, stats = "mean.diffs", binary = "std")$Balance[, c("Diff.Un", "Diff.Adj")])
    ids <- MatchIt::match.data(mm)$worker_id

    x <- agg_lon$lonely[agg_lon$condition == "control" & agg_lon$worker_id %in% ids]
    y <- agg_lon$lonely[agg_lon$condition == cond      & agg_lon$worker_id %in% ids]
    tt <- t.test(x, y); cd <- cohen.d(x, y)
    cat(sprintf("H2 loneliness (matched, %d pairs): Control M=%.2f vs %s M=%.2f, t(%.1f)=%.2f, p=%.4g, d=%.2f, 95%% CI [%.2f, %.2f]\n",
                length(ids) / 2, mean(x), cond, mean(y), tt$parameter, tt$statistic, tt$p.value, cd$estimate, cd$conf.int[1], cd$conf.int[2]))

    px <- agg_prf$pref[agg_prf$condition == "control" & agg_prf$worker_id %in% ids]
    py <- agg_prf$pref[agg_prf$condition == cond      & agg_prf$worker_id %in% ids]
    tt2 <- t.test(px, py); cd2 <- cohen.d(px, py)
    cat(sprintf("H4 preference (matched, %d pairs): Control M=%.2f vs %s M=%.2f, t(%.1f)=%.2f, p=%.4g, d=%.2f, 95%% CI [%.2f, %.2f]\n",
                length(ids) / 2, mean(px), cond, mean(py), tt2$parameter, tt2$statistic, tt2$p.value, cd2$estimate, cd2$conf.int[1], cd2$conf.int[2]))
}

run_psm("both_text")
run_psm("both_voice")

# Robustness of the matched-sample conclusions to the matching specification
# (SI Section B). The primary analysis uses a 0.2-SD-of-logit caliper (Austin 2011);
# here we re-estimate H2 (post-interaction loneliness) and H4 (AI preference) under
# alternative specifications: a tighter 0.1 caliper, an added 0.25-SD caliper on the
# key confounder (baseline AI preference), and Mahalanobis matching on that confounder
# within the caliper. The point is not to pick a spec but to show the effect does not
# hinge on one: every comparison stays significant across all of them.
cat("\n--- PSM robustness: H2/H4 across matching specifications ---\n")
psm_spec <- function(cond, ...) {
    s <- base[base$condition %in% c("control", cond), ]; s$treat <- as.numeric(s$condition == cond)
    m <- MatchIt::matchit(treat ~ age + female + education + ai_exp + lonely_before + pref_before,
                          data = s, method = "nearest", distance = "glm", link = "linear.logit",
                          ratio = 1, ...)
    ids <- MatchIt::match.data(m)$worker_id
    lx <- agg_lon$lonely[agg_lon$condition == "control" & agg_lon$worker_id %in% ids]
    ly <- agg_lon$lonely[agg_lon$condition == cond      & agg_lon$worker_id %in% ids]
    px <- agg_prf$pref[agg_prf$condition == "control" & agg_prf$worker_id %in% ids]
    py <- agg_prf$pref[agg_prf$condition == cond      & agg_prf$worker_id %in% ids]
    c(pairs = length(ids) / 2,
      H2_p = t.test(lx, ly)$p.value, H2_d = cohen.d(lx, ly)$estimate,
      H4_p = t.test(px, py)$p.value, H4_d = cohen.d(px, py)$estimate)
}
specs <- list("0.2 caliper (primary)"        = list(caliper = 0.2),
              "0.1 caliper"                   = list(caliper = 0.1),
              "0.2 + 0.25-SD caliper on pref" = list(caliper = c(0.2, pref_before = 0.25)),
              "Mahalanobis on pref in caliper"= list(caliper = 0.2, mahvars = ~ pref_before))
for (cond in c("both_text", "both_voice")) {
    for (nm in names(specs)) {
        r <- do.call(psm_spec, c(list(cond), specs[[nm]]))
        cat(sprintf("%-11s | %-30s pairs=%3.0f | H2 d=%.2f p=%.3g | H4 d=%.2f p=%.3g\n",
                    cond, nm, r["pairs"], abs(r["H2_d"]), r["H2_p"], abs(r["H4_d"]), r["H4_p"]))
    }
}

# Matched-sample replication (SI Section B, end of the PSM section): re-estimate the
# H1-H4 analyses on the matched samples, mirroring the "Full-sample replication"
# structure. H1/H3 involve only the interaction conditions and are not subject to
# chatbot-vs-control selection, but are re-estimated on the matched interaction
# participants for completeness.
cat("\n--- Matched-sample replication (H1-H4) ---\n")
bfi_r <- function(x, y) BayesFactor::extractBF(BayesFactor::ttestBF(x, y))$bf
bfp_r <- function(x, y) BayesFactor::extractBF(BayesFactor::ttestBF(x, y, paired = TRUE))$bf
match_ids <- function(cond) {
    s <- base[base$condition %in% c("control", cond), ]; s$treat <- as.numeric(s$condition == cond)
    m <- MatchIt::matchit(treat ~ age + female + education + ai_exp + lonely_before + pref_before,
                          data = s, method = "nearest", distance = "glm", link = "linear.logit",
                          ratio = 1, caliper = 0.2)
    MatchIt::match.data(m)$worker_id
}
ids_t <- match_ids("both_text"); ids_v <- match_ids("both_voice")

for (cond in c("both_text", "both_voice")) {
    ids <- if (cond == "both_text") ids_t else ids_v
    cat("\n=====", cond, "vs matched control (pairs =", length(ids) / 2, ") =====\n")
    ## H2 mixed model on the matched sample
    sub <- d_long_h2[d_long_h2$worker_id %in% ids & d_long_h2$condition %in% c("control", cond), ]
    sub$condition_f <- factor(sub$condition, levels = c("control", cond))
    m2 <- lmer(loneliness ~ condition_f * day_n + (1 | worker_id), data = sub)
    co <- summary(m2)$coefficients; ci <- confint(m2, method = "Wald")
    rn <- paste0("condition_f", cond)
    cat(sprintf("H2 model: condition b = %.2f [%.2f, %.2f], t = %.2f, p = %.3g | cond x day b = %.2f, p = %.3g\n",
                co[rn, "Estimate"], ci[rn, 1], ci[rn, 2], co[rn, "t value"], co[rn, "Pr(>|t|)"],
                co[paste0(rn, ":day_n"), "Estimate"], co[paste0(rn, ":day_n"), "Pr(>|t|)"]))
    ## H2 cumulative + daily
    x <- agg_lon$lonely[agg_lon$condition == "control" & agg_lon$worker_id %in% ids]
    y <- agg_lon$lonely[agg_lon$condition == cond & agg_lon$worker_id %in% ids]
    tt <- t.test(x, y); cd <- cohen.d(x, y)
    cat(sprintf("H2 cumulative: %.2f vs %.2f, t(%.1f)=%.2f, p=%.3g, d=%.2f [%.2f,%.2f], BF10=%.3g\n",
                mean(y), mean(x), tt$parameter, tt$statistic, tt$p.value, cd$estimate,
                cd$conf.int[1], cd$conf.int[2], bfi_r(x, y)))
    ps <- ds <- c()
    for (dy in DAYS) {
        xd <- d_long_h2$loneliness[d_long_h2$worker_id %in% ids & d_long_h2$condition == "control" & d_long_h2$day == dy]
        yd <- d_long_h2$loneliness[d_long_h2$worker_id %in% ids & d_long_h2$condition == cond & d_long_h2$day == dy]
        tt <- t.test(xd, yd); ps <- c(ps, tt$p.value); ds <- c(ds, cohen.d(xd, yd)$estimate)
    }
    cat(sprintf("H2 daily: d = %.2f-%.2f, max p = %.4g (%d/10 days p<.05)\n", min(ds), max(ds), max(ps), sum(ps < .05)))
    ## H4 cumulative + daily
    px <- agg_prf$pref[agg_prf$condition == "control" & agg_prf$worker_id %in% ids]
    py <- agg_prf$pref[agg_prf$condition == cond & agg_prf$worker_id %in% ids]
    tt <- t.test(px, py); cd <- cohen.d(px, py)
    cat(sprintf("H4 cumulative: %.2f vs %.2f, t(%.1f)=%.2f, p=%.3g, |d|=%.2f [%.2f,%.2f], BF10=%.3g\n",
                mean(py), mean(px), tt$parameter, tt$statistic, tt$p.value, abs(cd$estimate),
                abs(cd$conf.int[2]), abs(cd$conf.int[1]), bfi_r(px, py)))
    ps <- ds <- c()
    for (dy in DAYS) {
        xd <- d_pref$preference[d_pref$worker_id %in% ids & d_pref$condition == "control" & d_pref$day == dy]
        yd <- d_pref$preference[d_pref$worker_id %in% ids & d_pref$condition == cond & d_pref$day == dy]
        tt <- t.test(xd, yd); ps <- c(ps, tt$p.value); ds <- c(ds, abs(cohen.d(xd, yd)$estimate))
        cat(sprintf("  H4 day %2d: p=%.4f |d|=%.2f\n", dy, tt$p.value, abs(cohen.d(xd, yd)$estimate)))
    }
    cat(sprintf("H4 daily: |d| = %.2f-%.2f, max p = %.4g (%d/10 days p<.05)\n", min(ds), max(ds), max(ps), sum(ps < .05)))
    ## H1 cumulative within-session reduction, matched interaction participants
    tids <- ids[ids %in% base$worker_id[base$condition == cond]]
    agg <- d %>% filter(worker_id %in% tids) %>% group_by(worker_id) %>%
        summarise(lb = mean(lonely_before, na.rm = TRUE), la = mean(lonely_after, na.rm = TRUE))
    tt <- t.test(agg$lb, agg$la, paired = TRUE); cd <- cohen.d(agg$lb, agg$la, paired = TRUE)
    cat(sprintf("H1 cumulative (n=%d): %.2f -> %.2f, t(%d)=%.2f, p=%.3g, d=%.2f [%.2f,%.2f], BF10=%.3g\n",
                nrow(agg), mean(agg$lb), mean(agg$la), tt$parameter, tt$statistic, tt$p.value,
                cd$estimate, cd$conf.int[1], cd$conf.int[2], bfp_r(agg$lb, agg$la)))
}

## H3: voice vs text among matched interaction participants
tid <- ids_t[ids_t %in% base$worker_id[base$condition == "both_text"]]
vid <- ids_v[ids_v %in% base$worker_id[base$condition == "both_voice"]]
aggd <- d %>% filter(worker_id %in% c(tid, vid)) %>% group_by(worker_id, condition) %>%
    summarise(red = mean(lonely_before - lonely_after, na.rm = TRUE), .groups = "drop")
x <- aggd$red[aggd$condition == "both_text"]; y <- aggd$red[aggd$condition == "both_voice"]
tt <- t.test(x, y); cd <- cohen.d(x, y)
cat(sprintf("\nH3 matched (text n=%d, voice n=%d): %.2f vs %.2f, t(%.1f)=%.2f, p=%.3g, d=%.2f [%.2f,%.2f], BF10=%.3g\n",
            length(x), length(y), mean(x), mean(y), tt$parameter, tt$statistic, tt$p.value,
            cd$estimate, cd$conf.int[1], cd$conf.int[2], bfi_r(x, y)))


############################################################################
######### Comparison with Lasgaard et al. (2026) meta-analysis ############
############################################################################
# Manuscript (General Discussion): Study 2's contribution to the reported
# range "the immediate reductions in our studies (g = 0.38-0.59)".
#
# Same metric as Lasgaard et al. use for RCTs (Hedges' g on end-of-intervention
# scores, intervention vs. control; negative = reduction). The control
# condition rates loneliness once per day (no interaction), so the Day-10
# between-group contrast is each condition's post-interaction rating vs. the
# control's same-day rating -- the same contrast as H2. For completeness we
# also report the stricter durability contrast (Day-10 PRE-session vs.
# control), which is why the manuscript claims only *immediate* reductions.

hedges_j <- function(df) 1 - 3 / (4 * df - 1)
hedges_g_between <- function(x1, x0) {
    n1 <- length(x1); n0 <- length(x0); df <- n1 + n0 - 2
    sp <- sqrt(((n1 - 1) * var(x1) + (n0 - 1) * var(x0)) / df)
    g  <- (mean(x1) - mean(x0)) / sp * hedges_j(df)
    se <- sqrt((n1 + n0) / (n1 * n0) + g^2 / (2 * df))
    c(g = g, lo = g - qnorm(.975) * se, hi = g + qnorm(.975) * se)
}

x_ctl10 <- d_long$loneliness[d_long$condition == "control" & d_long$day == 10 &
                             d_long$timepoint == "before"]
for (cn in c("both_text", "both_voice")) {
    g_post <- hedges_g_between(
        d_long$loneliness[d_long$condition == cn & d_long$day == 10 &
                          d_long$timepoint == "after"], x_ctl10)
    g_pre <- hedges_g_between(
        d_long$loneliness[d_long$condition == cn & d_long$day == 10 &
                          d_long$timepoint == "before"], x_ctl10)
    cat(sprintf("%s Day-10 vs control: post-session g = %.2f [%.2f, %.2f]; pre-session g = %.2f [%.2f, %.2f]\n",
                cn, g_post[1], g_post[2], g_post[3], g_pre[1], g_pre[2], g_pre[3]))
}

lasgaard_s1_file <- "../lit/Supplemental Materials Table S1 AMP-2022-0748.xlsx"
if (file.exists(lasgaard_s1_file) && requireNamespace("readxl", quietly = TRUE)) {
    ts1 <- as.data.frame(readxl::read_excel(lasgaard_s1_file, sheet = 1, skip = 1))
    ts1 <- ts1[!is.na(ts1[["Study design"]]), ]
    es <- as.numeric(sub("^\\s*(-?[0-9.]+)\\s*\\[.*", "\\1",
                         gsub(",", ".", ts1[["Short-term effect [95% CI]"]])))
    rct_es <- es[ts1[["Study design"]] == "RCT" & !is.na(es)]
    for (cn in c("both_text", "both_voice")) {
        g <- hedges_g_between(
            d_long$loneliness[d_long$condition == cn & d_long$day == 10 &
                              d_long$timepoint == "after"], x_ctl10)[1]
        cat(sprintf("  %s Day-10 post g = %+.2f -> larger reduction than %.0f%% of %d RCT effects\n",
                    cn, g, 100 * mean(rct_es > g), length(rct_es)))
    }
} else {
    message("[reproducibility] Lasgaard Table S1 xlsx not found; skipping percentile comparison.")
}
