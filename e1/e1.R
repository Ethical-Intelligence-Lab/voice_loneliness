## clear workspace
rm(list = ls())

options(download.file.method = "libcurl")

if (!require(pacman)) { install.packages("pacman") }
pacman::p_load('ggplot2',
               'effsize',
               'tidyr',
               'lmerTest',
               'MatchIt',
               'cobalt',
               'interactions',
               'exact2x2',
               'reshape2',
               'pwr',
               'ltm',
               'sjstats',
               'dplyr',
               'quanteda',
               'topicmodels',
               'tidytext',
               'stopwords',
               'parallel',
               'stringr',
               'irr',
               'BayesFactor')

# Guarded so the script also runs headlessly (Rscript), where rstudioapi is unavailable;
# in that case the working directory is assumed to already be this file's directory.
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) #set working directory to current directory
}

# ---- Nature Research Journals figure settings (per NRJ "Guide to preparing final artwork") ----
# - Sans-serif typeface (Helvetica/Arial); all text between 5 pt and 7 pt.
# - Single-column = 88 mm, double-column = 180 mm; figures saved as vector PDF in RGB.
# Plots below overlay individual participant data points on every summary figure where feasible.
NATURE_FONT    <- "Helvetica"
NATURE_1COL_MM <- 88
NATURE_2COL_MM <- 180

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = NATURE_FONT) +
    theme(
      text            = element_text(size = base_size, family = NATURE_FONT),
      axis.text       = element_text(size = base_size - 1, colour = "black"),
      axis.title      = element_text(size = base_size),
      plot.title      = element_text(size = base_size, hjust = 0.5),
      legend.text     = element_text(size = base_size - 1),
      legend.title    = element_text(size = base_size),
      legend.position = "top",
      axis.line       = element_line(linewidth = 0.25),
      axis.ticks      = element_line(linewidth = 0.25)
    )
}

# Okabe-Ito colourblind-safe palette (distinguishable under red-green CVD), keyed by condition.
okabe_ito_conditions <- c(
  do_nothing          = "#999999",  # grey
  both_text           = "#0072B2",  # blue
  both_voice          = "#D55E00",  # vermillion
  ai_voice_human_text = "#009E73",  # bluish green
  ai_text_human_voice = "#CC79A7",  # reddish purple
  choose_your_own     = "#E69F00"   # orange
)

# plots/ is not committed (gitignored), so create it before the first ggsave() below.
dir.create("plots", showWarnings = FALSE)

d <- read.csv("data.csv")

# The first two rows of the Qualtrics export are metadata (question wording and
# ImportId JSON), not participants. Drop them explicitly so the session count below
# refers to people (the "session records collected" figure in LAB_LOG.md).
d <- d[-(1:2), ]
print(paste0("Session records collected: ", nrow(d)))

d <- d[d$Finished == "1", ]

############################ ATTENTION AND COMPREHENSION CHECKS ############################

# Attention check
size_before <- dim(d)[1]
d <- d[(d$att1 == "2" & d$att2 == "2"),]
print(paste0("Number of participants hired (before loneliness screening): ", dim(d)[1]))

d_screened_out <- d[d$screenedOutLowLoneliness == "true",]

d <- d[d$screenedOutLowLoneliness != "true",]

print(paste0("Number of participants hired (after loneliness screening): ", dim(d)[1]))

size_before <- dim(d)[1]

table(d$condition)

# Exclude participants based on comprehension check for all 3 conditions: 'Control', 'Prediction', 'Experience'
size_before <- dim(d)[1]
d <- d[d$comp_1 == "1",]
print(paste0("Exclusions from comprehension check: ", size_before - dim(d)[1]))

size_before <- dim(d)[1]
d <- d[d$ableToDoNothing != "2",]
print(paste0("Exclusions from ableToDoNothing: ", size_before - dim(d)[1]))

print(paste0("Number of participants after exclusions: ", dim(d)[1]))

table(d$condition)

######### Descriptive Stats #########

print(paste0("Age: ", mean(as.numeric(d$age), trim = 0, na.rm = TRUE))) ## mean age
print(paste0("Percentage of females: ", dim(d[d$gender == "2",])[1] / dim(d)[1]))
print(paste0("AI Experience: ", dim(d[d$ai_companion_exp == "1",])[1] / dim(d)[1]))

# Get number of participants in both experience and control conditions
print(table(d$condition))

conditions <- c("both_text", "both_voice", "ai_voice_human_text", "ai_text_human_voice", "choose_your_own", "do_nothing")

######### Descriptive Stats #########

# Convert all columns whose column name contains following strings to numeric: 'lone', 'comp'
cols <- colnames(d)
cols <- cols[grepl("lone", cols) | grepl("comp", cols) | grepl("social", cols) | grepl("close_friends", cols) | grepl("ai_familiarity", cols) | grepl("ai_capability", cols) | grepl("ai_companion_exp", cols)]
d[cols] <- lapply(d[cols], as.numeric)

# Some participants have these columns completed: "lonely_after_interac_1, lonely_after_interac_2, lonely_after_interac_3" and some have these columns completed: "lonely_after_nothing_1, lonely_after_nothing_2, lonely_after_nothing_3"
# Conbine these to lonely_after_1, lonely_after_2, lonely_after_3
d$lonely_after_1 <- ifelse(is.na(d$lonely_after_interac_1), d$lonely_after_nothing_1, d$lonely_after_interac_1)
d$lonely_after_2 <- ifelse(is.na(d$lonely_after_interac_2), d$lonely_after_nothing_2, d$lonely_after_interac_2)
d$lonely_after_3 <- ifelse(is.na(d$lonely_after_interac_3), d$lonely_after_nothing_3, d$lonely_after_interac_3)

# Now aggregate all 3 loneliness ratings into 1, by taking the mean
d$lonely_before <- rowMeans(d[, c(
    "lonely_before_1", 
    "lonely_before_2", 
    "lonely_before_3")])

d$lonely_after <- rowMeans(d[, c(
    "lonely_after_1", 
    "lonely_after_2", 
    "lonely_after_3")])

# Also rowMeans for 17 item social_anxiety scale, i.e., social_anxiety_1, social_anxiety_2, ..., social_anxiety_17
d$social_anxiety <- rowMeans(d[, c(
    "social_anxiety_1", 
    "social_anxiety_2", 
    "social_anxiety_3",
    "social_anxiety_4",
    "social_anxiety_5",
    "social_anxiety_6",
    "social_anxiety_7",
    "social_anxiety_8",
    "social_anxiety_9",
    "social_anxiety_10",
    "social_anxiety_11",
    "social_anxiety_12",
    "social_anxiety_13",
    "social_anxiety_14",
    "social_anxiety_15",
    "social_anxiety_16",
    "social_anxiety_17")])

# Cronbach's alpha for loneliness ratings
cronbach.alpha(d[, c("lonely_before_1", "lonely_before_2", "lonely_before_3")])
cronbach.alpha(d[, c("lonely_after_1", "lonely_after_2", "lonely_after_3")])
cronbach.alpha(d[, c("social_anxiety_1", 
                     "social_anxiety_2", 
                     "social_anxiety_3",
                     "social_anxiety_4",
                     "social_anxiety_5",
                     "social_anxiety_6",
                     "social_anxiety_7",
                     "social_anxiety_8",
                     "social_anxiety_9",
                     "social_anxiety_10",
                     "social_anxiety_11",
                     "social_anxiety_12",
                     "social_anxiety_13",
                     "social_anxiety_14",
                     "social_anxiety_15",
                     "social_anxiety_16",
                     "social_anxiety_17")])


###################### RESULTS ######################

d$loneliness_reduction <- d$lonely_before - d$lonely_after

# First, we will begin with a one-way ANOVA including all six conditions 
# to assess overall group differences in loneliness reduction.
# Conditions are: both_text, both_voice, do_nothing, ai_voice_human_text, ai_text_human_voice, choose_your_own
# Note: The 'do_nothing' condition is the control condition, while the
# 'choose_your_own' condition is the participant-choice condition.

anova_results <- aov(loneliness_reduction ~ condition, data = d)
summary(anova_results)
anova_stats(anova_results)

# Bayes factor for the omnibus effect of condition (JZS default prior).
# anovaBF() needs condition as an explicit factor and cannot take NAs, so we
# subset to complete cases on the two model variables first.
d_bf_omnibus <- d[!is.na(d$loneliness_reduction) & !is.na(d$condition),
                  c("loneliness_reduction", "condition")]
d_bf_omnibus$condition <- factor(d_bf_omnibus$condition)
anova_bf <- anovaBF(loneliness_reduction ~ condition, data = d_bf_omnibus)
print(anova_bf)


# Additionally,  we will run a 2 (AI modality: voice vs. text)
# × 2 (human modality: voice vs. text) ANOVA (excluding the control and 
# participant-choice conditions) to test for main effects and interaction 
# effects between AI and user communication modality.
d$ai_modality <- ifelse(d$condition %in% c("both_text", "ai_text_human_voice"), "Text", "Voice")
d$human_modality <- ifelse(d$condition %in% c("both_text", "ai_voice_human_text"), "Text", "Voice")

d_no_control_participant_choice <- d[!d$condition %in% c("do_nothing", "choose_your_own"),]

anova_modality <- aov(loneliness_reduction ~ ai_modality * human_modality, data = d_no_control_participant_choice)
summary(anova_modality)
anova_stats(anova_modality)

# Bayes factors for the 2x2 (Table 2a): each main effect is its main-only model
# vs. the intercept; the interaction is the full model divided by the both-main-
# effects model. The interaction BF is Monte Carlo-estimated, so we scope the
# seed locally (via .Random.seed save/restore) to keep it reproducible without
# perturbing the global RNG stream used later by the LDA topic model.
local({
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (!is.null(old_seed))
    assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  set.seed(1)
  dd <- d_no_control_participant_choice
  dd$ai_modality    <- factor(dd$ai_modality)
  dd$human_modality <- factor(dd$human_modality)
  dd <- dd[!is.na(dd$loneliness_reduction), ]
  bf_modality <- anovaBF(loneliness_reduction ~ ai_modality * human_modality, data = dd)
  print(bf_modality)  # [1] ai_modality, [2] human_modality (each vs. intercept)
  cat("2x2 interaction BF10 =",
      as.numeric(extractBF(bf_modality[4])$bf / extractBF(bf_modality[3])$bf), "\n")
})

# Second, we will run paired t-tests comparing baseline loneliness scores to post-interaction loneliness scores in all conditions.
for (condition in conditions) {
    print(paste0("------- Condition: ", condition, " -------"))
    x <- d[d$condition == condition, 'lonely_before']
    y <- d[d$condition == condition, 'lonely_after']

    vart <- var.test(x, y)
    ttest <- t.test(x, y, paired = TRUE)
    cohen.d_result <- cohen.d(y, x, paired=TRUE)

    print(paste0("Paired t-test for condition: ", condition))
    print(ttest)
    print(paste0("Cohen's d for condition: ", condition))
    print(cohen.d_result)

    # Print means
    print(paste0("Mean loneliness before: ", mean(x)))
    print(paste0("Mean loneliness after: ", mean(y)))

    ttest_bayes <- ttestBF(x, y, paired = TRUE)
    print(paste0("Bayes Factor for condition: ", condition))
    print(ttest_bayes)
}

# Third, we will run two-sided t-tests comparing the loneliness 
# difference (before minus after) between the doing nothing condition 
# and each of the other conditions. 

# Run the t-tests
t_test_results <- list()
for (condition in conditions) {
    if (condition == "do_nothing") {
        next  # Skip the control condition
    }

    print(paste0("------- Condition: ", condition, " -------"))
    x <- d[d$condition == "do_nothing", 'loneliness_reduction']
    y <- d[d$condition == condition, 'loneliness_reduction']

    t_test_results[[condition]] <- t.test(x, y, paired = FALSE)
    cohen.d_result <- cohen.d(y, x, pooled=TRUE)

    print(paste0("T-test for condition: ", condition))
    print(t_test_results[[condition]])
    print(paste0("Cohen's d for condition: ", condition))
    print(cohen.d_result)

    # Bayes Factor
    ttest_bayes <- ttestBF(x, y, paired = FALSE)
    print(paste0("Bayes Factor for condition: ", condition))
    print(ttest_bayes)
}


# Fourth, we will run two t-tests to test whether AI responses in voice (vs. text) 
# lead to a greater reduction in loneliness (H1). The comparisons are:
# •	AI voice-human text vs. AI text-human text.
# •	AI voice-human voice vs. AI text-human voice.
print("------- ai_voice_human_text vs. both_text -------")
x <- d[d$condition == "ai_voice_human_text", 'loneliness_reduction']
y <- d[d$condition == "both_text", 'loneliness_reduction']
vart <- var.test(x, y)
ttest_voice_text <- t.test(x, y, paired = FALSE, var.equal = vart$p.value > 0.05)
print(ttest_voice_text)

cohen.d(x, y)
ttestBF(x, y, paired = FALSE)


# AI voice-human voice vs. AI text-human voice
print("------- both_voice vs. ai_text_human_voice -------")
x <- d[d$condition == "both_voice", 'loneliness_reduction']
y <- d[d$condition == "ai_text_human_voice", 'loneliness_reduction']
vart <- var.test(x, y)
ttest_voice_voice <- t.test(x, y, paired = FALSE, var.equal = vart$p.value > 0.05)
print(ttest_voice_voice)

cohen.d(x, y)
ttestBF(x, y, paired = FALSE)


# Fifth, we will conduct two t-tests to test whether participants responding with voice (vs. text) 
# leads to a greater reduction in loneliness (H2). The comparisons are:
# •	AI text-human voice vs. AI text-human text.
# •	AI voice-human voice vs. AI voice-human text.
print("------- ai_text_human_voice vs. both_text -------")
x <- d[d$condition == "ai_text_human_voice", 'loneliness_reduction']
y <- d[d$condition == "both_text", 'loneliness_reduction']
vart <- var.test(x, y)
ttest_voice_text <- t.test(x, y, paired = FALSE, var.equal = vart$p.value > 0.05)
print(ttest_voice_text)

cohen.d(x, y)
ttestBF(x, y, paired = FALSE)

# AI voice-human voice vs. AI text-human voice
print("------- both_voice vs. ai_voice_human_text -------")
x <- d[d$condition == "both_voice", 'loneliness_reduction']
y <- d[d$condition == "ai_voice_human_text", 'loneliness_reduction']
vart <- var.test(x, y)
ttest_voice_voice <- t.test(x, y, paired = FALSE, var.equal = vart$p.value > 0.05)
print(ttest_voice_voice)

cohen.d(x, y)
ttestBF(x, y, paired = FALSE)

# Sixth, we will run exploratory t-tests comparing the participant-choice 
# condition against the other “modality fixed” conditions to investigate whether 
# allowing users to choose their preferred modality leads to greater loneliness reduction. 
# Specifically, we will compare loneliness reduction in the participant-choice condition 
# against each of the four fixed communication modality conditions.
for (condition in c("both_text", "both_voice", "ai_voice_human_text", "ai_text_human_voice")) {
    print(paste0("------- choose_your_own vs. ", condition, " -------"))
    x <- d[d$condition == "choose_your_own", 'loneliness_reduction']
    y <- d[d$condition == condition, 'loneliness_reduction']

    ttest <- t.test(x, y, paired = FALSE)
    print(ttest)

    cohen.d_result <- cohen.d(y, x, pooled=TRUE)
    print(paste0("Cohen's d for condition: choose_your_own vs. ", condition))
    print(cohen.d_result)

    ttestBF <- ttestBF(x, y, paired = FALSE)
    print(ttestBF)
}


# Next, we will identify the voice condition that led to the greatest reduction in loneliness 
# by selecting the one with the highest decrease. Additionally, we will run an ANOVA with 
# interaction terms to examine if age, AI familiarity, number of close friends, and social anxiety
# moderate the effect of the intervention on loneliness reduction. We examine these factors given 
# that they have been implicated in prior studies on AI adoption and loneliness.

# Calculate the mean loneliness reduction for each condition
mean_loneliness_reduction <- tapply(d$loneliness_reduction, d$condition, mean)
print(mean_loneliness_reduction)

# Identify the condition with the highest mean loneliness reduction
best_condition <- names(which.max(mean_loneliness_reduction))
print(paste0("Condition with highest loneliness reduction: ", best_condition))

# Print condition with the second highest loneliness reduction
second_best_condition <- names(sort(mean_loneliness_reduction, decreasing = TRUE)[2])
print(paste0("Condition with second highest loneliness reduction: ", second_best_condition))

# Now, we will run an ANOVA with interaction terms to examine if age, AI familiarity,
# number of close friends, and social anxiety moderate the effect of the intervention on loneliness reduction.
# age is read from the CSV as character (the Qualtrics export carries header rows), so coerce it to
# numeric here; otherwise aov() would treat it as a categorical factor with one level per age.
d$age <- as.numeric(d$age)
anova_moderation <- aov(loneliness_reduction ~ condition * age + condition * ai_familiarity_1 +
                        condition * close_friends + condition * social_anxiety,
                        data = d)
summary(anova_moderation)
anova_stats(anova_moderation)

# Direction of the significant covariate main effects (age, AI familiarity). Because none of the
# condition x covariate interactions is significant, refit as an additive model and report the
# unstandardized coefficients; these are the b [95% CI] values cited in the main text and in
# "Moderation by individual differences" (SI Section A). Both are positive: older participants and
# participants more familiar with AI show greater loneliness reduction.
lm_covariate_direction <- lm(loneliness_reduction ~ condition + age + ai_familiarity_1 +
                             close_friends + social_anxiety, data = d)
print(summary(lm_covariate_direction))
print(round(confint(lm_covariate_direction)[c("age", "ai_familiarity_1"), ], 2))

# We will also conduct ANCOVAs comparing conditions on loneliness reduction while controlling for individual
# difference variables (social anxiety, number of close friends, AI familiarity), to examine the
# robustness of our findings.
anova_loneliness_reduction <- aov(loneliness_reduction ~ condition + social_anxiety + 
                                  close_friends + ai_familiarity_1, data = d)
summary(anova_loneliness_reduction)
anova_stats(anova_loneliness_reduction)



############# EXPLORATORY ANALYSES ###################

##########################################################
############# Baseline Loneliness Moderation #############
##########################################################

# Do we find again that the magnitude of the reduction is a function of baseline loneliness?
# And does that relationship differ by condition (i.e., is baseline loneliness moderation different per condition)?

# Keep factor coding local to these models (so we don't unintentionally change earlier analyses)
d$condition_f <- factor(
  d$condition,
  levels = c(
    "do_nothing",
    "both_text",
    "both_voice",
    "ai_voice_human_text",
    "ai_text_human_voice",
    "choose_your_own"
  )
)

# Single model: do_nothing as reference; test whether any condition differs in the
# baseline-loneliness -> loneliness_reduction slope (i.e., condition interactions).
d$condition_f <- stats::relevel(d$condition_f, ref = "do_nothing")

lm_change_by_condition <- lm(loneliness_reduction ~ lonely_before * condition_f, data = d)
print(summary(lm_change_by_condition))

# Visualize the effect between lonely_before and loneliness reduction by condition (with color)
# Individual participant points are already shown (geom_point); means are summarised by the fitted lines.
plt_l <- ggplot(d, aes(x = lonely_before, y = loneliness_reduction, color = condition_f)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.5) +
  theme_nature() +
  ggtitle("Loneliness Reduction by Baseline Loneliness and Condition") +
  xlab("Baseline Loneliness") +
  ylab("Loneliness Reduction (Pre - Post)") +
  scale_color_manual(name = "Condition", values = okabe_ito_conditions)


plt_l

# save plt_l (double-column width, vector PDF)
ggsave("plots/loneliness_reduction_by_baseline_and_condition.pdf", plt_l,
       width = NATURE_2COL_MM, height = 120, units = "mm")

############################################################################
################### Modality change analysis ###############################
############################################################################

# NOTE: reads the raw app-database exports in database/, which are not committed to
# the public repo (participant identifiers). Guarded so the script runs from committed
# data; this modality analysis only executes when the database/ files are present.
if (file.exists("database/mturk_app_modalitychangeactions.csv") &&
    file.exists("database/mturk_app_participant.csv")) {

modality_change <- read.csv("database/mturk_app_modalitychangeactions.csv", stringsAsFactors = FALSE)
participants <- read.csv("database/mturk_app_participant.csv", stringsAsFactors = FALSE)

# Ensure expected types/labels
modality_change$rem_time <- suppressWarnings(as.numeric(modality_change$rem_time))
modality_change$action <- tolower(trimws(modality_change$action))

participants_filtered <- participants %>%
  dplyr::filter(worker_id %in% d[d$condition == 'choose_your_own', "worker_id"])

modality_change_filtered <- modality_change %>%
  dplyr::filter(user_id %in% participants_filtered$id)

# Time in each modality based on rem_time.
# Interpretation: each row indicates the modality selected at that rem_time (remaining seconds).
# For each user, we sort by decreasing rem_time and allocate time until the next switch (or 0).
modality_durations <- modality_change_filtered %>%
  dplyr::filter(!is.na(rem_time), rem_time >= 0, action %in% c("text", "voice")) %>%
  dplyr::arrange(user_id, dplyr::desc(rem_time), timestamp) %>%
  dplyr::group_by(user_id, rem_time) %>%
  # If there are multiple modality changes logged at the same rem_time,
  # keep only the latest one (by timestamp).
  dplyr::slice_tail(n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(user_id) %>%
  dplyr::arrange(dplyr::desc(rem_time), .by_group = TRUE) %>%
  dplyr::mutate(
    next_rem_time = dplyr::lead(rem_time, default = 0),
    seconds_in_segment = pmax(rem_time - next_rem_time, 0)
  ) %>%
  dplyr::summarise(
    seconds_text = sum(seconds_in_segment[action == "text"], na.rm = TRUE),
    seconds_voice = sum(seconds_in_segment[action == "voice"], na.rm = TRUE),
    total_seconds_observed = max(rem_time, na.rm = TRUE),
    num_modality_actions = dplyr::n(),
    num_modality_switches = pmax(dplyr::n() - 1L, 0L),
    .groups = "drop"
  )

# Add worker_id for merging back to d
user_worker <- participants_filtered %>%
  dplyr::select(user_id = id, worker_id)

modality_durations <- merge(modality_durations, user_worker, by = "user_id", all.x = TRUE)

# Merge modality change counts back to main data frame (keep original d intact)
d_cyo <- d[d$condition == 'choose_your_own', ]
d_cyo <- merge(d_cyo, modality_durations, by = "worker_id", all.x = TRUE)

mean(d_cyo$seconds_text, na.rm = TRUE)
mean(d_cyo$seconds_voice, na.rm = TRUE)
mean(d_cyo$num_modality_switches, na.rm = TRUE)
sd(d_cyo$num_modality_switches, na.rm = TRUE)

cor.test(d_cyo$num_modality_switches, d_cyo$loneliness_reduction)

# Optional: relate time-in-modality to loneliness reduction
lm_modality_time <- lm(loneliness_reduction ~ seconds_text, data = d_cyo)
summary(lm_modality_time)

# get people who spend more time in voice than text
d_cyo$prefers_voice <- ifelse(d_cyo$seconds_voice > d_cyo$seconds_text, 1, 0)
table(d_cyo$prefers_voice)

mean(d_cyo[d_cyo$prefers_voice == 1, c("seconds_voice")])
mean(d_cyo[d_cyo$prefers_voice == 0, c("seconds_text")])

# Compare loneliness reduction between those who prefer voice vs. text
t_test_prefers_voice <- t.test(
  d_cyo$loneliness_reduction[d_cyo$prefers_voice == 1],
  d_cyo$loneliness_reduction[d_cyo$prefers_voice == 0],
  paired = FALSE
)
print(t_test_prefers_voice)


# Compare social anxiety between those who prefer voice vs. text
t_test_social_anxiety <- t.test(
  d_cyo$social_anxiety[d_cyo$prefers_voice == 1],
  d_cyo$social_anxiety[d_cyo$prefers_voice == 0],
  paired = FALSE
)
print(t_test_social_anxiety)

# Compare number of close friends between those who prefer voice vs. text
t_test_close_friends <- t.test(
  d_cyo$close_friends[d_cyo$prefers_voice == 1],
  d_cyo$close_friends[d_cyo$prefers_voice == 0],
  paired = FALSE
)
print(t_test_close_friends)

# Compare AI familiarity between those who prefer voice vs. text
t_test_ai_familiarity <- t.test(
  d_cyo$ai_familiarity_1[d_cyo$prefers_voice == 1],
  d_cyo$ai_familiarity_1[d_cyo$prefers_voice == 0],
  paired = FALSE
)
print(t_test_ai_familiarity)

d_cyo$age <- as.numeric(d_cyo$age)
# compare age between those who prefer voice vs. text
t_test_age <- t.test(
  d_cyo$age[d_cyo$prefers_voice == 1],
  d_cyo$age[d_cyo$prefers_voice == 0],
  paired = FALSE
)
print(t_test_age)

cohen.d_age <- cohen.d(
  d_cyo$age[d_cyo$prefers_voice == 1],
  d_cyo$age[d_cyo$prefers_voice == 0],
  pooled = TRUE
)

cohen.d_age



# Compare d_cyo$prefers_voice == 1 with both_voice condition
t.test_prefers_voice_vs_both_voice <- t.test(
  d_cyo$loneliness_reduction[d_cyo$prefers_voice == 1],
  d$loneliness_reduction[d$condition == 'both_voice'],
  paired = FALSE
)
print(t.test_prefers_voice_vs_both_voice)

cohen.d_prefers_voice_vs_both_voice <- cohen.d(
  d_cyo$loneliness_reduction[d_cyo$prefers_voice == 1],
  d$loneliness_reduction[d$condition == 'both_voice'],
  pooled = TRUE
)
cohen.d_prefers_voice_vs_both_voice

# Compare d_cyo$prefers_voice == 0 with both_text condition
t.test_prefers_text_vs_both_text <- t.test(
  d_cyo$loneliness_reduction[d_cyo$prefers_voice == 0],
  d$loneliness_reduction[d$condition == 'both_text'],
  paired = FALSE
)
print(t.test_prefers_text_vs_both_text)

cohen.d_prefers_text_vs_both_text <- cohen.d(
  d_cyo$loneliness_reduction[d_cyo$prefers_voice == 0],
  d$loneliness_reduction[d$condition == 'both_text'],
  pooled = TRUE
)
cohen.d_prefers_text_vs_both_text

} else {
  message("[reproducibility] Skipping modality-change analysis: database/ exports not present (excluded from the public repo).")
}

###################### PLOTTING ######################

# Create a long format of the data, i.e., before and after loneliness ratings in separate rows
d_long <- melt(d, id.vars = c('condition', 'worker_id'), measure.vars = c('lonely_before', 'lonely_after'))

# Rename lonely_before as before and lonely_after as after.
d_long$variable <- ifelse(d_long$variable == "lonely_before", "Pre", "Post")
d_long$variable <- factor(d_long$variable, levels = c("Pre", "Post"))
d_long$condition <- factor(d_long$condition, levels = c('do_nothing', 'both_text', 'both_voice', 'ai_voice_human_text', 'ai_text_human_voice', 'choose_your_own'))

# Create more readable condition labels
d_long$condition_label <- factor(d_long$condition, 
                                levels = c('do_nothing', 'both_text', 'both_voice', 'ai_voice_human_text', 'ai_text_human_voice', 'choose_your_own'),
                                labels = c('Control\n(Do Nothing)', 'Both Text', 'Both Voice', 'AI Voice\nHuman Text', 'AI Text\nHuman Voice', 'Participant\nChoice'))

# Then change 'variable' column name to 'lonely'
d_long <- d_long %>% rename(lonely = variable)

plotter <- function(y_var='value', title='') {
  p1 <- ggplot(d_long, aes_string(x = 'condition_label', y = y_var, fill = 'lonely')) +
    geom_bar(position = position_dodge(width = 0.7), stat = "summary", fun = mean,
             width = 0.7, alpha = 1, colour = "black", linewidth = 0.3) +
    # Individual participant ratings overlaid on the condition means (jittered within each dodged bar)
    geom_point(
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
      size = 0.25, alpha = 0.15, shape = 16, colour = "black", show.legend = FALSE
    ) +
    stat_summary(
      fun.data = mean_cl_boot,
      geom = "errorbar",
      color = "black",
      width = 0.15,
      linewidth = 0.4,
      position = position_dodge(width = 0.7)
    ) +
    theme_nature() +
    ggtitle(title) +
    scale_fill_manual(values = c("#CCCCCC", "#404040"),
                     name = "Measurement Time:",
                     guide = guide_legend(reverse = FALSE)) +
    xlab("Condition") +
    ylab("Loneliness Rating") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

  return(p1)
}

p1 <- plotter()

# Create plots directory if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

ggsave("plots/loneliness_by_condition.pdf", p1, width = NATURE_2COL_MM, height = 120, units = "mm")
ggsave("plots/loneliness_by_condition.png", p1, dpi = 300, width = NATURE_2COL_MM, height = 120, units = "mm")


##########################################################
################# Engagement Analysis ####################
##########################################################

# NOTE: reads the raw message-level app-database export (database/), not committed to
# the public repo. Guarded so the script runs from committed data; the engagement
# metrics (num_messages, avg_words_per_message) only build when database/ is present.
if (file.exists("database/mturk_app_message.csv") &&
    file.exists("database/mturk_app_participant.csv")) {

messages <- read.csv("database/mturk_app_message.csv", stringsAsFactors = FALSE)
participants <- read.csv("database/mturk_app_participant.csv", stringsAsFactors = FALSE)

participants_filtered <- participants %>%
  dplyr::filter(worker_id %in% d$worker_id)

# Merge messages with participants by participants$id and messages$user_id
merged_data <- merge(messages, participants_filtered, by.x = "user_id", by.y = "id")

# Robust word count (no extra packages)
merged_data$content_clean <- trimws(ifelse(is.na(merged_data$content), "", merged_data$content))
merged_data$word_count <- ifelse(
  merged_data$content_clean == "",
  0L,
  lengths(strsplit(merged_data$content_clean, "\\s+"))
)

# Calculate engagement metrics per participant
engagement_metrics <- merged_data %>%
  dplyr::group_by(worker_id) %>%
  dplyr::summarise(
    num_messages = dplyr::n(),
    avg_words_per_message = mean(word_count, na.rm = TRUE),
    .groups = "drop"
  )

# Merge engagement metrics back to main data frame (keep original d intact)
d_eng <- merge(d, engagement_metrics, by = "worker_id", all.x = TRUE)

# Participants with no logged messages will have NA engagement metrics after the merge.
# Treat these as 0 for engagement.
d_eng$num_messages[is.na(d_eng$num_messages)] <- 0
d_eng$avg_words_per_message[is.na(d_eng$avg_words_per_message)] <- 0

# Predict loneliness reduction from engagement metrics
# (exclude do_nothing: no chatbot interaction by design)
d_eng_chatbot <- d_eng[d_eng$condition != 'do_nothing', ]

# Factor for condition (reference = both_text by default; change here if you want another baseline)
d_eng_chatbot$condition_f <- factor(
  d_eng_chatbot$condition,
  levels = c(
    "both_text",
    "both_voice",
    "ai_voice_human_text",
    "ai_text_human_voice",
    "choose_your_own"
  )
)

# Center predictors for cleaner interaction interpretation
d_eng_chatbot$num_messages_c <- d_eng_chatbot$num_messages - mean(d_eng_chatbot$num_messages, na.rm = TRUE)
d_eng_chatbot$avg_words_per_message_c <- d_eng_chatbot$avg_words_per_message - mean(d_eng_chatbot$avg_words_per_message, na.rm = TRUE)

# Main-effects model (no interactions)
lm_engagement_main <- lm(
  loneliness_reduction ~ condition_f + num_messages_c + avg_words_per_message_c,
  data = d_eng_chatbot
)
print(summary(lm_engagement_main))

# Print effect size of lm (eta-squared) and coefficient CIs. Both are reported in SI
# Section A; print() is required because this sits inside the database/ if-block, where
# only the block's final expression would otherwise auto-print.
library(lsr)
print(etaSquared(lm_engagement_main))
print(round(confint(lm_engagement_main), 4))

# Interaction model: do engagement effects differ by condition?
lm_engagement_interaction <- lm(
  loneliness_reduction ~ condition_f * (num_messages_c + avg_words_per_message_c),
  data = d_eng_chatbot
)
print(summary(lm_engagement_interaction))


# Anova of engagement and conditino
anova_engagement <- aov(
  loneliness_reduction ~ condition_f * (num_messages_c + avg_words_per_message_c),
  data = d_eng_chatbot
)

summary(anova_engagement)

# PLot condition_fboth_voice:avg_words_per_message_c 

# plot.points = TRUE overlays the individual participant observations on the fitted lines.
plt <- interact_plot(
  lm_engagement_interaction,
  pred = avg_words_per_message_c,
  modx = condition_f,
  colors = unname(okabe_ito_conditions[levels(d_eng_chatbot$condition_f)]),
  plot.points = TRUE,
  point.size = 0.6,
  interval = TRUE,
  int.type = "confidence",
  main.title = "Interaction of Avg. Words per Message and Condition on Loneliness Reduction",
  x.label = "Average Words per Message (centered)",
  y.label = "Loneliness Reduction (Pre - Post)"
) + theme_nature()

plt

# Save plt (double-column width, vector PDF)
ggsave("plots/interaction_avg_words_per_message_condition.pdf", plt,
       width = NATURE_2COL_MM, height = 120, units = "mm")


# Plot engagement as two separate plots, saved into one PDF figure
p3_messages <- ggplot(d_eng, aes(x = condition, y = num_messages)) +
  stat_summary(fun = mean, fun.args = list(na.rm = TRUE), geom = "bar", fill = "grey70", alpha = 0.9) +
  # Individual participants overlaid on the condition means
  geom_jitter(width = 0.15, height = 0, size = 0.3, alpha = 0.2, colour = "grey20") +
  stat_summary(
    fun.data = mean_cl_boot,
    geom = "errorbar",
    color = "black",
    width = 0.2,
    linewidth = 0.4
  ) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Average # Messages") +
  xlab("Condition") +
  ylab("Average Number of Messages") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

p3_words <- ggplot(d_eng, aes(x = condition, y = avg_words_per_message)) +
  stat_summary(fun = mean, fun.args = list(na.rm = TRUE), geom = "bar", fill = "grey70", alpha = 0.9) +
  # Individual participants overlaid on the condition means
  geom_jitter(width = 0.15, height = 0, size = 0.3, alpha = 0.2, colour = "grey20") +
  stat_summary(
    fun.data = mean_cl_boot,
    geom = "errorbar",
    color = "black",
    width = 0.2,
    linewidth = 0.4
  ) +
  theme_nature() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Average Words / Message") +
  xlab("Condition") +
  ylab("Average Words per Message") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Display both plots in the current device (no extra packages)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
print(p3_messages, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p3_words, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))

# Save both plots into a single PDF
if (!dir.exists("plots")) {
  dir.create("plots")
}

# 180 mm (double-column) wide, two panels stacked; pdf() expects inches (180 mm = 7.09 in).
grDevices::pdf("plots/engagement_by_condition.pdf", width = NATURE_2COL_MM / 25.4, height = 8.5)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 1)))
print(p3_messages, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p3_words, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
grDevices::dev.off()

} else {
  message("[reproducibility] Skipping engagement analysis: database/ message exports not present (excluded from the public repo).")
}


###################################################

# Finally, we will analyze the conversation data to assess the richness and depth of the chatbot interactions. 
# First, we will perform topic modeling on the conversation transcripts (both text and voice transcribed to text) 
# using Latent Dirichlet Allocation. This will allow us to identify dominant themes across participants and 
# assess the variety of conversation topics. 

# NOTE: rebuilds conversation transcripts from the raw app-database exports (database/),
# not committed to the public repo. The committed data_with_transcripts.csv already holds
# the output of this block. Guarded so the script runs from committed data; transcript
# rebuild + topic modelling only execute when database/ is present.
if (file.exists("database/mturk_app_participant.csv") &&
    file.exists("database/mturk_app_message.csv")) {

participants <- read.csv("database/mturk_app_participant.csv", stringsAsFactors = FALSE)
messages <- read.csv("database/mturk_app_message.csv", stringsAsFactors = FALSE)

############################
# MERGE PARTICIPANTS WITH d
############################

d <- d %>%
  left_join(participants, by = "worker_id")

d$condition <- d$condition.x

############################
# BUILD AI TRANSCRIPTS
############################

# Keep only messages belonging to participants in d
messages_filtered <- messages %>%
  filter(room_id %in% d$room_id)

# Ensure correct order
messages_filtered <- messages_filtered %>%
  arrange(room_id, datetime)

# Convert sender type to readable labels
messages_filtered <- messages_filtered %>%
  mutate(
    speaker = case_when(
      sender_type == "human" ~ "User",
      sender_type == "chatbot" ~ "Chatbot",
      TRUE ~ "Chatbot"
    ),
    formatted = paste0(speaker, ": ", content)
  )

# Collapse messages into one transcript per room
transcripts_all <- messages_filtered %>%
  group_by(room_id) %>%
  summarise(
    whole_text_transcript = paste(formatted, collapse = "\n"),
    .groups = "drop"
  )

# Only get human messages for the user text analysis
transcripts_user <- messages_filtered %>%
  filter(speaker == "User") %>%
  group_by(room_id) %>%
  summarise(
    user_text_transcript = paste(content, collapse = "\n"),
    .groups = "drop"
  )

############################
# MERGE TRANSCRIPTS INTO d
############################

d <- d %>%
  left_join(transcripts_all, by = "room_id") %>%
  left_join(transcripts_user, by = "room_id")

# Save file
write.csv(d, "data_with_transcripts.csv", row.names = FALSE)

### Topic modelling using transcripts_user

# Create plots directory if it doesn't exist
if (!dir.exists("plots")) {
  dir.create("plots")
}

# Clean text function
clean_text <- function(text) {
  text <- tolower(text)
  text <- gsub(
    "http[s]?://(?:[a-zA-Z]|[0-9]|[$-_@.&+]|[!*\\(\\),]|(?:%[0-9a-fA-F][0-9a-fA-F]))+",
    "",
    text
  )
  text <- gsub("\\S+@\\S+", "", text)
  text <- gsub("\\s+", " ", text)
  text <- trimws(text)
  return(text)
}

# Topic modelling function (LDA over a chosen text column)
topic_modelling <- function(
  data_topic,
  text_column = "user_text_transcript",
  check_perplexity = FALSE,
  condition = "all",
  k_topics = 10L,
  top_n_terms = 10L
) {
  explanations <- data_topic[[text_column]]
  cat("Original number of documents:", length(explanations), "\n")

  explanations <- clean_text(explanations)
  cat("After cleaning:", length(explanations), "\n")

  explanations <- explanations[!is.na(explanations) & explanations != ""]
  cat("After removing NAs/empty:", length(explanations), "\n")

  if (length(explanations) == 0) {
    cat("Error: No valid documents after preprocessing\n")
    return(NULL)
  }

  texts <- quanteda::corpus(explanations)
  toks <- quanteda::tokens(texts, remove_punct = TRUE, remove_numbers = TRUE)
  dfm_obj <- quanteda::dfm(toks)
  dfm_obj <- quanteda::dfm_remove(dfm_obj, stopwords::stopwords(source = "smart"))

  # Remove documents with no terms
  dfm_obj <- dfm_obj[rowSums(dfm_obj) > 0, ]
  cat("DFM dimensions after removing empty documents:", dim(dfm_obj), "\n")

  dtm <- quanteda::convert(dfm_obj, to = "topicmodels")

  if (nrow(dtm) == 0 || ncol(dtm) == 0) {
    cat("Error: Document-term matrix is empty after preprocessing\n")
    return(NULL)
  }

  # Optional perplexity sweep (kept for parity with the example; default off)
  if (isTRUE(check_perplexity)) {
    set.seed(123)
    train <- sample(rownames(dtm), floor(nrow(dtm) * 0.8))
    dtm_train <- dtm[rownames(dtm) %in% train, ]
    dtm_test <- dtm[!rownames(dtm) %in% train, ]

    p <- data.frame(k = 2:20, perplexity = NA_real_)

    calc_perplexity <- function(k) {
      perplexity_scores <- rep(NA_real_, 10)
      for (j in 1:10) {
        cat(paste0("*** ", k, ": ", j, " ***\n"))
        m <- topicmodels::LDA(dtm_train, method = "Gibbs", k = k, control = list(seed = 123))
        perplexity_scores[j] <- topicmodels::perplexity(m, dtm_test)
      }
      mean(perplexity_scores)
    }

    p$perplexity <- unlist(parallel::mclapply(p$k, calc_perplexity, mc.cores = parallel::detectCores()))
    p$k <- as.numeric(p$k)
    p$perplexity <- as.numeric(p$perplexity)

    plt <- ggplot(p, aes(x = k, y = perplexity)) +
      geom_point(size = 0.6) +
      geom_line(linewidth = 0.4) +
      theme_minimal(base_size = 7, base_family = NATURE_FONT) +
      theme(text = element_text(size = 7, family = NATURE_FONT)) +
      labs(x = "Number of Topics", y = "Perplexity") +
      theme(axis.line = element_line(colour = "black", linewidth = 0.25))

    print(plt)
    sstr <- gsub(" ", "_", condition)
    ggsave(paste0("plots/perplexity_", sstr, ".pdf"), plt, dpi = 300, width = 8, height = 6)
  }

  # Fit LDA
  set.seed(123)
  model <- topicmodels::LDA(dtm, method = "Gibbs", k = as.integer(k_topics), control = list(seed = 123))

  # Extract top terms per topic (for plotting + console output)
  topics <- tidytext::tidy(model, matrix = "beta") %>%
    dplyr::mutate(topic = as.integer(topic))
  top_terms <- topics %>%
    dplyr::group_by(topic) %>%
    dplyr::slice_max(beta, n = as.integer(top_n_terms)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(topic, -beta)

  fig <- top_terms %>%
    dplyr::mutate(
      topic_facet = factor(
        topic,
        levels = sort(unique(topic)),
        labels = paste("Topic", sort(unique(topic)))
      ),
      term = tidytext::reorder_within(term, beta, topic_facet)
    ) %>%
    ggplot(aes(beta, term)) +
    theme(text = element_text(size = 7, family = NATURE_FONT)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.25),
      axis.text = element_text(size = 6, family = NATURE_FONT),
      strip.text = element_text(size = 7, family = NATURE_FONT),
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5)
    ) +
    xlab("Topic Word Probability") +
    ylab("") +
    geom_col(show.legend = FALSE, fill = "#4f4f4f") +
    facet_wrap(~ topic_facet, scales = "free", ncol = 3) +
    tidytext::scale_y_reordered()

  print(fig)

  cat(paste0("Condition: ", condition, "\n"))
  for (i in seq_len(as.integer(k_topics))) {
    cat(paste0("Top terms for Topic ", i, ":\n"))
    top_terms_topic <- top_terms %>%
      dplyr::filter(topic == i) %>%
      dplyr::arrange(-beta) %>%
      dplyr::slice(1:as.integer(top_n_terms))
    print(top_terms_topic)
    cat("\n")
  }

  sstr <- gsub(" ", "_", condition)
  ggsave(paste0("plots/topics_", sstr, ".pdf"), fig, dpi = 300, width = 15, height = 10)

  return(model)
}

# Run topic modelling once with all conditions combined
d_topic <- d %>%
  dplyr::filter(!is.na(user_text_transcript) & user_text_transcript != "")

lda_model_all <- topic_modelling(
  data_topic = d_topic,
  text_column = "user_text_transcript",
  check_perplexity = FALSE,
  condition = "all_conditions",
  k_topics = 15L
)

} else {
  message("[reproducibility] Skipping transcript rebuild + topic modelling: database/ exports not present (excluded from the public repo).")
}

# Finally, to gain a sense of the emotional depth of the conversations, we will randomly select 200 participant 
# conversations and two of the authors will manually code them for emotional depth, using a binary coding scheme 
# developed for the current study.

if(FALSE) {
  data_with_llm <- read.csv("data_with_llm_classifications.csv")
  table(data_with_llm$label)

  # Randomly sample 200 conversations for manual coding
  set.seed(123)
  sampled_conversations <- data_with_llm %>%
    filter(!is.na(user_text_transcript) & user_text_transcript != "") %>%
    sample_n(size = min(200, nrow(.)))

  # Save sampled conversations for manual coding
  write.csv(sampled_conversations, "sampled_conversations_for_manual_coding_new.csv")
}

# Now read sampled_conversations_for_manual_coding.csv
manual_coding <- read.csv("sampled_conversations_for_manual_coding.csv", sep=";")
colnames(manual_coding)

# Get cohens kappa between classification_k and classification_z
kappa_result <- irr::kappa2(manual_coding[, c("classification_k", "classification_z")])
print(kappa_result)

# Now get the number of conversations coded as emotionally deep by each coder
num_deep_k <- sum(manual_coding$classification_k == 1, na.rm = TRUE)
num_deep_z <- sum(manual_coding$classification_z == 1, na.rm = TRUE)
n_coded <- sum(!is.na(manual_coding$classification_k))
cat(sprintf("Emotional content present (coder k): %d of %d conversations (%.1f%%)\n",
            num_deep_k, n_coded, 100 * num_deep_k / n_coded))
cat(sprintf("Emotional content present (coder z): %d of %d conversations (%.1f%%)\n",
            num_deep_z, n_coded, 100 * num_deep_z / n_coded))

# The emotional-depth regression combines the manual codes with engagement metrics
# (num_messages_c, avg_words_per_message_c) built in the engagement block from the
# database/ exports. The Cohen's kappa above runs from committed data; this regression
# only runs when d_eng_chatbot is available (i.e. database/ was present).
if (exists("d_eng_chatbot")) {

# Add num messages and avg words per message to manual_coding by merging with d_eng_chatbot
manual_coding <- manual_coding %>%
  left_join(d_eng_chatbot %>% select(worker_id, num_messages_c, avg_words_per_message_c, loneliness_reduction), by = "worker_id")


# We included this measure in a linear regression alongside engagement metrics, predicting loneliness decrease.
lm_emotional_depth <- lm(
  loneliness_reduction.x ~ classification_k + num_messages_c + avg_words_per_message_c,
  data = manual_coding
)

print(summary(lm_emotional_depth))

library(lsr)
# print() is required here: adding the confint() call below means etaSquared() is no
# longer the if-block's final expression and would otherwise stop auto-printing.
print(etaSquared(lm_emotional_depth))
print(round(confint(lm_emotional_depth), 4))

} else {
  message("[reproducibility] Skipping emotional-depth regression: engagement metrics (d_eng_chatbot) require database/ exports; Cohen's kappa above reproduces from committed data.")
}

############################################################################
######### Comparison with Lasgaard et al. (2026) meta-analysis ############
############################################################################
# Manuscript (General Discussion): "the immediate reductions in our studies
# (g = 0.38-0.59) were comparable to the pooled short-term effect ... (g ~ 0.50
# across 122 randomized controlled trials) and exceeded the study-level effects
# reported in the majority (51-66%) of those trials".
#
# Lasgaard et al. (2026, American Psychologist) compute RCT effects as
# standardized mean differences based on Hedges' g, from post-intervention
# means and SDs (intervention vs. control); negative values = loneliness
# reduction. We compute the identical statistic for Study 1: post-interaction
# loneliness in each AI condition vs. the do-nothing control. The small-sample
# correction is J = 1 - 3/(4*df - 1) (Hedges, 1981); CIs use qnorm(.975).
# We then locate each estimate in the distribution of the 122 RCT short-term
# effects parsed from their Table S1 (file provided by the authors' journal
# supplement; guarded so the script still runs without it).

hedges_j <- function(df) 1 - 3 / (4 * df - 1)
hedges_g_between <- function(x1, x0) {
    n1 <- length(x1); n0 <- length(x0); df <- n1 + n0 - 2
    sp <- sqrt(((n1 - 1) * var(x1) + (n0 - 1) * var(x0)) / df)
    g  <- (mean(x1) - mean(x0)) / sp * hedges_j(df)
    se <- sqrt((n1 + n0) / (n1 * n0) + g^2 / (2 * df))
    c(g = g, lo = g - qnorm(.975) * se, hi = g + qnorm(.975) * se)
}

x_ctl <- d$lonely_after[d$condition == "do_nothing"]
g_s1 <- t(sapply(setdiff(unique(d$condition), "do_nothing"), function(cn)
    hedges_g_between(d$lonely_after[d$condition == cn], x_ctl)))
print(round(g_s1, 2))
g_s1_pooled <- hedges_g_between(d$lonely_after[d$condition != "do_nothing"], x_ctl)
cat(sprintf("Study 1, all AI conditions pooled: g = %.2f [%.2f, %.2f]\n",
            g_s1_pooled[1], g_s1_pooled[2], g_s1_pooled[3]))

lasgaard_s1_file <- "../lit/Supplemental Materials Table S1 AMP-2022-0748.xlsx"
if (file.exists(lasgaard_s1_file) && requireNamespace("readxl", quietly = TRUE)) {
    ts1 <- as.data.frame(readxl::read_excel(lasgaard_s1_file, sheet = 1, skip = 1))
    ts1 <- ts1[!is.na(ts1[["Study design"]]), ]
    # effects appear as e.g. "-0.16 [0.25, -0.57]"; some rows use decimal commas
    es <- as.numeric(sub("^\\s*(-?[0-9.]+)\\s*\\[.*", "\\1",
                         gsub(",", ".", ts1[["Short-term effect [95% CI]"]])))
    rct_es <- es[ts1[["Study design"]] == "RCT" & !is.na(es)]
    cat(sprintf("Lasgaard Table S1: %d RCT short-term effects, median = %.2f (pooled in paper: -0.50)\n",
                length(rct_es), median(rct_es)))
    for (i in seq_len(nrow(g_s1)))
        cat(sprintf("  %-22s g = %+.2f -> larger reduction than %.0f%% of RCT effects\n",
                    rownames(g_s1)[i], g_s1[i, 1], 100 * mean(rct_es > g_s1[i, 1])))

    # Basis for the manuscript's qualifier that the meta-analytic estimates index
    # TRAIT loneliness whereas ours index momentary loneliness: tabulate the
    # instrument used by each of those RCTs. All are dispositional scales (UCLA,
    # De Jong Gierveld, LSDS, ...); none uses a momentary/EMA protocol.
    instr <- trimws(ts1[["Loneliness measure"]][ts1[["Study design"]] == "RCT" & !is.na(es)])
    fam <- ifelse(grepl("UCLA", instr, ignore.case = TRUE), "UCLA (trait scale)",
           ifelse(grepl("DJG|Jong", instr, ignore.case = TRUE), "De Jong Gierveld (trait scale)",
           ifelse(grepl("^Single-item$", instr), "Single-item (trait framing)",
                  "Other multi-item scale")))
    cat("Instruments used by these RCTs:\n")
    print(sort(table(fam), decreasing = TRUE))
} else {
    message("[reproducibility] Lasgaard Table S1 xlsx not found; skipping percentile comparison.")
}
