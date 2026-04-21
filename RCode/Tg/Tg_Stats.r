# =============================================================================
#  Tg_Results_NoCoh1 — Spatial Strategy Analysis
#  Days 2 & 4 | Strategy Use, Sex, Genotype (APP/+ vs WT)
# =============================================================================
# Required packages:
#   readxl, dplyr, tidyr, ggplot2, lme4, lmerTest, emmeans, rstatix, car
#
# Install if needed:
#   install.packages(c("readxl","dplyr","tidyr","ggplot2","lme4",
#                      "lmerTest","emmeans","rstatix","car","scales"))
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)   # p-values for lmer
library(emmeans)
library(rstatix)    # pipe-friendly t-tests / ANOVA helpers
library(car)        # Levene's test
library(scales)

# ---------------------------------------------------------------------------
# 0.  LOAD & PREPARE DATA
# ---------------------------------------------------------------------------

raw <- read_excel('/Users/miasponseller/Desktop/Lab/Rtrack/Tg/Tg_Results_NoCoh1.xlsx')

# ---- Recode / factor key variables ----------------------------------------
raw <- raw %>%
  rename(Day      = `_Day`,
         Trial    = `_Trial`,
         SubjectID = `_TargetID`) %>%
  mutate(
    Sex      = factor(Sex, levels = c("M", "F")),
    Genotype = factor(APP,  levels = c("WT", "APP/+")),  # WT as reference
    Day      = as.integer(Day),
    Age      = as.numeric(Age),
    # Strategy as ordered factor (1 = least spatial → 8 = most spatial)
    StrategyNum  = as.integer(strategy),
    StrategyName = factor(name, levels = c("thigmotaxis","circling",
                                           "random path","scanning",
                                           "chaining","directed search",
                                           "corrected path","direct path"))
  )

# ---- Subset to Days 2 & 4 --------------------------------------------------
d24 <- raw %>% filter(Day %in% c(2, 4))

cat("=== Data dimensions ===\n")
cat("Full dataset:", nrow(raw), "rows,", n_distinct(raw$SubjectID), "subjects\n")
cat("Days 2 & 4  :", nrow(d24), "rows\n\n")

# ---------------------------------------------------------------------------
# 1.  DESCRIPTIVE STATISTICS — Age & Sex per Strategy (Days 2 & 4)
# ---------------------------------------------------------------------------

cat("============================================================\n")
cat("SECTION 1: Descriptive Stats — Age & Sex by Strategy\n")
cat("============================================================\n\n")

# Age summary per strategy
age_by_strat <- d24 %>%
  group_by(StrategyName) %>%
  summarise(
    n        = n(),
    Age_Mean = round(mean(Age, na.rm = TRUE), 2),
    Age_SD   = round(sd(Age,   na.rm = TRUE), 2),
    Age_Min  = min(Age, na.rm = TRUE),
    Age_Max  = max(Age, na.rm = TRUE),
    .groups = "drop"
  )
cat("--- Age by Strategy ---\n")
print(age_by_strat)

# Sex counts per strategy
sex_by_strat <- d24 %>%
  count(StrategyName, Sex) %>%
  pivot_wider(names_from = Sex, values_from = n, values_fill = 0)
cat("\n--- Sex counts by Strategy ---\n")
print(sex_by_strat)

# Age summary across ALL strategies combined (Days 2 & 4)
cat("\n--- Age overall (Days 2 & 4) ---\n")
d24 %>%
  summarise(n = n(), Mean = round(mean(Age,na.rm=TRUE),2),
            SD = round(sd(Age,na.rm=TRUE),2),
            Min = min(Age,na.rm=TRUE), Max = max(Age,na.rm=TRUE)) %>%
  print()

# ---------------------------------------------------------------------------
# 2.  T-TESTS — Strategy Numeric Score: M vs F; APP vs WT
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 2: T-Tests on Strategy Score (StrategyNum)\n")
cat("============================================================\n\n")

# --- 2a. M vs F -------------------------------------------------------------
cat("--- 2a. Sex: Male vs Female (Days 2 & 4) ---\n")
tt_sex <- t.test(StrategyNum ~ Sex, data = d24, var.equal = FALSE)
print(tt_sex)

# Effect size (Cohen's d)
sex_d <- d24 %>%
  cohens_d(StrategyNum ~ Sex, var.equal = FALSE)
cat("Cohen's d:", round(sex_d$effsize, 3), "\n\n")

# --- 2b. APP/+ vs WT --------------------------------------------------------
cat("--- 2b. Genotype: APP/+ vs WT (Days 2 & 4) ---\n")
tt_app <- t.test(StrategyNum ~ Genotype, data = d24, var.equal = FALSE)
print(tt_app)

app_d <- d24 %>%
  cohens_d(StrategyNum ~ Genotype, var.equal = FALSE)
cat("Cohen's d:", round(app_d$effsize, 3), "\n\n")

# --- 2c. M vs F separately within each Day ----------------------------------
cat("--- 2c. Sex t-test: Day 2 vs Day 4 separately ---\n")
for (dy in c(2, 4)) {
  sub <- d24 %>% filter(Day == dy)
  tt  <- t.test(StrategyNum ~ Sex, data = sub, var.equal = FALSE)
  cat(sprintf("  Day %d: t(%.1f) = %.3f, p = %.4f\n",
              dy, tt$parameter, tt$statistic, tt$p.value))
}

# --- 2d. APP vs WT per Day --------------------------------------------------
cat("\n--- 2d. Genotype t-test: Day 2 vs Day 4 separately ---\n")
for (dy in c(2, 4)) {
  sub <- d24 %>% filter(Day == dy)
  tt  <- t.test(StrategyNum ~ Genotype, data = sub, var.equal = FALSE)
  cat(sprintf("  Day %d: t(%.1f) = %.3f, p = %.4f\n",
              dy, tt$parameter, tt$statistic, tt$p.value))
}

# ---------------------------------------------------------------------------
# 3.  ONE-WAY ANOVA — Strategy across Strategy groups (sanity / distribution)
#     Main question: do genotype / sex groups differ in StrategyNum?
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 3: One-Way ANOVAs — StrategyNum across Groups\n")
cat("============================================================\n\n")

# Helper: compute eta-squared from an aov summary (no extra packages needed)
eta2_from_aov <- function(aov_obj) {
  ss  <- summary(aov_obj)[[1]][["Sum Sq"]]
  # last element is residuals; everything before is effects
  eta2 <- ss[-length(ss)] / sum(ss)
  names(eta2) <- rownames(summary(aov_obj)[[1]])[-length(ss)]
  round(eta2, 4)
}

# --- 3a. Across Genotype ----------------------------------------------------
cat("--- 3a. One-way ANOVA: Genotype ---\n")
aov_geno <- aov(StrategyNum ~ Genotype, data = d24)
print(summary(aov_geno))
cat("Eta-squared:", eta2_from_aov(aov_geno), "\n\n")

# --- 3b. Across Sex ---------------------------------------------------------
cat("--- 3b. One-way ANOVA: Sex ---\n")
aov_sex <- aov(StrategyNum ~ Sex, data = d24)
print(summary(aov_sex))
cat("Eta-squared:", eta2_from_aov(aov_sex), "\n\n")

# --- 3c. Two-way ANOVA: Sex × Genotype (D2 & D4 combined) ------------------
cat("--- 3c. Two-way ANOVA: Sex × Genotype ---\n")
aov2 <- aov(StrategyNum ~ Sex * Genotype, data = d24)
print(summary(aov2))
cat("Eta-squared per term:\n")
print(eta2_from_aov(aov2))

# Levene's test for homogeneity of variance
cat("\nLevene's test (Sex × Genotype):\n")
print(leveneTest(StrategyNum ~ Sex * Genotype, data = d24))

# ---------------------------------------------------------------------------
# 4.  STRATEGY PROBABILITY — by Genotype, Sex, StrategyName
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 4: Strategy Probability (proportion of trials per group)\n")
cat("============================================================\n\n")

# Overall strategy probability
strat_prob_overall <- d24 %>%
  count(StrategyName) %>%
  mutate(Prob = round(n / sum(n), 4)) %>%
  arrange(desc(Prob))
cat("--- Overall ---\n")
print(strat_prob_overall)

# By Genotype
strat_prob_geno <- d24 %>%
  count(Genotype, StrategyName) %>%
  group_by(Genotype) %>%
  mutate(Prob = round(n / sum(n), 4)) %>%
  ungroup() %>%
  arrange(Genotype, desc(Prob))
cat("\n--- By Genotype ---\n")
print(strat_prob_geno)

# By Sex
strat_prob_sex <- d24 %>%
  count(Sex, StrategyName) %>%
  group_by(Sex) %>%
  mutate(Prob = round(n / sum(n), 4)) %>%
  ungroup() %>%
  arrange(Sex, desc(Prob))
cat("\n--- By Sex ---\n")
print(strat_prob_sex)

# By Sex × Genotype
strat_prob_full <- d24 %>%
  count(Sex, Genotype, StrategyName) %>%
  group_by(Sex, Genotype) %>%
  mutate(Prob = round(n / sum(n), 4)) %>%
  ungroup() %>%
  arrange(Sex, Genotype, desc(Prob))
cat("\n--- By Sex × Genotype ---\n")
print(strat_prob_full)

# ---------------------------------------------------------------------------
# 5.  IS THE SEX × GENOTYPE INTERACTION SIGNIFICANT?
#     Mixed ANOVA: StrategyNum ~ Sex * Genotype * Day + (1|SubjectID)
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 5: Sex × Genotype Interaction (mixed model)\n")
cat("============================================================\n\n")

# Aggregate to subject mean per day (removes trial-level noise)
subj_day <- d24 %>%
  group_by(SubjectID, Sex, Genotype, Day) %>%
  summarise(MeanStrat = mean(StrategyNum, na.rm = TRUE), .groups = "drop")

# Linear mixed model: fixed = Sex * Genotype * Day, random intercept = Subject
m_full <- lmer(MeanStrat ~ Sex * Genotype * Day + (1 | SubjectID),
               data = subj_day, REML = FALSE)
cat("--- Full mixed model: Sex * Genotype * Day ---\n")
print(anova(m_full))

# Test just Sex × Genotype interaction (collapsed across Day)
m_2way <- lmer(MeanStrat ~ Sex * Genotype + (1 | SubjectID),
               data = subj_day, REML = FALSE)
cat("\n--- Two-way interaction: Sex × Genotype ---\n")
print(anova(m_2way))

# Post-hoc: pairwise contrasts for Sex × Genotype
cat("\n--- Post-hoc: Sex × Genotype pairwise ---\n")
emm_sg <- emmeans(m_2way, ~ Sex * Genotype)
cat("Cell means:\n")
print(emm_sg)
cat("\nPairwise contrasts (Tukey):\n")
print(contrast(emm_sg, method = "pairwise", adjust = "tukey"))

# ---------------------------------------------------------------------------
# 6.  REPEATED MEASURES ANOVA
#     Strategy Use × Sex × Genotype  (subject = unit, Day as within factor)
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 6: Repeated Measures ANOVA — StrategyNum ~ Sex × Genotype × Day\n")
cat("============================================================\n\n")

# Using lmerTest for Type III ANOVA with Satterthwaite df
m_rm <- lmer(MeanStrat ~ Sex * Genotype * Day + (1 | SubjectID),
             data = subj_day, REML = TRUE)
cat("--- Type III ANOVA (Satterthwaite) ---\n")
print(anova(m_rm, type = 3))

# Marginal means
cat("\n--- Marginal means: Sex ---\n")
print(emmeans(m_rm, ~ Sex))

cat("\n--- Marginal means: Genotype ---\n")
print(emmeans(m_rm, ~ Genotype))

cat("\n--- Marginal means: Day ---\n")
print(emmeans(m_rm, ~ Day))

cat("\n--- Marginal means: Sex × Genotype ---\n")
print(emmeans(m_rm, ~ Sex * Genotype))

# ---------------------------------------------------------------------------
# 7.  PER-STRATEGY REPEATED MEASURES ANOVA
#     For each strategy TYPE: is probability ~ Sex × Genotype × Day significant?
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 7: Per-Strategy RM ANOVA (binary use of each strategy)\n")
cat("============================================================\n\n")

strategy_names <- levels(d24$StrategyName)

# Build subject × strategy × day binary use table
strat_use <- d24 %>%
  mutate(Dummy = 1) %>%
  group_by(SubjectID, Sex, Genotype, Day) %>%
  summarise(across(all_of(strategy_names),
                   ~ mean(StrategyName == cur_column()),
                   .names = "{.col}"),
            .groups = "drop")

rm_results <- list()

for (s in strategy_names) {
  formula <- as.formula(paste0("`", s, "` ~ Sex * Genotype * Day + (1 | SubjectID)"))
  tryCatch({
    m <- lmer(formula, data = strat_use, REML = TRUE)
    a <- anova(m, type = 3)
    rm_results[[s]] <- a
    cat(sprintf("\n=== Strategy: %s ===\n", s))
    print(a)
  }, error = function(e) {
    cat(sprintf("\n=== Strategy: %s — model failed: %s ===\n", s, e$message))
  })
}

# ---------------------------------------------------------------------------
# 8.  SEX DIFFERENCES — detailed breakdown
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 8: Sex Differences\n")
cat("============================================================\n\n")

cat("--- Mean StrategyNum by Sex × Day ---\n")
d24 %>%
  group_by(Sex, Day) %>%
  summarise(
    n    = n(),
    Mean = round(mean(StrategyNum, na.rm=TRUE), 3),
    SD   = round(sd(StrategyNum,   na.rm=TRUE), 3),
    SE   = round(SD / sqrt(n), 3),
    .groups = "drop"
  ) %>% print()

cat("\n--- T-test: StrategyNum M vs F per Day ---\n")
for (dy in c(2, 4)) {
  sub <- d24 %>% filter(Day == dy)
  tt  <- t.test(StrategyNum ~ Sex, data = sub, var.equal = FALSE)
  cd  <- cohens_d(sub, StrategyNum ~ Sex, var.equal = FALSE)
  cat(sprintf("Day %d: t(%.1f)=%.3f, p=%.4f, Cohen's d=%.3f\n",
              dy, tt$parameter, tt$statistic, tt$p.value, cd$effsize))
}

cat("\n--- Chi-square: Strategy choice frequency by Sex ---\n")
sex_tab <- table(d24$StrategyName, d24$Sex)
print(sex_tab)
chi_sex <- chisq.test(sex_tab)
print(chi_sex)

# ---------------------------------------------------------------------------
# 9.  GENOTYPE DIFFERENCES — APP/+ vs WT
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 9: Genotype Differences (APP/+ vs WT)\n")
cat("============================================================\n\n")

cat("--- Mean StrategyNum by Genotype × Day ---\n")
d24 %>%
  group_by(Genotype, Day) %>%
  summarise(
    n    = n(),
    Mean = round(mean(StrategyNum, na.rm=TRUE), 3),
    SD   = round(sd(StrategyNum,   na.rm=TRUE), 3),
    SE   = round(SD / sqrt(n), 3),
    .groups = "drop"
  ) %>% print()

cat("\n--- T-test: StrategyNum APP/+ vs WT per Day ---\n")
for (dy in c(2, 4)) {
  sub <- d24 %>% filter(Day == dy)
  tt  <- t.test(StrategyNum ~ Genotype, data = sub, var.equal = FALSE)
  cd  <- cohens_d(sub, StrategyNum ~ Genotype, var.equal = FALSE)
  cat(sprintf("Day %d: t(%.1f)=%.3f, p=%.4f, Cohen's d=%.3f\n",
              dy, tt$parameter, tt$statistic, tt$p.value, cd$effsize))
}

cat("\n--- Chi-square: Strategy choice frequency by Genotype ---\n")
geno_tab <- table(d24$StrategyName, d24$Genotype)
print(geno_tab)
chi_geno <- chisq.test(geno_tab)
print(chi_geno)

# ---------------------------------------------------------------------------
# 10.  OVERALL EFFECT OF ALZHEIMER'S ON STRATEGY USE
#      (Summarising APP effect across days, sex, and strategy type)
# ---------------------------------------------------------------------------

cat("\n============================================================\n")
cat("SECTION 10: Overall Effect of APP (Alzheimer's Model) on Strategy Use\n")
cat("============================================================\n\n")

# All 4 days, not just D2/D4
app_overall <- raw %>%
  group_by(Day, Genotype) %>%
  summarise(
    n         = n(),
    MeanStrat = round(mean(StrategyNum, na.rm=TRUE), 3),
    SD        = round(sd(StrategyNum,   na.rm=TRUE), 3),
    SE        = round(SD / sqrt(n), 3),
    .groups = "drop"
  )
cat("--- StrategyNum by Genotype across All Days ---\n")
print(app_overall)

# Full RM model across all 4 days
subj_allday <- raw %>%
  group_by(SubjectID, Sex, Genotype, Day) %>%
  summarise(MeanStrat = mean(StrategyNum, na.rm=TRUE), .groups="drop")

m_all <- lmer(MeanStrat ~ Genotype * Day + Sex + (1 | SubjectID),
              data = subj_allday, REML = TRUE)
cat("\n--- RM model (all days): Genotype × Day + Sex ---\n")
print(anova(m_all, type = 3))

cat("\n--- Post-hoc: Genotype × Day pairwise ---\n")
emm_gd <- emmeans(m_all, ~ Genotype * Day)
print(contrast(emm_gd, method = "pairwise", adjust = "tukey"))

# Proportion of "spatial" strategies (strategies 6–8: directed, corrected, direct)
cat("\n--- Proportion of spatial strategies (≥ directed search) by Genotype & Day ---\n")
raw %>%
  mutate(Spatial = StrategyNum >= 6) %>%
  group_by(Genotype, Day) %>%
  summarise(PctSpatial = round(mean(Spatial) * 100, 2), n = n(), .groups="drop") %>%
  print()

cat("\n=== Analysis complete. ===\n")