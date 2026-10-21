#code for thermal performance chapter :) 

library(tidyverse)
library(lubridate)
library(patchwork)
library(lme4)
library(car)
library(emmeans)
library(statmod) 
library(tweedie)  

#set contrasts for ANOVA TIII
options(contrasts = c("contr.sum", "contr.poly"))
set.seed(2025)

#create folder for the figures
dir.create("outputs", showWarnings = FALSE)

#95% confidence interval function
ci95 <- function(x) {             
  x <- na.omit(x); h <- qt(0.975, length(x) - 1) * sd(x) / sqrt(length(x))
  data.frame(y = mean(x), ymin = mean(x) - h, ymax = mean(x) + h)
}

#colour palette!!
pal       <- c(sericata = "#C87F17", vicina = "blue")

#species labels
sp_labels <- c(sericata = expression(italic("L. sericata")),
               vicina   = expression(italic("C. vicina")))

#same theme for all figs
theme_ch <- theme_classic(base_size = 16) +
  theme(axis.text = element_text(size = 16, colour = "black"),
        axis.title = element_text(size = 16, colour = "black"),
        legend.title = element_blank(), legend.text = element_text(size = 16),
        plot.tag = element_text(face = "bold", size = 22))


#Load the data
comp <- read_csv("comp_data.csv", na = c("N/A", "NA", ""), show_col_types = FALSE) %>%
  mutate(across(c(date_larvae_transferred_to_treatment, date_first_emerged, date_second_emerged),
                ~ suppressWarnings(dmy(.x))),
         density = case_when(control == "treatment"          ~ "Mixed",
                             str_detect(replicate, "^..10-") ~ "n = 10",
                             str_detect(replicate, "^..20-") ~ "n = 20"))

# one row per species per cage: larvae added, adults emerged, days to first adult
species_rows <- function(sp, code) {
  comp %>% transmute(
    replicate, temperature, control, density, species = sp,
    initial = case_when(control == "treatment"                       ~ 10,
                        str_detect(replicate, paste0("^", code, "10")) ~ 10,
                        str_detect(replicate, paste0("^", code, "20")) ~ 20,
                        TRUE ~ 0),
    success = .data[[paste0("final_number_observed_", sp)]],
    emerge  = case_when(species_first_emerged %in% c(sp, "both") ~ date_first_emerged,
                        species_second_emerged == sp            ~ date_second_emerged),
    days_to_first = as.numeric(emerge - date_larvae_transferred_to_treatment))
}

long <- bind_rows(species_rows("sericata", "LS"), species_rows("vicina", "CV")) %>%
  filter(initial > 0) %>%
  mutate(species = factor(species, levels = c("sericata", "vicina")),
         prop_success = success / initial,
         temp_f  = factor(temperature),
         density = factor(density, levels = c("Mixed", "n = 10", "n = 20")))

# body size: one row per cage x species x sex (cage means)
body <- comp %>%
  select(replicate, temperature, control, density,
         sericata_male = average_body_length_males_lucilia,
         sericata_female = average_body_length_females_lucilia,
         vicina_male = average_body_length_males_vicina,
         vicina_female = average_body_length_females_vicina) %>%
  pivot_longer(sericata_male:vicina_female, names_to = c("species", "sex"),
               names_sep = "_", values_to = "thorax_mm", values_drop_na = TRUE) %>%
  mutate(species = factor(species, levels = c("sericata", "vicina")),
         sex     = factor(sex, levels = c("male", "female")),
         temp_f  = factor(temperature),
         density = factor(density, levels = c("Mixed", "n = 10", "n = 20")))

fec <- read_csv("fecundity.csv", show_col_types = FALSE) %>%
  transmute(temperature = parse_number(treatment),
            species = factor(if_else(str_detect(treatment, "ls"), "sericata", "vicina"),
                             levels = c("sericata", "vicina")),
            offspring = lifetime_fecundity, pi_cage = performance_index) %>%
  mutate(temp_f = factor(temperature))

long_mix <- filter(long, control == "treatment")   # mixed-species temperature cages (i.e., no single species)
body_mix <- filter(body, control == "treatment")


#stats for the thermal performance part

# Emergence: binomial GLM - including random effect for cage led to a singularity warning (explained no variance in the data)
m_emerg <- glm(cbind(success, initial - success) ~ species * temp_f, family = binomial, data = long_mix)
Anova(m_emerg, type = "III")                               
emm_emerg <- emmeans(m_emerg, ~ species | temp_f)
summary(emm_emerg, type = "response"); pairs(emm_emerg)   

# Development time: linear model (30 °C dropped, no C. vicina emerged) - same as above RE the singularity warning
time_dat <- long_mix %>% filter(temperature != 30) %>% droplevels()
m_time <- lm(days_to_first ~ species * temp_f, data = time_dat)
Anova(m_time, type = "III")                               
emm_time <- emmeans(m_time, ~ species | temp_f)
emm_time; pairs(emm_time)                                 

# Body size: LMM with cage random intercept. I specify to use the F statistic because this will give X2 by default, however, F will corrects the degrees of freedom for uncertainty (just a bit more conservative in this case)
body_dat <- body_mix %>% filter(temperature != 30) %>% droplevels()
m_body <- lmer(thorax_mm ~ species * temp_f * sex + (1 | replicate), data = body_dat)
Anova(m_body, type = "III", test.statistic = "F")         
emm_body <- emmeans(m_body, ~ species | temp_f * sex)
emm_body; pairs(emm_body)                                

# Offspring: Tweedie GLM, power parameter chosen by profile likelihood
best_p <- tweedie.profile(offspring ~ species * temp_f, data = fec, p.vec = seq(1.1, 1.9, by = 0.05),
                          link.power = 0, do.plot = FALSE)$p.max
best_p
m_eggs <- glm(offspring ~ species * temp_f, data = fec,
              family = tweedie(var.power = best_p, link.power = 0))
Anova(m_eggs, type = "III")                                # Table 4.1


##PERFORMANCE INDEX (with bootstrap)
# PI = ln(F * S) / D from species x temperature means (0 when F = 0). - fecundity cages and emergence cages resampled separately, 10,000 times.
pi_fun <- function(F, S, D) { FS <- mean(F) * mean(S); if (FS > 0) log(FS) / mean(D) else 0 }
rs     <- function(x) x[sample.int(length(x), replace = TRUE)]

pi_summary <- expand_grid(species = c("sericata", "vicina"), temperature = sort(unique(fec$temperature))) %>%
  pmap_dfr(function(species, temperature) {
    F <- fec$offspring[fec$species == species & fec$temperature == temperature]
    e <- filter(long_mix, species == !!species, temperature == !!temperature)
    S <- e$prop_success; D <- na.omit(e$days_to_first)
    boot <- replicate(10000, pi_fun(rs(F), rs(S), rs(D)))
    tibble(species, temperature, PI = pi_fun(F, S, D),
           lo = if (length(F) > 1) quantile(boot, 0.025) else PI,   # n = 1 female: no CI
           hi = if (length(F) > 1) quantile(boot, 0.975) else PI)
  }) %>%
  mutate(species = factor(species, levels = c("sericata", "vicina")))
pi_summary


# mixed v species cages at 22 degrees :) (for the supplementary stuff)
# One-way ANOVAs (emergence, development, offspring, PI); two-way density x sex
# ANOVA (Type I) for body size; Tukey HSD where significant.
long_22 <- filter(long, temperature == 22)
body_22 <- filter(body, temperature == 22)
fec22 <- read_csv("fecun_sum_22.csv", show_col_types = FALSE) %>%
  transmute(species = factor(if_else(str_detect(tolower(treatment), "ls"), "sericata", "vicina"),
                             levels = c("sericata", "vicina")),
            density = factor(case_when(str_detect(treatment, "20") ~ "n = 20",
                                       str_detect(treatment, "10") ~ "n = 10",
                                       TRUE                        ~ "Mixed"),
                             levels = c("Mixed", "n = 10", "n = 20")),
            offspring = lifetime_fecundity, pi_cage = performance_index)

for (sp in c("vicina", "sericata")) {
  cat("\n==========", sp, "==========\n")
  d <- filter(long_22, species == sp); b <- filter(body_22, species == sp); f <- filter(fec22, species == sp)
  print(summary(aov(prop_success  ~ density,       data = d)))
  print(summary(aov(days_to_first ~ density,       data = d)))
  print(summary(aov(thorax_mm     ~ density * sex, data = b)))
  print(summary(aov(offspring     ~ density,       data = f)))
  print(summary(aov(pi_cage       ~ density,       data = f)))
}
TukeyHSD(aov(prop_success  ~ density, data = filter(long_22, species == "sericata")))
TukeyHSD(aov(days_to_first ~ density, data = filter(long_22, species == "sericata")))
TukeyHSD(aov(thorax_mm ~ density * sex, data = filter(body_22, species == "sericata")), "density")


#ACUTE STUFF
kd <- read.csv("thermal_tolerance_data.csv") %>%
  mutate(sex = factor(sex, levels = c("M", "F")), species = factor(species),
         collection_date = factor(collection_date))

kd_models <- map(c(cold = "cold", heat = "heat"), function(x)
  lmer(log(coma_time) ~ species * sex + thorax_length + (1 | collection_date),
       data = filter(kd, cold.heat == x)))

map(kd_models, Anova, type = "III")                                           
kd_emm <- map(kd_models, ~ emmeans(.x, ~ species * sex))
map(kd_emm, summary, type = "response")                                        
map(kd_emm, pairs, type = "response")                                          
with(as.data.frame(VarCorr(kd_models$cold)), vcov[1] / sum(vcov))             


###FIGUREEEEEES

line_panel <- function(d, y, lo, hi, ylab, breaks = c(15, 19, 22, 25, 30)) {
  ggplot(d, aes(temperature, {{ y }}, colour = species)) +
    geom_ribbon(aes(ymin = {{ lo }}, ymax = {{ hi }}, fill = species), alpha = 0.2, colour = NA) +
    geom_line(linewidth = 1.2) + geom_point(size = 6) +
    scale_colour_manual(values = pal, labels = sp_labels) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(breaks = breaks) +
    labs(x = "Temperature (°C)", y = ylab) + theme_ch
}

##emergence (A), development time (B); mean ± SE
emerg_sum <- long_mix %>% group_by(species, temperature) %>%
  summarise(p = sum(success) / sum(initial), se = sqrt(p * (1 - p) / sum(initial)), .groups = "drop")
time_sum <- long_mix %>% filter(!is.na(days_to_first)) %>% group_by(species, temperature) %>%
  summarise(m = mean(days_to_first), se = sd(days_to_first) / sqrt(n()), .groups = "drop")

fig4.3 <- (line_panel(emerg_sum, p, p - se, p + se, "Emergence proportion") +
             scale_y_continuous(limits = c(0, 1), labels = scales::percent) |
           line_panel(time_sum, m, m - se, m + se, "Days to first emergence")) +
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
fig4.3
ggsave("outputs/Fig4.3_emergence_development.png", fig4.3, width = 14, height = 6.5, dpi = 600)

## Fig 4.4: body size (A), offspring (B), performance index (C)
body_sum <- body_mix %>% group_by(species, sex, temperature) %>%
  summarise(m = mean(thorax_mm), se = replace_na(sd(thorax_mm) / sqrt(n()), 0), .groups = "drop") %>%
  mutate(grp = factor(paste(species, sex), levels = c("sericata female", "sericata male",
                                                      "vicina female", "vicina male")))
grp_labels <- c(expression(italic("L. sericata") * " (F)"), expression(italic("L. sericata") * " (M)"),
                expression(italic("C. vicina") * " (F)"),   expression(italic("C. vicina") * " (M)"))

fig4.4A <- ggplot(body_sum, aes(temperature, m, colour = grp, linetype = grp)) +
  geom_ribbon(aes(ymin = m - se, ymax = m + se, fill = grp), alpha = 0.3, colour = NA) +
  geom_line(linewidth = 1.2) + geom_point(size = 6) +
  scale_colour_manual(values = unname(pal[c(1, 1, 2, 2)]), labels = grp_labels) +
  scale_linetype_manual(values = c("dotted", "solid", "dotted", "solid"), labels = grp_labels) +
  scale_fill_manual(values = unname(pal[c(1, 1, 2, 2)]), guide = "none") +
  scale_x_continuous(breaks = c(15, 19, 22, 25, 30)) +
  labs(x = "Temperature (°C)", y = "Thorax length (mm)") + theme_ch +
  theme(legend.position = "bottom", legend.key.width = unit(4, "lines")) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
fig4.4A

#eggs
eggs_sum <- fec %>% group_by(species, temperature) %>%
  summarise(m = mean(offspring), se = replace_na(sd(offspring) / sqrt(n()), 0), .groups = "drop")

fig4.4 <- (fig4.4A |
           (line_panel(eggs_sum, m, m - se, m + se, "Mean offspring laid", c(15, 19, 22, 25)) +
            line_panel(pi_summary, PI, lo, hi, "Performance index", c(15, 19, 22, 25)) +
            plot_layout(guides = "collect"))) +
  plot_layout(widths = c(1, 2)) + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")
fig4.4
ggsave("outputs/Fig4.4_bodysize_offspring_performance.png", fig4.4, width = 18, height = 7.5, dpi = 600)

## Fig 4.4 (acute): knockdown EMMs ± 95% CI
kd_df <- map_dfr(kd_emm, ~ as.data.frame(summary(.x, type = "response")), .id = "assay") %>%
  rename(any_of(c(lower.CL = "asymp.LCL", upper.CL = "asymp.UCL"))) %>%
  mutate(panel = if_else(assay == "cold", "A", "B"))

fig_kd <- ggplot(kd_df, aes(species, response / 60, colour = species, shape = sex)) +
  geom_errorbar(aes(ymin = lower.CL / 60, ymax = upper.CL / 60), width = 0.2, linewidth = 1,
                position = position_dodge(0.4)) +
  geom_point(size = 5, stroke = 1.2, position = position_dodge(0.4)) +
  facet_wrap(~ panel, scales = "free") +
  scale_colour_manual(values = c(Calliphora_vicina = "blue", Lucilia_sericata = "#C87F17"), guide = "none") +
  scale_shape_manual(values = c(M = 16, F = 17), labels = c(M = "Male", F = "Female")) +
  scale_x_discrete(labels = c(Calliphora_vicina = "C. vicina", Lucilia_sericata = "L. sericata")) +
  labs(x = NULL, y = "Time to knockdown (min)") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(face = "italic", colour = "black"),
        axis.text.y = element_text(colour = "black"), legend.title = element_blank(),
        strip.text = element_text(face = "bold", size = 22, hjust = 0),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
fig_kd
ggsave("outputs/Fig4.4_knockdown.png", fig_kd, width = 10, height = 6, dpi = 600)


##supp figss
## thermal curves
cage_data <- bind_rows(
  long_mix %>% transmute(species, temperature, trait = "emergence", value = prop_success),
  long_mix %>% transmute(species, temperature, trait = "days",      value = days_to_first),
  body_mix %>% transmute(species, temperature, trait = paste0(sex, "_thorax"), value = thorax_mm),
  fec      %>% transmute(species, temperature, trait = "offspring", value = offspring),
  fec      %>% transmute(species, temperature, trait = "pi",        value = pi_cage)
) %>% filter(!is.na(value))

quad_curve <- function(d) {
  grid <- tibble(temperature = seq(min(d$temperature), max(d$temperature), length.out = 200))
  grid$fit <- predict(lm(value ~ I(temperature - 20) + I((temperature - 20)^2), data = d), grid)
  grid
}
curves <- bind_rows(list("All temperatures" = cage_data,
                         "Excluding 22 °C"  = filter(cage_data, temperature != 22)), .id = "fit_type") %>%
  group_by(fit_type, trait, species) %>% filter(n_distinct(temperature) >= 3) %>%
  group_modify(~ quad_curve(.x)) %>% ungroup() %>%
  mutate(fit = case_when(str_detect(trait, "thorax") ~ fit,
                         trait == "emergence"         ~ pmin(pmax(fit, 0), 1),
                         TRUE                         ~ pmax(fit, 0)))
means <- cage_data %>% group_by(trait, species, temperature) %>%
  summarise(value = mean(value), .groups = "drop") %>%
  mutate(rearing = if_else(temperature == 22, "Climate room (22 °C)", "Environmental chamber"))

tpc_panel <- function(key, ylab) {
  ggplot(mapping = aes(temperature, colour = species)) +
    geom_jitter(data = filter(cage_data, trait == key), aes(y = value),
                width = 0.25, height = 0, alpha = 0.25, size = 3, show.legend = FALSE) +
    geom_line(data = filter(curves, trait == key), aes(y = fit, linetype = fit_type), linewidth = 1.2) +
    geom_point(data = filter(means, trait == key), aes(y = value, shape = rearing),
               fill = "white", size = 6, stroke = 1.3) +
    scale_colour_manual(values = pal, labels = sp_labels) +
    scale_shape_manual(values = c("Environmental chamber" = 16, "Climate room (22 °C)" = 21)) +
    scale_linetype_manual(values = c("All temperatures" = "solid", "Excluding 22 °C" = "dashed")) +
    scale_x_continuous(breaks = c(15, 19, 22, 25, 30)) +
    guides(colour = guide_legend(order = 1), shape = guide_legend(order = 2),
           linetype = guide_legend(order = 3, override.aes = list(colour = "black"))) +
    labs(x = "Temperature (°C)", y = ylab) + theme_ch
}

figA3.1 <- tpc_panel("emergence", "Emergence proportion") + scale_y_continuous(labels = scales::percent) +
  tpc_panel("days", "Days to first emergence") +
  tpc_panel("female_thorax", "Female thorax length (mm)") +
  tpc_panel("male_thorax", "Male thorax length (mm)") +
  tpc_panel("offspring", "Mean offspring laid") +
  tpc_panel("pi", "Performance index") +
  plot_layout(ncol = 3, guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.text = element_text(size = 14))
figA3.1

ggsave("outputs/FigA3.1_thermal_performance_curves.png", figA3.1, width = 20, height = 13, dpi = 600)

## Figs A3.2-A3.3
cage_panel <- function(d, var, ylab, ...) {
  d %>% mutate(species = factor(species, levels = c("vicina", "sericata"))) %>%
    ggplot(aes(density, {{ var }}, colour = species, ...)) +
    stat_summary(fun.data = ci95, geom = "errorbar", width = 0.15, linewidth = 0.8,
                 position = position_dodge(0.5)) +
    stat_summary(fun = mean, geom = "point", size = 5, position = position_dodge(0.5)) +
    facet_wrap(~ species, labeller = as_labeller(c(vicina = "C. vicina", sericata = "L. sericata"))) +
    scale_colour_manual(values = pal, guide = "none") +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 16) +
    theme(panel.border = element_rect(fill = NA, linewidth = 1.2), axis.text = element_text(colour = "black"),
          strip.background = element_blank(), strip.text = element_text(face = "bold.italic"),
          legend.position = "bottom", legend.title = element_blank(),
          plot.tag = element_text(face = "bold", size = 22))
}

figA3.2 <- cage_panel(long_22, prop_success, "Emergence proportion (22 °C)") +
  scale_y_continuous(labels = scales::percent) + coord_cartesian(ylim = c(0, 1)) +
  cage_panel(long_22, days_to_first, "Days to first emergence (22 °C)") +
  plot_layout(ncol = 1) + plot_annotation(tag_levels = "A")
figA3.2
ggsave("outputs/FigA3.2_22C_emergence_development.png", figA3.2, width = 8, height = 10, dpi = 600)

figA3.3 <- (cage_panel(mutate(body_22, sex = str_to_title(sex)), thorax_mm, "Thorax length (mm; 22 °C)",
                       shape = sex) + scale_shape_manual(values = c(Male = 16, Female = 17)) +
            cage_panel(fec22, offspring, "Mean offspring laid (22 °C)")) /
  (plot_spacer() + cage_panel(fec22, pi_cage, "Performance index (22 °C)") + plot_spacer() +
     plot_layout(widths = c(1, 2, 1))) +
  plot_annotation(tag_levels = "A")
figA3.3
ggsave("outputs/FigA3.3_22C_body_offspring_performance.png", figA3.3, width = 14, height = 11, dpi = 600)
