library(janitor)

#load data and create output directory
infile <- "comp_data.csv"
outdir <- "outputs"

#LOAD SOME HELPER FUNCTIONS 
#species name normalisation:
norm_sp <- function(x) tolower(trimws(as.character(x)))

#emergence date normalisation
#this function translates cage-level emergence info into species level dates and adds columns for the day that each species emerged
derive_emergence_dates <- function(dat) {
  dat %>%
    mutate(
      s1 = norm_sp(species_first_emerged),
      s2 = norm_sp(species_second_emerged),
      
      first_emerge_sericata = case_when(
        s1 %in% c("sericata","both") ~ date_first_emerged,
        s2 == "sericata"             ~ date_second_emerged,
        TRUE ~ as_date(NA)
      ),
      
      first_emerge_vicina = case_when(
        s1 %in% c("vicina","both") ~ date_first_emerged,
        s2 == "vicina"             ~ date_second_emerged,
        TRUE ~ as_date(NA)
      )
    )
}

#initial larval counts 
add_initial_counts <- function(dat, treat_init = 10L) {
  dat %>%
    mutate(
      is_control = tolower(control) == "control",
      is_treat   = tolower(control) == "treatment",
      
      ls_n = if_else(is_control & str_detect(replicate, "^LS"),
                     as.integer(str_extract(replicate, "(?<=^LS)\\d+")),
                     NA_integer_),
      
      cv_n = if_else(is_control & str_detect(replicate, "^CV"),
                     as.integer(str_extract(replicate, "(?<=^CV)\\d+")),
                     NA_integer_),
      
      initial_sericata = case_when(
        is_control & !is.na(ls_n) ~ ls_n,
        is_treat                 ~ treat_init,
        TRUE ~ 0L
      ),
      
      initial_vicina = case_when(
        is_control & !is.na(cv_n) ~ cv_n,
        is_treat                 ~ treat_init,
        TRUE ~ 0L
      )
    )
}

#long format conversion (per species per cage)
make_long <- function(dat) {
  bind_rows(
    dat %>%
      transmute(
        replicate, temperature, control = factor(control),
        species = "sericata",
        initial = initial_sericata,
        success = final_number_observed_sericata,
        days_to_first = days_to_first_sericata
      ),
    dat %>%
      transmute(
        replicate, temperature, control = factor(control),
        species = "vicina",
        initial = initial_vicina,
        success = final_number_observed_vicina,
        days_to_first = days_to_first_vicina
      )
  ) %>%
    mutate(
      species = factor(species, levels = c("sericata","vicina")),
      prop_success = if_else(initial > 0, success / initial, NA_real_)
    )
}

#wilson confidence interval - jitter instead? this also calculates mean probability (I don't end up plotting the CI)
wilson_ci <- function(x, n, z = 1.96) {
  p <- x / n
  denom  <- 1 + z^2 / n
  center <- (p + z^2 / (2*n)) / denom
  half   <- z * sqrt((p*(1 - p) + z^2/(4*n)) / n) / denom
  
  tibble(
    prob = p,
    ci_low  = pmax(0, center - half),
    ci_high = pmin(1, center + half)
  )
}

#LOAD, CLEAN, AND PREP THE DATA
dat <- read_csv(infile, na = c("N/A","NA","")) %>%
  clean_names()

date_cols <- c(
  "date_offered_meat",
  "date_larvae_transferred_to_treatment",
  "date_first_pupa",
  "date_first_emerged",
  "date_second_emerged"
)

#convert date columns (from date_cols) to a proper date format - dd-mm-yyyy and converts them from character --> date
#converts dates into development times (days)
dat <- dat %>%
  mutate(across(all_of(date_cols), ~ suppressWarnings(dmy(.x)))) %>%
  derive_emergence_dates() %>%
  mutate(
    days_to_pupa = as.numeric(date_first_pupa - date_larvae_transferred_to_treatment),
    days_to_first_sericata =
      as.numeric(first_emerge_sericata - date_larvae_transferred_to_treatment),
    days_to_first_vicina =
      as.numeric(first_emerge_vicina - date_larvae_transferred_to_treatment)
  ) %>%
  add_initial_counts()

#common theme for all of the figures
common_theme_fig2 <- theme_classic(base_size = 16) +
  theme(
    axis.line         = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks.length = unit(4, "pt"),
    axis.text         = element_text(size = 16, colour = "black"),
    axis.title        = element_text(size = 16, colour = "black"),
    legend.position   = NULL,
    legend.title = element_blank(),
    legend.text       = element_text(size = 16, face = "italic"),
    legend.key.height = unit(1.4, "lines"),
    legend.box.margin = margin(0, 0, 0, -50),
    plot.margin       = margin(10, 12, 10, 10)
  )

#########################################################################
##########################FIGURE 1#######################################
#########################################################################
#EMERGENCE SUCCESS
long <- make_long(dat)

agg_treat <- long %>%
  dplyr::filter(control == "treatment", initial > 0) %>%
  group_by(species, temperature) %>%
  summarise(
    successes = sum(success, na.rm = TRUE),
    trials    = sum(initial, na.rm = TRUE),
    n_reps    = n_distinct(replicate),
    .groups = "drop"
  ) %>%
  bind_cols(wilson_ci(.$successes, .$trials))

#calculate se
agg_treat <- agg_treat %>%
  mutate(
    se = sqrt(prob * (1 - prob) / trials)
  )

#Colours
pal_species <- c(sericata = "#C87F17", vicina = "blue")

# Plot treatments only with jitter
p_treat <- ggplot() +
  # raw cage-level data as jitter
  geom_jitter(
    data = long,
    aes(
      x = temperature,
      y = prop_success,
      colour = species
    ),
    width = 0.25,    # horizontal staggering
    height = 0,
    alpha = 0.2,
    size = 4
  ) +
  # species means
  geom_point(
    data = agg_treat,
    aes(
      x = temperature,
      y = prob,
      colour = species
    ),
    size = 6,
    shape = 16
  ) +
  # species mean lines
  geom_line(
    data = agg_treat,
    aes(
      x = temperature,
      y = prob,
      colour = species,
      group = species
    ),
    linewidth = 1.2
  ) +
  scale_color_manual(
    values = pal_species,
    labels = c(sericata = "L. sericata", vicina = "C. vicina"),
    name = NULL
  ) +
  scale_x_continuous(breaks = sort(unique(agg_treat$temperature))) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Emergence proportion"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        linewidth = 1.2,
        shape = 16,
        linetype = "solid"
      ),
      keywidth = unit(2, "lines")  
    )
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    legend.position   = "right",
    legend.text       = element_text(size = 16, face = "italic"),
    legend.title      = element_text(size = 13),
    legend.box.margin = margin(0, 0, 0, -50),
    legend.key.height = unit(1.5, "lines"),
    plot.margin       = margin(10, 14, 10, 10),
    plot.title        = element_text(hjust = 0.5, size = 14.5)
  )
p_treat


# Plot treatments only with no jitter
p_treat_nojitter <- ggplot() +
  geom_point(
    data = agg_treat,
    aes(
      x = temperature,
      y = prob,
      colour = species
    ),
    size = 6,
    shape = 16
  ) +
  geom_ribbon(
    data = agg_treat,
    aes(
      x = temperature,
      ymin = prob - se,
      ymax = prob + se,
      fill = species,
      group = species,
    ),
    alpha = 0.2,
    colour = NA,
    show.legend = FALSE
  ) +
  # species mean lines
  geom_line(
    data = agg_treat,
    aes(
      x = temperature,
      y = prob,
      colour = species,
      group = species
    ),
    linewidth = 1.2
  ) +
scale_color_manual(
  values = pal_species,
  labels = c(sericata = "L. sericata", vicina = "C. vicina"),
  name = NULL
) +

scale_fill_manual(
  values = pal_species,
  guide = "none"
) +
  scale_x_continuous(breaks = sort(unique(agg_treat$temperature))) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Emergence proportion"
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        linewidth = 1.2,
        shape = 16,
        linetype = "solid"
      ),
      keywidth = unit(2, "lines")  
    )
  ) + common_theme_fig2

p_treat_nojitter

ggsave(
  filename = file.path(outdir, "emergence_success_only_treatments.png"),
  plot = p_treat_nojitter, width = 7.5, height = 6, dpi = 600
)

#-----------------------------
# EMERGENCE TIMING (treatment only)
#-----------------------------

# summary stats (means only; CI calculated but not plotted)
timing_summary <- long %>%
  dplyr::filter(control == "treatment", initial > 0) %>%
  group_by(species, temperature) %>%
  summarise(
    n = n(),
    mean_days = mean(days_to_first, na.rm = TRUE),
    sd_days   = sd(days_to_first, na.rm = TRUE),
    se_days   = sd_days / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    ci_low  = mean_days - 1.96 * se_days,
    ci_high = mean_days + 1.96 * se_days
  )

# labels & colours
species_labels <- c(
  sericata = expression(italic("L. sericata")),
  vicina   = expression(italic("C. vicina"))
)

species_colors <- c(
  sericata = "#C87F17",
  vicina   = "blue"
)

# plot
emergence_time_nojitter <- ggplot() +
  
  # species means
  geom_point(
    data = timing_summary,
    aes(
      x = temperature,
      y = mean_days,
      colour = species
    ),
    shape = 19,
    size = 6
  ) +
  
  # SE ribbon (NEW)
  geom_ribbon(
    data = timing_summary,
    aes(
      x = temperature,
      ymin = mean_days - se_days,
      ymax = mean_days + se_days,
      fill = species,
      group = species
    ),
    alpha = 0.2,
    colour = NA,
    show.legend = FALSE
  ) +
  
  # species mean lines
  geom_line(
    data = timing_summary,
    aes(
      x = temperature,
      y = mean_days,
      group = species,
      colour = species
    ),
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = species_colors,
    labels = species_labels,
    breaks = c("sericata", "vicina"),
    name = NULL
  ) +
  
  # fill matches colour, but hidden from legend
  scale_fill_manual(
    values = species_colors,
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = sort(unique(timing_summary$temperature))
  ) +
  
  labs(
    x = "Temperature (°C)",
    y = "Days to first emergence"
  ) +
  
  common_theme_fig2

emergence_time  # warning expected due to missing 30 °C data

#plot no jitter
emergence_time_nojitter <- ggplot() +
  
  # species means
  geom_point(
    data = timing_summary,
    aes(
      x = temperature,
      y = mean_days,
      colour = species
    ),
    shape = 19,
    size = 6
  ) +
  # SE ribbon (NEW)
  geom_ribbon(
    data = timing_summary,
    aes(
      x = temperature,
      ymin = mean_days - se_days,
      ymax = mean_days + se_days,
      fill = species,
      group = species
    ),
    alpha = 0.2,
    colour = NA,
    show.legend = FALSE
  ) +
  
  # species mean lines
  geom_line(
    data = timing_summary,
    aes(
      x = temperature,
      y = mean_days,
      group = species,
      colour = species
    ),
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = species_colors,
    labels = species_labels,
    breaks = c("sericata", "vicina"),
    name = NULL
  ) +
  
  # fill matches colour, but hidden from legend
  scale_fill_manual(
    values = species_colors,
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = sort(unique(timing_summary$temperature))
  ) +
  
  labs(
    x = "Temperature (°C)",
    y = "Days to first emergence"
  ) +
  
  common_theme_fig2

emergence_time_nojitter

ggsave(
  filename = file.path(outdir, "emergence_timing.png"),
  plot = emergence_time,
  width = 7.5,
  height = 6,
  dpi = 600
)

library(cowplot)

################PLOT BOTH TOGETHER
pA <- p_treat_nojitter +
  theme(legend.position = "none",
        plot.margin = margin(10, 10, 10, 10))   # TLBR

pB <- emergence_time_nojitter +
  theme(legend.position = "none",
        plot.margin = margin(10, 10, 10, 10))

## 2) Build a horizontal, centered legend from one plot
legend_bottom <- cowplot::get_legend(
  p_treat_nojitter +
    guides(
      colour = guide_legend(
        nrow = 1, byrow = TRUE,
        keywidth = unit(2, "lines") # <- increase to spread items (try 16–30 pt)
      )
    ) +
    theme(
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.justification = "center",
      legend.box.just      = "center",
      legend.margin     = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.text       = element_text(
        size = 16, face = "italic",
        margin = margin(r = 15)     # <- extra padding after first label
      ),
      legend.key.height = unit(1.2, "lines")
      # optional: add a bit here too → legend.spacing.x = unit(6, "pt")
    )
)

## 3) Two panels with equal panel areas (align + axis)
panels_AB <- cowplot::plot_grid(
  pA, pB,
  ncol = 2,
  align = "hv",
  axis  = "tblr",
  rel_widths = c(1, 1),
  labels = c("A", "B"),
  label_size = 22, label_fontface = "bold",
  label_x = -0.01, label_y = 1.01
)

## 4) Stack panels over legend and center the legend row
fig_AB <- cowplot::plot_grid(
  panels_AB,
  legend_bottom,
  ncol = 1,
  rel_heights = c(1, 0.1)  
)

fig_AB

ggsave(
  file.path(outdir, "Figure_AB_bottom_legend.png"),
  plot = fig_AB, dpi = 600
)



####################################################################
##############FIGURE 2#############################################
###################################################################

#Body size

# Tidy to long: one row per cage × species × sex
body_long <- dat %>%
  transmute(
    replicate, temperature, control = factor(control), 
    lucilia_male   = as.numeric(average_body_length_males_lucilia),
    lucilia_female = as.numeric(average_body_length_females_lucilia),
    vicina_male    = as.numeric(average_body_length_males_vicina),
    vicina_female  = as.numeric(average_body_length_females_vicina)
  ) %>%
  pivot_longer(
    cols = c(lucilia_male, lucilia_female, vicina_male, vicina_female),
    names_to = "key", values_to = "thorax_mm"
  ) %>%
  mutate(
    species = case_when(
      str_detect(key, "^lucilia_") ~ "sericata",
      str_detect(key, "^vicina_")  ~ "vicina",
      TRUE ~ NA_character_
    ),
    sex = if_else(str_detect(key, "_male$"), "male", "female"),
    species = factor(species, levels = c("sericata","vicina")),
    sex     = factor(sex, levels = c("male","female"))
  ) %>%
  dplyr::select(-key)

# Treatments only
body_treat <- body_long %>%
  dplyr::filter(tolower(control) == "treatment", !is.na(thorax_mm))

##### SUMMARY (means only) #####
summary_body <- body_treat %>%
  group_by(species, sex, temperature) %>%
  summarise(
    n       = n(),
    mean_mm = mean(thorax_mm, na.rm = TRUE),
    sd_mm   = sd(thorax_mm,   na.rm = TRUE),
    se_mm   = sd_mm / sqrt(n),
    .groups = "drop"
  )

##### Species × sex combo factor #####
summary_body <- summary_body %>%
  mutate(
    grp = factor(
      paste(species, sex, sep = "_"),
      levels = c("sericata_female","sericata_male",
                 "vicina_female","vicina_male")
    )
  )

body_treat <- body_treat %>%
  mutate(grp = factor(paste(species, sex, sep = "_"),
                      levels = levels(summary_body$grp)))

##### Aesthetics #####
col_by_grp <- c(
  sericata_female = "#C87F17",
  sericata_male   = "#C87F17",
  vicina_female   = "blue",
  vicina_male     = "blue"
)

shape_by_grp <- c(
  sericata_female = 16,
  sericata_male   = 16,
  vicina_female   = 16,
  vicina_male     = 16
)

labels_by_grp <- c(
  sericata_female = expression(italic("L. sericata") * " (F)   "),
  sericata_male   = expression(italic("L. sericata") * " (M)"),
  vicina_female   = expression(italic("C. vicina") * " (F)"),
  vicina_male     = expression(italic("C. vicina") * " (M)")
)

##### PLOT #####
p_body <- ggplot(
  summary_body,
  aes(
    x = temperature,
    y = mean_mm,
    colour   = grp,
    shape    = grp,
    linetype = sex,
    group    = interaction(species, sex)
  )
) +
  # raw cage-level data FIRST (goes behind)
  geom_jitter(
    data = body_treat,
    aes(x = temperature, y = thorax_mm, colour = grp, shape = grp),
    width = 0.12,
    height = 0,
    alpha = 0.3,
    size  = 4,
    inherit.aes = FALSE
  ) +
  # species × sex means
  geom_line(linewidth = 1.2) +
  geom_point(size = 6) +
  scale_color_manual(values = col_by_grp, labels = labels_by_grp, name = "Species × Sex") +
  scale_shape_manual(values = shape_by_grp, labels = labels_by_grp, name = "Species × Sex") +
  scale_linetype_manual(values = c(male = "solid", female = "dotted"), guide = "none") +
  guides(
    colour = guide_legend(
      override.aes = list(shape = NA, linewidth = 1.2, linetype = c("solid","dotted","solid","dotted"))
    ),
    shape = "none"
  ) +
  scale_x_continuous(breaks = sort(unique(summary_body$temperature))) +
  labs(x = "Temperature (°C)", y = "Thorax length (mm)") +
  theme_classic(base_size = 16) +
  theme(
    axis.line         = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks        = element_line(linewidth = 0.6, colour = "black"),
    axis.ticks.length = unit(4, "pt"),
    legend.title      = element_blank(),
    legend.position   = "right",
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.text       = element_text(size = 16),
    legend.spacing.y  = unit(6, "pt"),
    legend.key.height = unit(1.4, "lines"),
    plot.margin       = margin(10, 12, 10, 10)
  )

p_body
###PLOT NO JITTER

##### PLOT #####
p_body_nojitter <- ggplot(
  summary_body,
  aes(
    x = temperature,
    y = mean_mm,
    colour   = grp,
    shape    = grp,
    linetype = sex,
    group    = interaction(species, sex)
  )
) +
  
  # -----------------------------
# SE RIBBON (SAFE: no legend impact)
# -----------------------------
geom_ribbon(
  data = summary_body,
  aes(
    x = temperature,
    ymin = mean_mm - se_mm,
    ymax = mean_mm + se_mm,
    fill = grp,
    group = interaction(species, sex)
  ),
  alpha = 0.3,
  colour = NA,
  show.legend = FALSE
) +
  
  # species × sex means
  geom_line(linewidth = 1.2) +
  geom_point(size = 6) +
  
  scale_color_manual(
    values = col_by_grp,
    labels = labels_by_grp,
    name = "Species × Sex"
  ) +
  scale_fill_manual(
    values = col_by_grp,
    guide = "none"
  ) +
  
  scale_shape_manual(
    values = shape_by_grp,
    labels = labels_by_grp,
    name = "Species × Sex"
  ) +
  
  scale_linetype_manual(
    values = c(male = "solid", female = "dotted"),
    guide = "none"
  ) +
  
  guides(
    colour = guide_legend(
      override.aes = list(
        shape = NA,
        linewidth = 1.2,
        linetype = c("dotted","solid","dotted","solid")
      )
    ),
    shape = "none"
  ) +
  
  scale_x_continuous(breaks = sort(unique(summary_body$temperature))) +
  
  labs(
    x = "Temperature (°C)",
    y = "Thorax length (mm)"
  ) +
  
  common_theme_fig2
p_body_nojitter

#eggs
# read data
eggs_summary <- read.csv("fecun_sum_full.csv")

eggs_summary <- eggs_summary %>%
  mutate(
    species = ifelse(grepl("ls", treatment), "L. sericata", "C. vicina"),
    temperature = as.numeric(substr(treatment, 1, 2))
  )

# species order
eggs_summary$species <- factor(
  eggs_summary$species,
  levels = c("L. sericata", "C. vicina")
)
# Summarise mean eggs per species per temperature
mean_eggs <- eggs_summary %>%
  group_by(temperature, species) %>%
  summarise(
    mean_eggs = mean(lifetime_fecundity, na.rm = TRUE),
    sd_eggs   = sd(lifetime_fecundity, na.rm = TRUE),
    n         = n(),
    se_eggs   = sd_eggs / sqrt(n),
    .groups = "drop"
  )

#so that the ribbon still looks nice for 25
mean_eggs <- mean_eggs %>%
  mutate(
    se_eggs = ifelse(n == 1, 0, se_eggs)
  )

# Keep species order
mean_eggs$species <- factor(mean_eggs$species, levels = c("L. sericata", "C. vicina"))

# Plot as line graph of means
p_total_eggs <- ggplot(mean_eggs, aes(x = temperature, y = mean_eggs, colour = species)) +
  
  # SE ribbon
  geom_ribbon(
    aes(
      ymin = mean_eggs - se_eggs,
      ymax = mean_eggs + se_eggs,
      fill = species,
      group = species
    ),
    alpha = 0.3,
    colour = NA
  )+
  
  geom_line(linewidth = 1.2) +
  geom_point(size = 6) +
  
  scale_color_manual(
    values = c("L. sericata" = "#C87F17", "C. vicina" = "blue"),
    name = NULL
  ) +
  
  scale_fill_manual(
    values = c("L. sericata" = "#C87F17", "C. vicina" = "blue"),
    guide = "none"
  ) +
  
  scale_x_continuous(breaks = sort(unique(mean_eggs$temperature))) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.05))) +
  
  labs(
    x = "Temperature (°C)",
    y = "Mean offspring laid"
  ) +
  
  common_theme_fig2
p_total_eggs

#-----------------------------
# PERFORMANCE INDEX
#-----------------------------

# read data
eggs_summary <- read.csv("fecun_sum_full.csv")

eggs_summary <- eggs_summary %>%
  mutate(
    species = ifelse(grepl("ls", treatment), "L. sericata", "C. vicina"),
    temperature = as.numeric(substr(treatment, 1, 2))
  )

# species order
eggs_summary$species <- factor(
  eggs_summary$species,
  levels = c("L. sericata", "C. vicina")
)

# summary stats (means only; no error bars plotted)
eggs_summary_summary <- eggs_summary %>%
  group_by(temperature, species) %>%
  summarise(
    avg_P = mean(performance_index, na.rm = TRUE),
    n     = n(),
    sd_P  = sd(performance_index, na.rm = TRUE),
    se_P  = sd_P / sqrt(n),
    .groups = "drop"
  )

eggs_summary_summary$species <- factor(
  eggs_summary_summary$species,
  levels = c("L. sericata", "C. vicina")
)

#-----------------------------
# Plot
#-----------------------------
alt <- ggplot() +
  
  # SE ribbon (behind everything)
  geom_ribbon(
    data = eggs_summary_summary,
    aes(
      x = temperature,
      ymin = avg_P - se_P,
      ymax = avg_P + se_P,
      fill = species,
      group = species
    ),
    alpha = 0.2,
    colour = NA
  ) +
  
  # raw replicate-level data
  geom_jitter(
    data = eggs_summary,
    aes(
      x = temperature,
      y = performance_index,
      colour = species
    ),
    width = 0.25,
    height = 0,
    size = 4,
    alpha = 0.3
  ) +
  
  # species means
  geom_point(
    data = eggs_summary_summary,
    aes(
      x = temperature,
      y = avg_P,
      colour = species
    ),
    size = 6
  ) +
  
  # mean lines
  geom_line(
    data = eggs_summary_summary,
    aes(
      x = temperature,
      y = avg_P,
      colour = species,
      group = species
    ),
    linewidth = 1.2
  ) +
  
  scale_color_manual(
    values = c(
      "C. vicina"   = "blue",
      "L. sericata" = "#C87F17"
    ),
    name = NULL
  ) +
  
  scale_fill_manual(
    values = c(
      "C. vicina"   = "blue",
      "L. sericata" = "#C87F17"
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = c(15, 19, 22, 25)
  ) +
  
  labs(
    x = "Temperature (°C)",
    y = "Performance index"
  ) +
  
  common_theme_fig2

alt

#so that 25 tapers nicely
eggs_summary_summary <- eggs_summary_summary %>%
  mutate(
    se_P = ifelse(n == 1, 0, se_P)
  )
#-----------------------------
# Plot no jitter
#-----------------------------
alt_nojitter <- ggplot() +
  
  # SE ribbon (behind everything)
  geom_ribbon(
    data = eggs_summary_summary,
    aes(
      x = temperature,
      ymin = avg_P - se_P,
      ymax = avg_P + se_P,
      fill = species,
      group = species
    ),
    alpha = 0.3,
    colour = NA
  ) +
  
  # species means
  geom_point(
    data = eggs_summary_summary,
    aes(
      x = temperature,
      y = avg_P,
      colour = species
    ),
    size = 6
  ) +
  
  # species mean lines
  geom_line(
    data = eggs_summary_summary,
    aes(
      x = temperature,
      y = avg_P,
      colour = species,
      group = species
    ),
    linewidth = 1.2,
    linetype = "solid"
  ) +
  
  scale_color_manual(
    values = c(
      "C. vicina"   = "blue",
      "L. sericata" = "#C87F17"
    ),
    name = NULL
  ) +
  
  scale_fill_manual(
    values = c(
      "C. vicina"   = "blue",
      "L. sericata" = "#C87F17"
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(
    breaks = c(15, 19, 22, 25)
  ) +
  
  labs(
    x = "Temperature (°C)",
    y = "Performance index"
  ) +
  
  common_theme_fig2

alt_nojitter

ggsave(
  filename = file.path(outdir, "performance_index_one_plot.png"),
  plot = alt, width = 7.5, height = 6, dpi = 600
)

library(cowplot)
library(ggplot2)

pC <- p_body_nojitter + theme(legend.position = "none", plot.margin = margin(10, 10, 10, 10))
pD <- p_total_eggs + theme(legend.position = "none", plot.margin = margin(10, 10, 10, 10)) 
pE <- PREDICTED + theme(legend.position = "none", plot.margin = margin(10, 10, 10, 10))

# --- 1) Remove legends from all plots ---
pC_noleg <- pC + theme(legend.position = "none", plot.margin = margin(10,10,10,10))
pD_noleg <- pD + theme(legend.position = "none", plot.margin = margin(10,10,10,10))
pE_noleg <- pE + theme(legend.position = "none", plot.margin = margin(10,10,10,10))

# --- 2) Extract legend from C (lines only, 2-row stacked) ---
legend_C <- cowplot::get_legend(
  pC +
    guides(colour = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        linetype = c("dotted", "solid", "dotted", "solid"),
        shape = 16,
        size = 5
      )
    )) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.key.width = unit(2, "lines"),
      legend.key.height = unit(0.5, "lines"),
      legend.spacing.x = unit(3, "lines"),
      legend.spacing.y = unit(1, "lines"),
      legend.text = element_text(size = 16, face = "italic")
    )
)

# --- 3) Extract legend from D (horizontal, shared for D+E) ---
legend_DE <- cowplot::get_legend(
  pD +
    guides(colour = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linetype = c("solid"),
        shape = 16,
        size = 5
      ),
      keywidth = unit(2, "lines"),
    )) +
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.justification = "center",
      legend.key.width = unit(3, "lines"),
      legend.key.height = unit(1.5, "lines"),
      legend.spacing.x = unit(2, "lines"),
      legend.text = element_text(size = 16, face = "italic")
    )
)

# --- 4) Stack D+E horizontally with shared legend ---
DE_panels <- cowplot::plot_grid(
  pD_noleg, pE_noleg,
  ncol = 2, align = "v", axis = "tblr"
)

DE_col <- cowplot::plot_grid(
  DE_panels,
  legend_DE,
  ncol = 1,
  rel_heights = c(1, 0.15)  # match C's total height
)

DE_col <- ggdraw(DE_col) +
  draw_plot_label(
    label = c("B", "C"),
    x = c(0, 0.5),  # adjust to center over D and E
    y = 1,
    fontface = "bold",
    size = 22,
    hjust = 0
  )

# --- 5) Stack C panel with its legend ---
C_col <- cowplot::plot_grid(
  pC_noleg,
  legend_C,
  ncol = 1,
  rel_heights = c(1, 0.15)  # smaller number -> legend closer to graph
)

C_col <- ggdraw(C_col) +
  draw_plot_label("A", x = 0, y = 1, fontface = "bold", size = 24, hjust = 0)

# --- 6) Combine C_col and DE_col into one row (panels same height) ---
fig_row <- cowplot::plot_grid(
  C_col, DE_col,
  ncol = 2,
  rel_widths = c(1, 2),
  align = "v",
  axis = "tblr"
)

# --- 7) Add labels A/B/C above each panel ---
fig_row <- ggdraw(fig_row) +
  draw_plot_label(
    label = c("A", "B", "C"),
    x = c(0.005, 0.34, 0.67),  # horizontal positions
    y = 1.05,                   # move higher above panels
    hjust = 0, vjust = 0,
    fontface = "bold",
    size = 22
  )

# --- 8) Display figure ---
fig_row

# --- 9) Save figure ---
ggsave(
  filename = file.path(outdir, "Figure_C_DE_E_two_legends.png"),
  plot = fig_row,
  dpi = 600
)


#############################################################################
########################################STATISTICAL ANALYSIS#################
#############################################################################
library(lme4)
library(car)
library(glmmTMB)

#make long
long_treat <- long %>%
  dplyr::filter(control == "treatment", initial > 0)


#####################
###EMERGENCE SUCCESS##
#######################
#set contrasts for type iii anova tests
options(contrasts = c("contr.sum", "contr.poly"))

#does temperature affect emergence success differently for L. sericata and C. vicina? 
#binomial GLMM (successes out of trials, tests species x temperature interaction)
#treat temperature as a factor
emergence_success_model <- glm(
  cbind(success, initial - success) ~ species * factor(temperature),
  family = binomial,
  data = long_treat
)
version(car)
Anova(emergence_success_model, type = "III")

#check residuals 
library(DHARMa)
sim_res <- simulateResiduals(fittedModel = emergence_success_model)
plot(sim_res)
testDispersion(sim_res)

#no evidence of overdispersion and data fits the model well

## emmeans to calculate emergence probability at 95% CI
library(emmeans)

emm <- emmeans(
  emergence_success_model,
  ~ species | factor(temperature)
)

summary(emm, type = "response")

#post-hoc
pairs(emm)  # species contrasts at each temperature

#plot to see if it matches the one i made earlier 
plot(emm, comparisons = TRUE, type = "response")

#############################
###EMERGENCE TIMING STATS###
##############################
#remove the 30 degrees from the model because no C. vicina emerged
df2 <- subset(long_treat, temperature != 30)

m_time <- lm(
  days_to_first ~ species * factor(temperature),
  data = long_treat
)

Anova(m_time, type = "III") 

#check residuals 
sim_res <- simulateResiduals(fittedModel = m_time)
plot(sim_res)
testDispersion(sim_res)

summary(m_time)

#emmeans for emergence time
emm_time <- emmeans(m_time, ~ species | factor(temperature))
pairs(emm_time)
summary(emm_time)


###BODY SIZE
#remove the 30 degrees from the model because no C. vicina emerged
df3 <- subset(body_treat, temperature != 30)
m_body <- lm(
  thorax_mm ~ species * factor(temperature) * sex,
  data = body_treat
)

#check residuals 
sim_res <- simulateResiduals(fittedModel = m_body)
plot(sim_res)
testDispersion(sim_res)

summary(m_body)
Anova(m_body, type = "III")


#emmeans for body size
emm_body <- emmeans(m_body, ~ species | factor(temperature) | sex)
pairs(emm_body)
summary(emm_body)

#eggs
m_eggs <- glmmTMB(
  lifetime_fecundity ~ species * factor(temperature),
  family = nbinom1(),
  data = eggs_summary
)

#check residuals 
sim_res <- simulateResiduals(fittedModel = m_eggs)
plot(sim_res)
testDispersion(sim_res)

summary(m_eggs)

#emmeans for body size
emm_eggs <- emmeans(m_eggs, ~ species | factor(temperature))
pairs(emm_eggs)
summary(emm_eggs)

####INCLUDING CONTROLS

# Only 22°C
long_22 <- long %>%
  filter(temperature == 22)

#remove the rows with zero (when we pivot to long form it included rows for the species not present in the single species replicates so need to remove them)
long_22 <- long_22 %>%
  filter(!is.na(days_to_first))

# Extract density from control replicates
long_22 <- long_22 %>%
  mutate(
    density = case_when(
      control == "treatment" ~ "treatment",
      str_detect(replicate, "10") ~ "n10",
      str_detect(replicate, "20") ~ "n20",
      TRUE ~ NA_character_
    )
  )

long_22 <- long_22 %>%
  mutate(
    group = case_when(
      species == "vicina"   & density == "treatment" ~ "C. vicina treatment",
      species == "vicina"   & density == "n10" ~ "C. vicina n = 10 density",
      species == "vicina"   & density == "n20" ~ "C. vicina n = 20 density",
      species == "sericata" & density == "treatment" ~ "L. sericata treatment",
      species == "sericata" & density == "n10" ~ "L. sericata n = 10 density",
      species == "sericata" & density == "n20" ~ "L. sericata n = 20 density"
    )
  )

long_22$group <- factor(
  long_22$group,
  levels = c(
    "C. vicina treatment",
    "C. vicina n = 10 density",
    "C. vicina n = 20 density",
    "L. sericata treatment",
    "L. sericata n = 10 density",
    "L. sericata n = 20 density"
  )
)

long_22 <- long_22 %>%
  mutate(
    colour_group = case_when(
      group == "C. vicina treatment" ~ "vicina_treat",
      group == "C. vicina n = 10 density" ~ "vicina_control_10",
      group == "C. vicina n = 20 density" ~ "vicina_control_20",
      group == "L. sericata treatment" ~ "sericata_treat",
      group == "L. sericata n = 10 density" ~ "sericata_control_10",
      group == "L. sericata n = 20 density" ~ "sericata_control_20"
    )
  )

summary_22 <- long_22 %>%
  group_by(group, colour_group) %>%
  summarise(
    mean_prop = mean(prop_success, na.rm = TRUE),
    n = n(),
    sd = sd(prop_success, na.rm = TRUE),
    se = sd / sqrt(n),
    ci_low = mean_prop - 1.96 * se,
    ci_high = mean_prop + 1.96 * se,
    .groups = "drop"
  )

summary_22$group <- factor(
  summary_22$group,
  levels = c(
    "C. vicina treatment",
    "C. vicina n = 10 density",
    "C. vicina n = 20 density",
    "L. sericata treatment",
    "L. sericata n = 10 density",
    "L. sericata n = 20 density"
  )
)

# ---- Optional: Simplify x-axis labels for clarity ----
summary_22 <- summary_22 %>%
  mutate(
    x_label = case_when(
      str_detect(group, "treatment") ~ "Mixed",
      str_detect(group, "n = 10") ~ "n=10",
      str_detect(group, "n = 20") ~ "n=20"
    ),
    species = case_when(
      str_detect(group, "vicina") ~ "C. vicina",
      str_detect(group, "sericata") ~ "L. sericata"
    )
  )

# ---- Plot ----
library(ggplot2)
library(scales)

p_density_facet <- ggplot(summary_22, aes(x = x_label, y = mean_prop, colour = species)) +
  # Error bars
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.7) +
  # Points
  geom_point(size = 6, stroke = 0.7) +
  # Y-axis
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  # Color mapping
  scale_color_manual(values = c("C. vicina" = "blue", "L. sericata" = "#C87F17")) +
  # Labels
  labs(x = NULL, y = "Emergence proportion (22 °C)") +
  # Facet by species side by side
  facet_wrap(~ species, scales = "free_x", nrow = 1) +
  # Theme
  theme_minimal(base_size = 20) +
  theme(
    # Axis lines
    axis.line = element_line(colour = "black", linewidth = 0),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),  
    axis.ticks.length = unit(0.25, "cm"),  
    
    # Axis text
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(colour = "black"),
    
    # Panel background and border
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    
    # Strip background and text
    strip.background = element_blank(),   # removes the box
    strip.text = element_text(face = "bold.italic", colour = "black", size = 20),
    
    # Space between panels
    panel.spacing = unit(1, "lines"),    # small space between facets
    
    # Remove grids and legend
    panel.grid = element_blank(),
    legend.position = "none"
  )
p_density_facet


#days to emergence 

summary_time_22 <- long_22 %>%
  group_by(group, colour_group) %>%
  summarise(
    mean_days = mean(days_to_first, na.rm = TRUE),
    n = n(),
    sd = sd(days_to_first, na.rm = TRUE),
    se = sd / sqrt(n),
    ci_low = mean_days - 1.96 * se,
    ci_high = mean_days + 1.96 * se,
    .groups = "drop"
  )

summary_time_22$group <- factor(
  summary_time_22$group,
  levels = levels(summary_22$group)
)

# ---- Add simplified x-axis labels and species column ----
summary_time_22 <- summary_time_22 %>%
  mutate(
    x_label = case_when(
      str_detect(group, "treatment") ~ "Mixed",
      str_detect(group, "n = 10") ~ "n=10",
      str_detect(group, "n = 20") ~ "n=20"
    ),
    species = case_when(
      str_detect(group, "vicina") ~ "C. vicina",
      str_detect(group, "sericata") ~ "L. sericata"
    )
  )

# ---- Plot ----
p_time_22_facet <- ggplot(summary_time_22, aes(x = x_label, y = mean_days, colour = species)) +
  
  # Error bars
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.7) +
  
  # Points
  geom_point(size = 6, stroke = 0.7) +
  
  # Y-axis
  scale_y_continuous(limits = c(15, max(summary_time_22$ci_high, na.rm = TRUE)),
                     labels = scales::number_format(accuracy = 1)) +
  
  # Color mapping
  scale_color_manual(values = c("C. vicina" = "blue", "L. sericata" = "#C87F17")) +
  
  # Labels
  labs(x = NULL, y = "Days to first emergence (22 °C)") +
  
  # Facet by species side by side
  facet_wrap(~ species, scales = "free_x", nrow = 1) +
  
  # Theme
  theme_minimal(base_size = 20) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(colour = "black"),
    
    # Axis lines
    axis.line = element_line(colour = "black", linewidth = 0),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.25, "cm"),
    
    # Panel background and border
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    
    # Facet strips (no box)
    strip.background = element_blank(),  # removes the box
    strip.text = element_text(face = "bold.italic", colour = "black", size = 20),
    
    # Space between panels
    panel.spacing = unit(1, "lines"),
    
    # Remove grids and legend
    panel.grid = element_blank(),
    legend.position = "none"
  )

p_time_22_facet


library(patchwork)
plot_together <- (p_density_facet + p_time_22_facet) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 26)
  )
plot_together
ggsave(
  filename = file.path("fig.S1.png"),
  plot = plot_together,
  dpi = 600
)


###########FIG S2###################
# =========================
# LOAD DATA
# =========================
df22 <- read.csv("fecun_sum_22.csv")

# =========================
# CLEAN + GROUP VARIABLES
# =========================
df22_clean <- df22 %>%
  mutate(
    species = case_when(
      str_detect(treatment, "ls|LS") ~ "L. sericata",
      str_detect(treatment, "cv|CV") ~ "C. vicina"
    ),
    
    group = case_when(
      treatment %in% c("22ls","22cv") ~ "treatment",
      str_detect(treatment, "10") ~ "n10",
      str_detect(treatment, "20") ~ "n20"
    ),
    
    x_group = case_when(
      species == "C. vicina" & group == "treatment" ~ "C. vicina treatment",
      species == "L. sericata" & group == "treatment" ~ "L. sericata treatment",
      species == "C. vicina" & group == "n10" ~ "C. vicina n = 10 density",
      species == "C. vicina" & group == "n20" ~ "C. vicina n = 20 density",
      species == "L. sericata" & group == "n10" ~ "L. sericata n = 10 density",
      species == "L. sericata" & group == "n20" ~ "L. sericata n = 20 density"
    )
  )

# =========================
# PERFORMANCE SUMMARY
# =========================
summary_perf_22 <- df22_clean %>%
  group_by(x_group) %>%
  summarise(
    mean_perf = mean(performance_index, na.rm = TRUE),
    sd = sd(performance_index, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_low = mean_perf - 1.96 * se,
    ci_high = mean_perf + 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(
    x_label = case_when(
      str_detect(x_group, "treatment") ~ "Mixed",
      str_detect(x_group, "n = 10") ~ "n=10",
      str_detect(x_group, "n = 20") ~ "n=20"
    ),
    species = case_when(
      str_detect(x_group, "vicina") ~ "C. vicina",
      str_detect(x_group, "sericata") ~ "L. sericata"
    )
  )

# =========================
# EGGS SUMMARY
# =========================
summary_eggs_22 <- df22_clean %>%
  group_by(x_group) %>%
  summarise(
    mean_eggs = mean(lifetime_fecundity, na.rm = TRUE),
    sd = sd(lifetime_fecundity, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    ci_low = mean_eggs - 1.96 * se,
    ci_high = mean_eggs + 1.96 * se,
    .groups = "drop"
  ) %>%
  mutate(
    x_label = case_when(
      str_detect(x_group, "treatment") ~ "Mixed",
      str_detect(x_group, "n = 10") ~ "n=10",
      str_detect(x_group, "n = 20") ~ "n=20"
    ),
    species = case_when(
      str_detect(x_group, "vicina") ~ "C. vicina",
      str_detect(x_group, "sericata") ~ "L. sericata"
    )
  )


# =========================
# COMMON THEME
# =========================
common_theme <- theme_minimal(base_size = 16) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    
    axis.line = element_line(colour = "black", linewidth = 0),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.25, "cm"),
    
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    
    strip.background = element_blank(),
    strip.text = element_text(face = "bold.italic", size = 16),
    
    panel.spacing = unit(1, "lines"),
    panel.grid = element_blank(),
    
    legend.position = "none"
  )

# =========================
# PLOTS
# =========================

# Performance
p_perf_22 <- ggplot(summary_perf_22,
                    aes(x = x_label, y = mean_perf, colour = species)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.6) +
  geom_point(size = 5) +
  scale_color_manual(values = c("C. vicina" = "blue", "L. sericata" = "#C87F17")) +
  labs(x = NULL, y = "Performance index (22 °C)") +
  facet_wrap(~ species, scales = "free_x", nrow = 1) +
  common_theme

# Eggs
p_eggs_22 <- ggplot(summary_eggs_22,
                    aes(x = x_label, y = mean_eggs, colour = species)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.15, linewidth = 0.6) +
  geom_point(size = 5) +
  scale_color_manual(values = c("C. vicina" = "blue", "L. sericata" = "#C87F17")) +
  labs(x = NULL, y = "Mean lifetime eggs/larvae laid (22 °C)") +
  facet_wrap(~ species, scales = "free_x", nrow = 1) +
  common_theme


p_perf_22
p_eggs_22

#body size is a bit weird
# 22 only + controls
body_22 <- body_long %>%
  dplyr::filter(tolower(temperature) == "22", !is.na(thorax_mm))

#extract density and create new column to summarise by
body_22 <- body_22 %>%
  mutate(
    density = case_when(
      str_detect(replicate, "10") ~ "n=10",
      str_detect(replicate, "20") ~ "n=20",
      TRUE ~ NA_character_
    )
  )
body_22 <- body_22 %>%
  mutate(
    group_label = case_when(
      control == "treatment" ~ paste0(species, "-treatment"),
      control == "control"   ~ paste0(density, "-", species),
      TRUE ~ NA_character_
    )
  )

##### SUMMARY (means only) #####
summary_body_22 <- body_22 %>%
  group_by(species, sex, group_label) %>%
  summarise(
    mean_mm = mean(thorax_mm, na.rm = TRUE),
    n = n(),
    sd = sd(thorax_mm, na.rm = TRUE),
    se = sd / sqrt(n),
    ci_low = mean_mm - 1.96 * se,
    ci_high = mean_mm + 1.96 * se,
    .groups = "drop"
  )

#####BODY SIZE
dodge <- position_dodge(width = 0.5)

# Create a new factor for facets with labels as expressions
summary_body_22 <- summary_body_22 %>%
  mutate(
    species_facet = factor(
      species,
      levels = c("vicina", "sericata"),  # Flip order
      labels = c(expression("C. vicina"), 
                 expression("L. sericata"))
    )
  )

#make the x-axis look better
summary_body_22 <- summary_body_22 %>%
  mutate(
    x_label = case_when(
      group_label == "sericata-treatment" ~ "Mixed",
      group_label == "vicina-treatment" ~ "Mixed",
      str_detect(group_label, "n=10") ~ "n = 10",
      str_detect(group_label, "n=20") ~ "n = 20",
      TRUE ~ NA_character_
    )
  )

# Plot
p_body_22 <- ggplot(
  summary_body_22,
  aes(
    x = x_label,
    y = mean_mm,
    colour = species,
    shape  = sex,
    group  = interaction(species, sex)
  )
) +
  
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.15,
    linewidth = 0.6,
    position = dodge,
    show.legend = FALSE
  ) +
  
  geom_point(
    size = 5,
    stroke = 0.7,
    position = dodge,
    show.legend = TRUE
  ) +
  
  scale_color_manual(values = c(sericata = "#C87F17", vicina = "blue"), guide = "none") +
  
  scale_shape_manual(
    values = c(male = 16, female = 17),
    labels = c(male = "Male", female = "Female"),  # Capitalized
    name = NULL,
    guide = guide_legend(
      keyheight = unit(2, "lines"),   # height of each legend item
      default.unit = "lines",
      label.position = "right",
      title.position = "top",
      label.theme = element_text(size = 16)
    )
  ) +
  
  facet_wrap(
    ~ species_facet,
    nrow = 1
  ) +
  
  labs(
    x = NULL,
    y = "Thorax length (mm; 22 °C)"
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 16, face = "bold.italic"),
    
    axis.line = element_line(linewidth = 0),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 1.2
    ),
    
    legend.position = "right"
  ) 
p_body_22



library(patchwork)
library(cowplot)

# 1. Extract the sex legend from body plot
legend_body <- cowplot::get_legend(
  p_body_22 +
    guides(shape = guide_legend(title = NULL)) +
    theme(legend.position = "bottom",
          legend.justification = c(0.12, 1.5),
          legend.direction = "horizontal",
          legend.spacing.x = unit(1, "cm"),
          legend.text = element_text(size = 16))
)

# 2. Remove legends from all plots
p_body_22_clean <- p_body_22 + theme(legend.position = "none")
p_eggs_22_clean <- p_eggs_22 + theme(legend.position = "none")
p_perf_22_clean <- PREDICTED + theme(legend.position = "none")

# 3. Combine the 3 plots side by side
plots_row <- p_body_22_clean + p_eggs_22_clean + p_perf_22_clean + 
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 26))

# 4. Combine the row of plots with the legend underneath body plot
final_controls <- cowplot::plot_grid(
  plots_row,       # top row with 3 plots
  legend_body,     # bottom row with legend
  ncol = 1,
  rel_heights = c(1, 0.1)  # adjust height for legend row
)

final_controls

ggsave(
  filename = file.path("fig.S2.png"),
  plot = final_controls,
  dpi = 600
)


#TWEEDIE MODEL FOR EGGS
# Load libraries
library(tidyverse)
library(tweedie)
library(statmod)
library(car)
# -----------------------------
# 1. DATA PREP
# -----------------------------
total_eggs_by_treatment <- total_eggs_by_treatment |>
  mutate(temperature = factor(temperature))

# -----------------------------
# 2. OPTIMAL TWEEDIE PARAMETER
# -----------------------------
library(tweedie)

profile <- tweedie.profile(
  total_eggs ~ species * temperature,
  data = total_eggs_by_treatment,
  p.vec = seq(1.1, 1.9, by = 0.05),
  link.power = 0
)

best_p <- profile$p.max
cat("Optimal power parameter:", best_p, "\n")

# -----------------------------
# 3. FIT FINAL MODEL
# -----------------------------
model_tweedie <- glm(
  total_eggs ~ species * temperature,
  data = total_eggs_by_treatment,
  family = tweedie(var.power = best_p, link.power = 0)
)

summary(model_tweedie)

# Type III ANOVA
car::Anova(model_tweedie, type = 3)

# -----------------------------
# 4. PREDICTION GRID
# -----------------------------
pred_grid <- expand.grid(
  species = unique(total_eggs_by_treatment$species),
  temperature = sort(unique(total_eggs_by_treatment$temperature))
)

# -----------------------------
# 5. PREDICTIONS (CORRECT SCALE)
# -----------------------------
pred <- predict(
  model_tweedie,
  newdata = pred_grid,
  type = "link",
  se.fit = TRUE
)

# Back-transform
pred_grid$fit <- exp(pred$fit)

# Correct 95% CI (IMPORTANT FIX)
pred_grid$lower <- exp(pred$fit - 1.96 * pred$se.fit)
pred_grid$upper <- exp(pred$fit + 1.96 * pred$se.fit)

# -----------------------------
# 6. PLOT
# -----------------------------
ggplot(pred_grid,
       aes(x = as.numeric(as.character(temperature)),
           y = fit,
           group = species)) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = species),
    alpha = 0.15,
    colour = NA
  ) +
  
  geom_line(aes(colour = species), linewidth = 1.2) +
  geom_point(aes(colour = species), size = 4) +
  
  scale_colour_manual(
    values = c(
      "calliphora vicina" = "blue",
      "lucilia sericata"  = "#C87F17"
    ),
    labels = c(
      "calliphora vicina" = expression(italic("C. vicina")),
      "lucilia sericata"  = expression(italic("L. sericata"))
    ),
    name = NULL
  ) +
  
  scale_fill_manual(
    values = c(
      "calliphora vicina" = "blue",
      "lucilia sericata"  = "#C87F17"
    ),
    guide = "none"
  ) +
  
  scale_x_continuous(breaks = c(15, 19, 22, 25)) +
  
  labs(
    x = "Temperature (°C)",
    y = "Predicted total eggs"
  ) +
  
  common_theme_fig2 +
  theme(
    legend.background = element_blank(),
    legend.key = element_blank()
  )

# -----------------------------
# 7. EMMEANS (CLEAN VERSION)
# -----------------------------
library(emmeans)

emm <- emmeans(model_tweedie, ~ species * temperature)

# Regrid to response scale properly
emm_resp <- regrid(emm)

emm_resp

# Pairwise comparisons (on response scale)
pairs(emmeans(model_tweedie, ~ species | temperature) |> regrid())


##########################################TO GET PREDICTED PERFORMANCE 
# 1. Load data
eggs_long <- read_csv("eggs.csv") |>
  rename(count = eggs, ID = cage_id) |>
  mutate(
    Species     = factor(species_code),
    Temperature = factor(temperature),
    Treatment   = as.character(Treatment),
    IsControl   = str_detect(Treatment, regex("control", ignore_case = TRUE))
  )

emerge <- read_csv("comp_data.csv") 

# 2. Reshape and clean emergence data
emerge_long <- emerge %>%
  pivot_longer(
    cols = c(final_number_observed_sericata, final_number_observed_vicina),
    names_to = "species",
    values_to = "emerged"
  ) %>%
  mutate(
    species = dplyr::recode(species,
                            final_number_observed_sericata = "sericata",
                            final_number_observed_vicina   = "vicina"),
    species = factor(species)
  )

emerge_long <- emerge_long %>%
  mutate(
    initial_larvae = case_when(
      str_detect(replicate, "^LS10") & species == "sericata" ~ 10,
      str_detect(replicate, "^LS20") & species == "sericata" ~ 20,
      str_detect(replicate, "^CV10") & species == "vicina"   ~ 10,
      str_detect(replicate, "^CV20") & species == "vicina"   ~ 20,
      TRUE ~ 10   
    ),
    failed = initial_larvae - emerged
  )

# 3. Filter to competition only and process dates
dat <- emerge_long %>%
  filter(control == "treatment")

date_cols <- c("date_offered_meat",
               "date_larvae_transferred_to_treatment",
               "date_first_pupa",
               "date_first_emerged",
               "date_second_emerged")

dat <- dat |>
  mutate(across(
    any_of(date_cols),
    ~ parse_date_time(.x, orders = c("dmy", "dmY")) |> as.Date()
  ))

dat <- dat |>
  mutate(
    species                = tolower(trimws(species)),
    species_first_emerged    = tolower(trimws(as.character(species_first_emerged))),
    species_second_emerged   = tolower(trimws(as.character(species_second_emerged))),
    species_emergence_date = case_when(
      species_first_emerged  == species ~ date_first_emerged,
      species_second_emerged == species ~ date_second_emerged,
      species_first_emerged  == "both"  ~ date_first_emerged,
      species_second_emerged == "both"  ~ date_second_emerged,
      TRUE                              ~ as.Date(NA)
    ),
    days_to_emergence = as.numeric(
      difftime(species_emergence_date, date_larvae_transferred_to_treatment, units = "days")
    )
  )

temp_order <- sort(unique(dat$temperature))
dat <- dat |>
  mutate(
    temperature_f = factor(temperature, levels = temp_order, ordered = TRUE),
    species_f     = factor(species)
  )

# 4. Strip down dataset to clean variables
cleaned <- dat |>
  mutate(replicate = as.integer(str_extract(replicate, "\\d+$"))) |>
  select(temperature, species_f, replicate, days_to_emergence, emerged, initial_larvae, failed)

# 5. Integrate Egg Data for Performance Calculations
total_eggs_by_treatment <- eggs_long |>
  group_by(Treatment, temperature, species) |>
  summarise(total_eggs = sum(count, na.rm = TRUE), .groups = "drop")

mean_eggs_by_temperature <- total_eggs_by_treatment |>
  group_by(temperature, species) |>
  summarise(mean_eggs = mean(total_eggs), .groups = "drop")

mean_eggs_by_temperature$temperature <- as.factor(mean_eggs_by_temperature$temperature)

cleaned_final_all <- cleaned |>
  mutate(species_key = case_when(
    species_f == "sericata" ~ "lucilia sericata",
    species_f == "vicina"   ~ "calliphora vicina"
  )) |>
  left_join(mean_eggs_by_temperature, by = c("temperature", "species_key" = "species")) |>
  select(-species_key)

cleaned_final_all <- cleaned_final_all |>
  mutate(
    S = emerged / initial_larvae,
    performance = log((mean_eggs * S) / days_to_emergence),
    performance = ifelse(is.na(performance) | is.infinite(performance), 0, performance)
  ) |> 
  filter(temperature != 30) |> 
  mutate(performance = abs(performance))

# 6. Fit GLM Model and Extract Predictions
model <- glm(performance ~ species_f * factor(temperature), data = cleaned_final_all, family = gaussian)
Anova(model, type=3)
summary(model)

#check residuals 
library(DHARMa)
sim_res <- simulateResiduals(fittedModel = model)
plot(sim_res)
testDispersion(sim_res)

pred_data <- expand.grid(
  species_f = unique(cleaned_final_all$species_f),
  temperature = unique(cleaned_final_all$temperature)
)

preds <- predict(model, newdata = pred_data, se.fit = TRUE)

pred_data <- pred_data |>
  mutate(
    predicted = preds$fit,
    se = preds$se.fit,
    lower = predicted - 1.96 * se,
    upper = predicted + 1.96 * se,
    temperature = factor(temperature)
  )

# 7. Generate Final Plot
# Note: Ensure 'common_theme_fig2' is already defined in your global environment.
PREDICTED <- ggplot() +
  # Confidence interval ribbon
  geom_ribbon(
    data = pred_data,
    aes(
      x = temperature,
      ymin = lower,
      ymax = upper,
      fill = species_f,
      group = species_f
    ),
    alpha = 0.3,
    colour = NA
  ) +
  # Predicted mean lines
  geom_line(
    data = pred_data,
    aes(
      x = temperature,
      y = predicted,
      colour = species_f,
      group = species_f
    ),
    linewidth = 1.2
  ) +
  # Predicted means
  geom_point(
    data = pred_data,
    aes(
      x = temperature,
      y = predicted,
      colour = species_f
    ),
    size = 6
  ) +
  scale_colour_manual(
    values = c(
      "sericata" = "#C87F17",
      "vicina"   = "blue"
    ),
    labels = c(
      "sericata" = "L. sericata",
      "vicina"   = "C. vicina"
    ),
    name = NULL
  ) +
  scale_fill_manual(
    values = c(
      "sericata" = "#C87F17",
      "vicina"   = "blue"
    ),
    guide = "none"
  ) +
  scale_x_discrete(
    limits = c("15", "19", "22", "25")
  ) +
  labs(
    x = "Temperature (°C)",
    y = "Predicted performance; ln(F*S/D)"
  ) +
  common_theme_fig2

# Render plot
PREDICTED
ggsave("predicted.png", plot = PREDICTED, dpi = 500)


#CONTROLS################################################################
#s1

# Only 22°C
long_22 <- long %>%
  filter(temperature == 22)

# Extract density from control replicates
long_22 <- long_22 %>%
  mutate(
    density = case_when(
      control == "treatment" ~ "treatment",
      str_detect(replicate, "10") ~ "n10",
      str_detect(replicate, "20") ~ "n20",
      TRUE ~ NA_character_
    )
  )

long_22 <- long_22 %>%
  mutate(
    group = case_when(
      species == "vicina"   & density == "treatment" ~ "C. vicina treatment",
      species == "vicina"   & density == "n10" ~ "C. vicina n = 10 density",
      species == "vicina"   & density == "n20" ~ "C. vicina n = 20 density",
      species == "sericata" & density == "treatment" ~ "L. sericata treatment",
      species == "sericata" & density == "n10" ~ "L. sericata n = 10 density",
      species == "sericata" & density == "n20" ~ "L. sericata n = 20 density"
    )
  )

long_22$group <- factor(
  long_22$group,
  levels = c(
    "C. vicina treatment",
    "C. vicina n = 10 density",
    "C. vicina n = 20 density",
    "L. sericata treatment",
    "L. sericata n = 10 density",
    "L. sericata n = 20 density"
  )
)

long_22 <- long_22 %>%
  mutate(
    colour_group = case_when(
      group == "C. vicina treatment" ~ "vicina_treat",
      group == "C. vicina n = 10 density" ~ "vicina_control_10",
      group == "C. vicina n = 20 density" ~ "vicina_control_20",
      group == "L. sericata treatment" ~ "sericata_treat",
      group == "L. sericata n = 10 density" ~ "sericata_control_10",
      group == "L. sericata n = 20 density" ~ "sericata_control_20"
    )
  )

summary_22 <- long_22 %>%
  group_by(group, colour_group) %>%
  summarise(
    mean_prop = mean(prop_success, na.rm = TRUE),
    n = n(),
    sd = sd(prop_success, na.rm = TRUE),
    se = sd / sqrt(n),
    ci_low = mean_prop - 1.96 * se,
    ci_high = mean_prop + 1.96 * se,
    .groups = "drop"
  )

summary_22$group <- factor(
  summary_22$group,
  levels = c(
    "C. vicina treatment",
    "C. vicina n = 10 density",
    "C. vicina n = 20 density",
    "L. sericata treatment",
    "L. sericata n = 10 density",
    "L. sericata n = 20 density"
  )
)

# ---- Optional: Simplify x-axis labels for clarity ----
summary_22 <- summary_22 %>%
  mutate(
    x_label = case_when(
      str_detect(group, "treatment") ~ "Mixed",
      str_detect(group, "n = 10") ~ "n=10",
      str_detect(group, "n = 20") ~ "n=20"
    ),
    species = case_when(
      str_detect(group, "vicina") ~ "C. vicina",
      str_detect(group, "sericata") ~ "L. sericata"
    )
  )

# ---- Plot ----
library(ggplot2)
library(scales)

p_density_facet <- ggplot(summary_22, aes(x = x_label, y = mean_prop, colour = species)) +
  # Error bars
  geom_errorbar(aes(ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.15, linewidth = 0.7) +
  # Points
  geom_point(size = 6, stroke = 0.7) +
  # Y-axis
  scale_y_continuous(limits = c(0,1), labels = scales::percent_format(accuracy = 1)) +
  # Color mapping
  scale_color_manual(values = c("C. vicina" = "blue", "L. sericata" = "#C87F17")) +
  # Labels
  labs(x = NULL, y = "Emergence proportion (22 °C)") +
  # Facet by species side by side
  facet_wrap(~ species, scales = "free_x", nrow = 1) +
  # Theme
  theme_minimal(base_size = 20) +
  theme(
    # Axis lines
    axis.line = element_line(colour = "black", linewidth = 0),
    axis.ticks = element_line(colour = "black", linewidth = 0.5),  
    axis.ticks.length = unit(0.25, "cm"),  
    
    # Axis text
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(colour = "black"),
    
    # Panel background and border
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    
    # Strip background and text
    strip.background = element_blank(),   # removes the box
    strip.text = element_text(face = "bold.italic", colour = "black", size = 20),
    
    # Space between panels
    panel.spacing = unit(1, "lines"),    # small space between facets
    
    # Remove grids and legend
    panel.grid = element_blank(),
    legend.position = "none"
  )
p_density_facet


#days to emergence 

summary_time_22 <- long_22 %>%
  group_by(group, colour_group) %>%
  summarise(
    mean_days = mean(days_to_first, na.rm = TRUE),
    n = n(),
    sd = sd(days_to_first, na.rm = TRUE),
    se = sd / sqrt(n),
    ci_low = mean_days - 1.96 * se,
    ci_high = mean_days + 1.96 * se,
    .groups = "drop"
  )

summary_time_22$group <- factor(
  summary_time_22$group,
  levels = levels(summary_22$group)
)

# ---- Add simplified x-axis labels and species column ----
summary_time_22 <- summary_time_22 %>%
  mutate(
    x_label = case_when(
      str_detect(group, "treatment") ~ "Mixed",
      str_detect(group, "n = 10") ~ "n=10",
      str_detect(group, "n = 20") ~ "n=20"
    ),
    species = case_when(
      str_detect(group, "vicina") ~ "C. vicina",
      str_detect(group, "sericata") ~ "L. sericata"
    )
  )

# ---- Plot ----
p_time_22_facet <- ggplot(summary_time_22, aes(x = x_label, y = mean_days, colour = species)) +
  
  # Error bars
  geom_errorbar(aes(ymin = mean_days - se, ymax = mean_days + se),
                width = 0.15, linewidth = 0.7) +
  
  # Points
  geom_point(size = 6, stroke = 0.7) +
  
  # Y-axis
  scale_y_continuous(limits = c(15, max(summary_time_22$ci_high, na.rm = TRUE)),
                     labels = scales::number_format(accuracy = 1)) +
  
  # Color mapping
  scale_color_manual(values = c("C. vicina" = "blue", "L. sericata" = "#C87F17")) +
  
  # Labels
  labs(x = NULL, y = "Days to first emergence (22 °C)") +
  
  # Facet by species side by side
  facet_wrap(~ species, scales = "free_x", nrow = 1) +
  
  # Theme
  theme_minimal(base_size = 20) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(colour = "black"),
    
    # Axis lines
    axis.line = element_line(colour = "black", linewidth = 0),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.ticks.length = unit(0.25, "cm"),
    
    # Panel background and border
    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    
    # Facet strips (no box)
    strip.background = element_blank(),  # removes the box
    strip.text = element_text(face = "bold.italic", colour = "black", size = 20),
    
    # Space between panels
    panel.spacing = unit(1, "lines"),
    
    # Remove grids and legend
    panel.grid = element_blank(),
    legend.position = "none"
  )

p_time_22_facet


library(patchwork)
plot_together <- (p_density_facet / p_time_22_facet) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(face = "bold", size = 26)
  )
plot_together
ggsave(
  filename = file.path("contrl_emerge.png"),
  plot = plot_together,
  dpi = 600
)


#FIGS2

# ==============================================================================
# LOAD DATA
# ==============================================================================
df22 <- read.csv("fecun_sum_22.csv")

# ==============================================================================
# CLEAN + GROUP VARIABLES (Fecundity/Performance Data)
# ==============================================================================
df22_clean <- df22 %>%
  mutate(
    species = case_when(
      str_detect(treatment, "ls|LS") ~ "L. sericata",
      str_detect(treatment, "cv|CV") ~ "C. vicina"
    ),
    group = case_when(
      treatment %in% c("22ls","22cv") ~ "treatment",
      str_detect(treatment, "10") ~ "n10",
      str_detect(treatment, "20") ~ "n20"
    ),
    x_group = case_when(
      species == "C. vicina" & group == "treatment" ~ "C. vicina treatment",
      species == "L. sericata" & group == "treatment" ~ "L. sericata treatment",
      species == "C. vicina" & group == "n10" ~ "C. vicina n = 10 density",
      species == "C. vicina" & group == "n20" ~ "C. vicina n = 20 density",
      species == "L. sericata" & group == "n10" ~ "L. sericata n = 10 density",
      species == "L. sericata" & group == "n20" ~ "L. sericata n = 20 density"
    )
  )

# ==============================================================================
# SUMMARIES
# ==============================================================================
summary_perf_22 <- df22_clean %>%
  group_by(x_group) %>%
  summarise(
    mean_perf = mean(performance_index, na.rm = TRUE),
    sd = sd(performance_index, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    x_label = case_when(
      str_detect(x_group, "treatment") ~ "Mixed",
      str_detect(x_group, "n = 10") ~ "n = 10",
      str_detect(x_group, "n = 20") ~ "n = 20"
    ),
    species = case_when(
      str_detect(x_group, "vicina") ~ "C. vicina",
      str_detect(x_group, "sericata") ~ "L. sericata"
    ),
    x_label = fct_relevel(x_label, "Mixed", "n = 10", "n = 20"),
    species = fct_relevel(species, "C. vicina", "L. sericata")
  )

summary_eggs_22 <- df22_clean %>%
  group_by(x_group) %>%
  summarise(
    mean_eggs = mean(lifetime_fecundity, na.rm = TRUE),
    sd = sd(lifetime_fecundity, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(
    x_label = case_when(
      str_detect(x_group, "treatment") ~ "Mixed",
      str_detect(x_group, "n = 10") ~ "n = 10",
      str_detect(x_group, "n = 20") ~ "n = 20"
    ),
    species = case_when(
      str_detect(x_group, "vicina") ~ "C. vicina",
      str_detect(x_group, "sericata") ~ "L. sericata"
    ),
    x_label = fct_relevel(x_label, "Mixed", "n = 10", "n = 20"),
    species = fct_relevel(species, "C. vicina", "L. sericata")
  )

# ==============================================================================
# CLEAN + GROUP VARIABLES (Body Size Data)
# ==============================================================================
dat <- read.csv("comp_data.csv", na.strings = c("N/A","NA",""))
names(dat) <- tolower(gsub("\\.", "_", names(dat))) 

body_long <- dat %>%
  transmute(
    replicate, temperature, control = factor(control),
    lucilia_male   = as.numeric(average_body_length_males_lucilia),
    lucilia_female = as.numeric(average_body_length_females_lucilia),
    vicina_male    = as.numeric(average_body_length_males_vicina),
    vicina_female  = as.numeric(average_body_length_females_vicina)
  ) %>%
  pivot_longer(
    cols = c(lucilia_male, lucilia_female, vicina_male, vicina_female),
    names_to = "key", values_to = "thorax_mm"
  ) %>%
  mutate(
    species = case_when(
      str_detect(key, "^lucilia_") ~ "sericata",
      str_detect(key, "^vicina_")  ~ "vicina",
      TRUE ~ NA_character_
    ),
    sex = if_else(str_detect(key, "_male$"), "male", "female"),
    species = factor(species, levels = c("sericata","vicina")),
    sex     = factor(sex, levels = c("male","female"))
  ) %>%
  dplyr::select(-key)

body_22 <- body_long %>%
  dplyr::filter(tolower(temperature) == "22", !is.na(thorax_mm))

body_22_plot_data <- body_22 %>%
  mutate(
    species_label = case_when(
      species == "sericata" ~ "L. sericata",
      species == "vicina"   ~ "C. vicina",
      TRUE                  ~ as.character(species)
    ),
    density_group = case_when(
      control == "treatment"      ~ "Mixed",   
      str_detect(replicate, "10") ~ "n = 10",  
      str_detect(replicate, "20") ~ "n = 20",  
      TRUE                        ~ "Unknown"  
    ),
    density_group = fct_relevel(density_group, "Mixed", "n = 10", "n = 20"),
    species_label = fct_relevel(species_label, "C. vicina", "L. sericata"),
    sex = str_to_title(sex)
  ) %>%
  group_by(species_label, density_group, sex) %>%
  summarise(
    mean_thorax = mean(thorax_mm, na.rm = TRUE),
    sd_thorax   = sd(thorax_mm, na.rm = TRUE),
    n           = n(),
    se_thorax   = sd_thorax / sqrt(n),
    .groups     = "drop"
  )

dodge <- position_dodge(width = 0.5)

# ==============================================================================
# BASELINE INDIVIDUAL PLOTS WITH LOCAL LEGENDS
# ==============================================================================

# A. Body Size Plot
p_body_22 <- ggplot(
  body_22_plot_data,
  aes(
    x      = density_group,                  
    y      = mean_thorax,                  
    colour = species_label,                 
    shape  = sex,
    group  = interaction(species_label, sex)
  )
) +
  geom_errorbar(aes(ymin = mean_thorax - se_thorax, ymax = mean_thorax + se_thorax), width = 0.15, linewidth = 1.2, position = dodge) +
  geom_point(size = 5.5, stroke = 0.8, position = dodge) +
  scale_color_manual(values = c("L. sericata" = "#C87F17", "C. vicina" = "blue"), 
                     labels = c("C. vicina" = "*C. vicina*", "L. sericata" = "*L. sericata*"), 
                     name = NULL) +
  scale_shape_manual(values = c("Male" = 15, "Female" = 17), name = NULL) +
  facet_wrap(~ species_label, ncol = 2, scales = "fixed") +
  labs(x = NULL, y = "Thorax length (mm; 22 °C)", tag = "A") +
  guides(
    shape = guide_legend(order = 1),
    colour = guide_legend(order = 2, override.aes = list(shape = 16, size = 6))
  )

# B. Eggs Plot (Turned on legends, configured styling parameters matching plot A)
p_eggs_22 <- ggplot(
  summary_eggs_22,
  aes(x = x_label, y = mean_eggs, colour = species, group = species)
) +
  geom_errorbar(aes(ymin = mean_eggs - se, ymax = mean_eggs + se), width = 0.15, linewidth = 1.2) +
  geom_point(size = 5.5) +
  scale_color_manual(values = c("C. vicina" = "blue", "L. sericata" = "#C87F17"), 
                     labels = c("C. vicina" = "*C. vicina*", "L. sericata" = "*L. sericata*"),
                     name = NULL) +
  facet_wrap(~ species, ncol = 2, scales = "fixed") +
  labs(x = NULL, y = "Mean lifetime eggs laid (22 °C)", tag = "B") +
  guides(
    colour = guide_legend(override.aes = list(shape = 16, size = 6))
  )

# C. Performance Index Plot (Turned on legends, configured styling parameters matching plot A)
p_perf_22 <- ggplot(
  summary_perf_22,
  aes(x = x_label, y = mean_perf, colour = species, group = species)
) +
  geom_errorbar(aes(ymin = mean_perf - se, ymax = mean_perf + se), width = 0.15, linewidth = 1.2) +
  geom_point(size = 5.5) +
  scale_color_manual(values = c("C. vicina" = "blue", "L. sericata" = "#C87F17"), 
                     labels = c("C. vicina" = "*C. vicina*", "L. sericata" = "*L. sericata*"),
                     name = NULL) +
  facet_wrap(~ species, ncol = 2, scales = "fixed") +
  labs(x = NULL, y = "Performance index (22 °C)", tag = "C") +
  guides(
    colour = guide_legend(override.aes = list(shape = 16, size = 6))
  )

# ==============================================================================
# UNIFIED THEME CONFIGURATION (Large, Bold, & Readable)
# ==============================================================================
unified_theme <- theme_classic(base_size = 20) +
  theme(
    # Bounding Boxes & Lines
    axis.line = element_line(linewidth = 0),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
    axis.ticks = element_line(colour = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.25, "cm"),
    
    # Text Dimensions
    axis.text.x = element_text(size = 16, color = "black", margin = margin(t = 8)),
    axis.text.y = element_text(size = 16, color = "black", margin = margin(r = 8)),
    axis.title.y = element_text(size = 18, color = "black", margin = margin(r = 12)),
    
    # Facet Label Formatting
    strip.background = element_blank(),
    strip.text = element_text(size = 18, face = "bold.italic", color = "black", margin = margin(b = 10)),
    
    # Subplot Tag Sizing (A, B, C)
    plot.tag = element_text(face = "bold", size = 24, margin = margin(b = 13)),
    
    # Grid Polish
    panel.spacing = unit(1.5, "lines"),
    
    # Legend Uniformity Options (Now targets individual layout tracking base rules)
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.direction = "horizontal",
    legend.margin = margin(t = 5, b = 5, r = 0, l = 0),
    legend.text = element_markdown(size = 18), 
    legend.spacing.x = unit(0.3, "cm"),
    legend.key.height = unit(1, "lines")
  )

# ==============================================================================
# GRID ASSEMBLY VIA PATCHWORK
# ==============================================================================

# Row 1: Contains Subplots A and B side-by-side
row1 <- p_body_22 + p_eggs_22 + plot_layout(ncol = 2)

# Row 2: Centers Subplot C by putting it inside a 3-column split grid
row2 <- plot_spacer() + p_perf_22 + plot_spacer() + 
  plot_layout(ncol = 3, widths = c(1, 2, 1))

# Final Stack preserving decentralized local guides
final_controls_with_legend <- (row1) / (row2) + 
  plot_layout(
    guides = "keep",  # This retains unique legends fixed under each individual chart block
    heights = c(1, 1)    
  ) & 
  unified_theme       

# ==============================================================================
# RENDER FINAL COMPOSITE
# ==============================================================================
final_controls_with_legend

ggsave("control_fig2.png", plot = final_controls_with_legend, dpi = 500)




body_22 <- body_22 %>%
  mutate(
    # 1. Create the categorical density groups
    density_group = case_when(
      control == "treatment"      ~ "Mixed",
      str_detect(replicate, "10") ~ "n = 10",
      str_detect(replicate, "20") ~ "n = 20",
      TRUE                        ~ "Unknown"
    ),
    # 2. Convert to a factor and fix the ordering (Mixed -> n = 10 -> n = 20)
    density_group = fct_relevel(density_group, "Mixed", "n = 10", "n = 20")
  )


anova_vicina <- aov(thorax_mm ~ density_group * sex, data = filter(body_22, species == "vicina"))
summary(anova_vicina)

anova_sericata <- aov(thorax_mm ~ density_group * sex, data = filter(body_22, species == "sericata"))
summary(anova_sericata)
TukeyHSD(anova_sericata, "sex")

# 1. Filter the dataset for C. vicina and explicitly convert 'group' to a factor
vicina_eggs_data <- df22_clean %>%
  filter(species == "C. vicina") %>%
  mutate(group = as.factor(group))

# 2. Run the ANOVA using the clean 'group' variable
anova_vicina_eggs <- aov(performance_index ~ x_group, data = vicina_eggs_data)
summary(anova_vicina_eggs)

# 3. View the overall ANOVA summary
summary(anova_vicina_eggs)

sericata_eggs_data <- df22_clean %>%
  filter(species == "L. sericata") %>%
  mutate(group = as.factor(group))

# 2. Run the ANOVA using the clean 'group' variable
anova_sericata_eggs <- aov(performance_index ~ x_group, data = sericata_eggs_data)
summary(anova_sericata_eggs)

# 3. View the overall ANOVA summary
summary(anova_vicina_eggs)

# 1. Filter for C. vicina and ensure density is treated as a factor
vicina_emergence_data <- long_22 %>%
  filter(species == "vicina") %>%
  mutate(density = as.factor(density))

# 2. Run the One-Way ANOVA on the success proportions
anova_vicina_emergence <- aov(days_to_first ~ density, data = vicina_emergence_data)

# 3. View the ANOVA results table
summary(anova_vicina_emergence)

# 1. Filter for C. sericata and ensure density is treated as a factor
sericata_emergence_data <- long_22 %>%
  filter(species == "sericata") %>%
  mutate(density = as.factor(density))

# 2. Run the One-Way ANOVA on the success proportions
anova_sericata_emergence <- aov(days_to_first ~ density, data = sericata_emergence_data)

# 3. View the ANOVA results table
summary(anova_sericata_emergence)
TukeyHSD(anova_sericata_emergence, "density")


