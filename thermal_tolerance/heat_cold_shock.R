#cold + heat tolerance
setwd("/nesi/project/uow03920/11_thermal_performance_feb")
df <- read.csv("thermal_tolerance_data.csv")

# Convert to factors
df$sex <- factor(df$sex, levels = c("M", "F"))
df$species <- factor(df$species)
df$cold.heat <- factor(df$cold.heat, levels = c("cold", "heat"))

# Summarise mean and CI
summary_df <- df %>%
  group_by(species, sex, cold.heat) %>%
  summarise(
    mean_coma = mean(coma_time, na.rm = TRUE),
    se = sd(coma_time, na.rm = TRUE)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    lower = mean_coma - 1.96*se,
    upper = mean_coma + 1.96*se
  )

####STATS!!!!!!!
##cold, linear model with log transform
library(lme4)
cold <- subset(df, cold.heat == "cold")
m_cold <- lm(
  log(coma_time) ~ species * sex + thorax_length,
  data = cold
)

#set contrasts for type iii anova tests
options(contrasts = c("contr.sum", "contr.poly"))
#anova for cold
Anova(m_cold, type = "III")

#check residuals, fit, etc.
library(DHARMa)
sim_res <- simulateResiduals(m_cold)
plot(sim_res)
testDispersion(sim_res)
testOutliers(sim_res)

#post hoc emmeans
library(emmeans)
emmeans(m_cold, pairwise~ species * sex)

#back transform cold
em_cold <- emmeans(m_cold, ~ species * sex)
em_cold_resp <- summary(em_cold, type = "response")
em_cold_resp
contrast(em_cold, method = "pairwise", type = "response")


#heat
heat <- subset(df, cold.heat == "heat")
m_heat <- lm(
  log(coma_time) ~ species * sex + thorax_length,
  data = heat
)
summary(m_heat)
emmeans(m_heat, pairwise~ species * sex)

Anova(m_heat, type = "III")

#check
sim_res <- simulateResiduals(m_heat)
plot(sim_res)
testDispersion(sim_res)
testOutliers(sim_res)

#back transform heat
em_heat <- emmeans(m_heat, ~ species * sex)
em_heat_resp <- summary(em_heat, type = "response")
em_heat_resp
contrast(em_heat, method = "pairwise", type = "response")

#put the model results all into one data frame
em_cold_df <- as.data.frame(summary(em_cold, type = "response")) %>%
  mutate(cold.heat = "cold")

em_heat_df <- as.data.frame(summary(em_heat, type = "response")) %>%
  mutate(cold.heat = "heat")

em_df <- bind_rows(em_cold_df, em_heat_df)


#plot without jitter
em_df <- em_df %>%
  mutate(
    facet_label = factor(
      cold.heat,
      levels = c("cold", "heat"),
      labels = c("A", "B")
    )
  )

plot_emm <- ggplot(
  em_df,
  aes(
    x = species,
    y = response /60,
    color = sex
  )
) +
  geom_point(
    size = 5,
    position = position_dodge(width = 0.4)
  ) +
  geom_errorbar(
    aes(
      ymin = lower.CL / 60,
      ymax = upper.CL / 60
    ),
    width = 0.2,
    size = 1,
    position = position_dodge(width = 0.4)
  ) +
  facet_wrap(~facet_label, scales = "free") +
  scale_color_manual(
    values = c("M" = "#16A660", "F" = "#A6165C"),
    labels = c("M" = "Male", "F" = "Female")
  ) +
  scale_x_discrete(
    labels = c(
      Calliphora_vicina = "C. vicina",
      Lucilia_sericata = "L. sericata"
    )
  ) +
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(face = "italic", color = "black", size = 18),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_text(face = "bold", color = "black"),
    legend.text = element_text(color = "black", size = 16),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    x = NULL,
    y = "Coma time (min)",
    color = NULL
  )
plot_emm
    
ggsave(
  filename = file.path("cold_heat_tolerance.png"),
  plot = plot_emm, dpi = 600
)

plot_emm <- ggplot(
  em_df,
  aes(
    x = species,
    y = response / 60,
    colour = cold.heat,
    shape = sex
  )
) +
  geom_point(
    size = 5,
    position = position_dodge(width = 0.4),
    stroke = 1.2   # makes open circles clearer
  ) +
  geom_errorbar(
    aes(
      ymin = lower.CL / 60,
      ymax = upper.CL / 60
    ),
    width = 0.2,
    size = 1,
    position = position_dodge(width = 0.4)
  ) +
  facet_wrap(~facet_label, scales = "free") +
  
  # Cold vs warm colours
  scale_colour_manual(
    values = c(
      "cold" = "#6380F2",   # cool blue
      "heat" = "#D11B40"    # warm red
    )
  ) +
  
  # Filled vs open circles for sex
  scale_shape_manual(
    values = c(
      "M" = 16,  # filled circle
      "F" = 1    # open circle
    ),
    labels = c("M" = "Male", "F" = "Female")
  ) +
  
  scale_x_discrete(
    labels = c(
      Calliphora_vicina = "C. vicina",
      Lucilia_sericata = "L. sericata"
    )
  ) +
  
  theme_minimal(base_size = 18) +
  theme(
    axis.text.x = element_text(face = "italic", color = "black", size = 18),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 16),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    x = NULL,
    y = "Coma time (min)",
    colour = "Treatment",
    shape = NULL
  ) + guides(colour = "none")

plot_emm


ggsave(
  filename = file.path("cold_heat_tolerance.png"),
  plot = plot_emm, dpi = 600
)
plot_emm <- ggplot(
  em_df,
  aes(
    x = species,
    y = response / 60,
    colour = species,
    shape = sex
  )
) +
  geom_point(
    size = 5,
    position = position_dodge(width = 0.4),
    stroke = 1.2
  ) +
  geom_errorbar(
    aes(
      ymin = lower.CL / 60,
      ymax = upper.CL / 60
    ),
    width = 0.2,
    size = 1,
    position = position_dodge(width = 0.4)
  ) +
  facet_wrap(~facet_label, scales = "free") +
  
  # Species colours (NO LEGEND)
  scale_colour_manual(
    values = c(
      Calliphora_vicina = "blue",
      Lucilia_sericata = "#C87F17"
    ),
    guide = "none"
  ) +
  
  # Sex shapes (THIS is the only legend)
  scale_shape_manual(
    values = c(
      "M" = 16,
      "F" = 17
    ),
    labels = c(
      "M" = "Male",
      "F" = "Female"
    )
  ) +
  
  scale_x_discrete(
    labels = c(
      Calliphora_vicina = "C. vicina",
      Lucilia_sericata = "L. sericata"
    )
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(face = "italic", color = "black", size = 16),
    axis.text.y = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(color = "black", size = 14),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", color = "black", size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  labs(
    x = NULL,
    y = "Coma time (min)",
    colour = NULL,
    shape = NULL
  )

plot_emm

ggsave(
  filename = file.path("cold_heat_tolerance.png"),
  plot = plot_emm,
  dpi = 600
)
