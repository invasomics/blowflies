
#Required packages used across all analysis
library(readr)
library(lme4)
library(MuMIn)
library(glmmTMB)
library(effects)
library(brglm)
library(brglm2)
library(gtools)
library(emmeans)
library(car)
library(dplyr)
library(AER)
library(betareg)
library(ggplot2)
library(mass)
library(multcomp)
library(lmtest)
library(statmod)
library(cowplot)
library(patchwork)
library(ggplot2)



####### Reading in the data for intraspecific data analysis#####
Data <- read.csv("Nunn_et_al_2026_data.csv")
library(dplyr)

Data %>%
  group_by(Density, Treatment) %>%
  summarise(n_unique_groups = n_distinct(Group))

Intra_comp_data <- Data |> filter(Treatment == "Fl" | Treatment == "Fs") |> #filter just for intraspecific competition data
  mutate(Density = as.factor(Density)) |> #ensure density is a factor, not continuous
  mutate(species_sex = paste(Species, Sex, sep = "_")) #create species_sex column for later facets in ggplot

Inter_comp_data <- Data #no filter, retain all data including intraspecific which is needed for the control


# Fit a linear model
model1 <- lm(Thorax.length ~ Density + Species + Sex + Density:Species + Sex:Species, data = Intra_comp_data)
Anova(model1, type=3)
summary(model1)

#### Perform post-hoc analysis####
emmeans_model1 <- emmeans(model1, pairwise ~ Density * Species)
emmeans_model1

####Sex specific means######
emmeans(model1, ~ Species * Sex * Density)
df <- as.data.frame(emmeans(model1, ~ Species * Sex * Density))

#####Calculating percentage differences between treatments####
df %>%
  arrange(Species, Sex, Density) %>%
  group_by(Species, Sex) %>%
  mutate(
    prev_emmean = lag(emmean),
    pct_decrease = (emmean - prev_emmean) / prev_emmean * 100
  )

####Overall species levels emmeans contrasts between densities
emm <- emmeans(model1, ~ Density | Species)

within_species_contrasts <- contrast(
  emm,
  method = "pairwise",
  adjust = "tukey"
)

within_species_contrasts

# EMMs for plotting
pred_lm1 <- emmeans(model1, c("Density", "Species", "Sex")) %>%
  as_tibble() %>%
  mutate(species_sex = paste(Species, Sex, sep = "_"))

########## make a named mapping of the original levels -> plotmath strings#####
label_map <- c(
  "Stygia_f"  = "italic('C. stygia')~'(F)'",
  "Stygia_m"  = "italic('C. stygia')~'(M)'",
  "Lucilia_f" = "italic('L. sericata')~'(F)'",
  "Lucilia_m" = "italic('L. sericata')~'(M)'"
)
##### add a new column with those strings, and set factor levels to control order####
pred_lm1 <- pred_lm1 %>%
  mutate(
    species_sex_lab = label_map[as.character(species_sex)],
    species_sex_lab = factor(species_sex_lab, levels = label_map) # preserves order
  )


#######Generate Figure 2 A (Intraspecific competition and thorax length)  ########
intra_comp_plot_2 <- ggplot(pred_lm1,
                            aes(x = Density,
                                y = emmean,
                                colour = species_sex,
                                fill = species_sex)) +
  scale_fill_manual(values = c("#90EE90", "#90EE90", "#FFD700", "#FFD700")) +
  scale_colour_manual(values = c("#90EE90", "#90EE90", "#FFD700", "#FFD700")) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = emmean - SE,
      ymax = emmean + SE
    ),
    width = 0.15
  ) +
  geom_line(
    aes(group = Sex),
    linetype = "dotted",
    colour = "black"
  ) +
  facet_grid(
    rows = vars(species_sex_lab),
    scales = "free_y",
    labeller = label_parsed
  ) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  theme_bw() +
  xlab("Density") +
  ylab("Body size (mm)") +
  theme(
    legend.position = "none",
    strip.text.y.right = element_text(angle = 90)   # flipped vertical text
  )

intra_comp_plot_2

# -----------------------------
#####Load and create required survival stats#####
# -----------------------------
Surv <- read.csv("Nunn_et_al_2026_survival.csv") %>%
  mutate(
    Treatment = as.factor(Treatment),
    Density   = as.factor(Density),
    Group     = as.factor(Group),
    Date      = as.factor(Date),
    prop_stygia_survived  = ifelse(Total.stygia == 0, 0, Surviving.Styg / Total.stygia),
    prop_lucilia_survived = ifelse(Total.Lucilia == 0, 0, Surviving.luc / Total.Lucilia),
    Lfail = Total.Lucilia - Surviving.luc,
    Sfail = Total.stygia  - Surviving.Styg,
    total_surviving_per_density =
      (Surviving.luc + Surviving.Styg) / as.numeric(as.character(Density))
  )


# -----------------------------
# Intraspecific competition subset
# -----------------------------
Intra_comp_Surv <- Surv %>%
  filter(Treatment %in% c("Fl", "Fs")) %>%
  mutate(
    Sp = case_when(
      Treatment == "Fs" ~ "C.stygia",
      Treatment == "Fl" ~ "Lucilia spp"
    )
  ) %>%
  droplevels()

# Beta regression correction
Intra_comp_Surv$total_surviving_per_density <- 
  (Intra_comp_Surv$total_surviving_per_density * (nrow(Intra_comp_Surv) - 1) + 0.5) / nrow(Intra_comp_Surv)

# -----------------------------
# Beta regression model
# -----------------------------
m_full <- betareg(total_surviving_per_density ~ Density + Sp + Density:Sp,
                  data = Intra_comp_Surv)

m_null <- betareg(total_surviving_per_density ~ 1,
                  data = Intra_comp_Surv)

lrtest(m_null, m_full)
summary(m_full)
Anova(m_full, type = 3)

# -----------------------------
# EMMs + within-species contrasts
# -----------------------------
emm_full <- emmeans(m_full, pairwise ~ Density * Sp)
emm_full

# Clean within-species contrasts
within_species <- contrast(
  emmeans(m_full, ~ Density | Sp),
  method = "pairwise",
  adjust = "tukey"
)

within_species

# -----------------------------
# Predictions for plotting
# -----------------------------
pred_lm2 <- emmeans(m_full, c("Density", "Sp")) %>%
  as_tibble()

# -----------------------------
# #### Plotting intraspecific emergence ( Figure 2B)########
# -----------------------------
intra_comp_plot_emergence <- ggplot(pred_lm2,  aes(x = Density,  y = emmean,  colour = Sp,  fill = Sp)) +
  scale_fill_manual(values = c("#FFD700", "#90EE90")) +
  scale_colour_manual(values = c("#FFD700", "#90EE90")) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.15) +
  geom_line(aes(group = Sp), linetype = "dotted", colour = "black") +
  facet_grid(rows = vars(Sp), scales = "fixed") +
  theme_bw() +
  xlab("Density") +
  ylab("Proportion of surviving larvae") +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(color = "gray90"),
    strip.background = element_blank(),
    strip.text = element_blank()
  )

# -----------------------------
# Combine with body-size plot
# -----------------------------
big_theme <- theme(
  text = element_text(size = 15),
  axis.title = element_text(size = 16),
  axis.text  = element_text(size = 16)
)

p1 <- intra_comp_plot_2 + big_theme
p2 <- intra_comp_plot_emergence + big_theme

plot_grid(
  p1, p2,
  ncol = 2,
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold"
)




###### Interspecific comparisons ##########

#### Creating the isolated species datasets####
Data$Date      <- as.factor(Data$Date)
Data$Species   <- as.factor(Data$Species)
Data$Density   <- as.factor(Data$Density)
Data$Treatment <- as.factor(Data$Treatment)
Data$Group     <- as.factor(Data$Group)
Data$Sex       <- as.factor(Data$Sex)


Lucilia.3 <- Data[Data$Species == 'Lucilia', ]
str(Lucilia.3)
Lucilia.3$Treatment <- droplevels(Lucilia.3$Treatment)
levels(Lucilia.3$Treatment)

Stygia.2.D <- Data[Data$Species == 'Stygia', ]
Stygia.2.D $Treatment <- droplevels(Stygia.2.D $Treatment)
levels(Stygia.2.D$Treatment)

#######Subsetting Lucilia by Sex ########
Luc.male   <- subset(Lucilia.3, Sex == "m")
Luc.female <- subset(Lucilia.3, Sex == "f")


######Fitting and testing the Lucilia male model #######
Luc.male$Treatment <- factor(Luc.male$Treatment,
                             levels = c("Fl", "50/50", "70s", "70l"))
LM_male <- lm(Thorax.length ~ Density + Treatment:Density, data = Luc.male)

emmeans(LM_male, ~ Treatment | Density)
hist(resid(LM_male))
plot(density(resid(LM_male, type="response")))
plot(density(resid(LM_male, type="pearson")))

summary(LM_male)

plot(allEffects(LM_male))
pairs(emmeans(LM_male, ~Treatment:Density))

######Extracting effects########
effects_male <- data.frame(summary(LM_male)$coef)
effects_male$upper <- effects_male$Estimate + 1.96*effects_male$Std..Error
effects_male$lower <- effects_male$Estimate - 1.96*effects_male$Std..Error

######## Keep only interaction terms#######
effects_male_int.l <- effects_male[grep(":", rownames(effects_male)), ]

### Assign the correct density and treatment#######
effects_male_int.l$Comp <- c("50/50","50/50","70s","70s","70l","70l")
effects_male_int.l$Dens <- c("10","50","10","50","10","50")
effects_male_int.l$sig  <- ifelse(effects_male_int.l$Pr...t.. <= 0.05, "sig", "nonsig")
effects_male_int.l


###### plotting lucilia male interspecific effects (Figure C)########
Plot_Male.Luc <- ggplot(effects_male_int.l,
                        aes(x = Estimate, y = Comp, xmin = lower, xmax = upper)) +
  theme_bw() +
  geom_errorbar(height = 0) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "black") +
  geom_point(aes(fill = sig, size = sig), shape = 21) +
  scale_fill_manual(values = c("sig"="black","nonsig"="white")) +
  scale_size_manual(values = c("sig"=2,"nonsig"=2)) +
  facet_grid(. ~ Dens) +
  xlab("Effect size") + ylab("Competition") +  scale_y_discrete(limits = c("50/50","70s","70l")) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = "#90EE90")
  )







######## Lucilia female interspecific effects ###########
Luc.female$Treatment <- factor(Luc.female$Treatment,
                               levels = c("Fl", "50/50", "70s", "70l"))
LM_female <- lm(Thorax.length ~ Density + Treatment:Density, data = Luc.female)

emmeans(LM_female, ~ Treatment | Density)
summary(LM_female)

###### Extracting Lucilia female effects######
effects_female.luc <- data.frame(summary(LM_female)$coef)
effects_female.luc $upper <- effects_female.luc $Estimate + 1.96*effects_female.luc $Std..Error
effects_female.luc $lower <- effects_female.luc $Estimate - 1.96*effects_female.luc $Std..Error

#Extracting interaction effects######
effects_female_int.L <- effects_female.luc [grep(":", rownames(effects_female.luc )), ]

### Assign the correct density and treatment#######
effects_female_int.L $Comp <- c("50/50","50/50","70s","70s","70l","70l")
effects_female_int.L $Dens <- c("10","50","10","50","10","50")
effects_female_int.L $sig  <- ifelse(effects_female_int.L$Pr...t.. <= 0.05, "sig", "nonsig")
effects_female_int.L 


########Plotting Lucilia female Interspecific effects plot (Figure 3D) #######
Plot_Female.luc <- ggplot(effects_female_int.L ,
                          aes(x = Estimate, y = Comp, xmin = lower, xmax = upper)) +
  theme_bw() +
  geom_errorbar(height = 0) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "black") +
  geom_point(aes(fill = sig, size = sig), shape = 21) +
  scale_fill_manual(values = c("sig"="black","nonsig"="white")) +
  scale_size_manual(values = c("sig"=2,"nonsig"=2)) +
  facet_grid(. ~ Dens) +
  xlab("Effect size") + 
  ylab("Competition") +
  scale_y_discrete(limits = c("50/50","70s","70l")) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = "#90EE90")
  )









########## Generating isolated stygia interspecific effects#########
Styg.male   <- subset(Stygia.2.D, Sex == "m")
Styg.female <- subset(Stygia.2.D, Sex == "f")

######### Stygia male isolated model ###########
Styg.male$Treatment <- factor(Styg.male$Treatment,
                              levels = c("Fs", "50/50", "70s", "70l"))

LM_Styg_male <- lm(Thorax.length ~ Density + Treatment:Density, data = Styg.male)
emmeans(LM_Styg_male, ~ Treatment | Density)

summary(LM_Styg_male)
plot(allEffects(LM_Styg_male))
hist(resid(LM_Styg_male))
plot(density(resid(LM_Styg_male, type="response")))
plot(density(resid(LM_Styg_male, type="pearson")))

###################### Emmeans tests#############################
pairs(emmeans(LM_Styg_male, ~Treatment))
pairs(emmeans(LM_Styg_male.t, ~Treatment:Density))

##### Extracting effects ########
effects_male <- data.frame(summary(LM_Styg_male)$coef)
effects_male$upper <- effects_male$Estimate + 1.96 * effects_male$Std..Error
effects_male$lower <- effects_male$Estimate - 1.96 * effects_male$Std..Error

# Keep only interaction terms and assing corrects densities
effects_male_int <- effects_male[grep(":", rownames(effects_male)), ]
effects_male_int$Comp <- c("50/50","50/50","70s","70s","70l","70l")
effects_male_int$Dens <- c("10","50","10","50","10","50")
effects_male_int$sig <- ifelse(effects_male_int$Pr...t.. <= 0.05, "sig", "nonsig")
effects_male_int




########Plotting stygia male interspecific effects plot (Figure 3A) #######
Styg_Male_Plot <- ggplot(effects_male_int,
                         aes(x = Estimate, y = Comp, xmin = lower, xmax = upper)) +
  theme_bw() +
  geom_errorbar(height = 0) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "black") +
  geom_point(aes(fill = sig, size = sig), shape = 21) +
  scale_fill_manual(values = c("sig"="black","nonsig"="white")) +
  scale_size_manual(values = c("sig"=2,"nonsig"=2)) +
  facet_grid(. ~ Dens) +
  xlab("") +
  ylab("Species composition") +
  scale_y_discrete(
    limits = c("50/50","70s","70l"),
    labels = c("50/50" = "50%S 50%L",
               "70s"   = "70% S",
               "70l"   = "70% L")
  ) +
  
  # >>> X‑axis ticks every 0.2 <<<
  scale_x_continuous(breaks = seq(-1, 1, by = 0.2)) +
  
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = "#FFD700")
  )








############ Stygia female interspecific effects ###################
Styg.female$Treatment <- factor(Styg.female$Treatment,
                                levels = c("Fs", "50/50", "70s", "70l"))
LM_Styg_female <- lm(Thorax.length ~ Density + Treatment:Density, data = Styg.female)

emmeans(LM_Styg_female, ~ Treatment | Density)
summary(LM_Styg_female)
plot(allEffects(LM_Styg_female))


hist(resid(LM_Styg_female))
plot(density(resid(LM_Styg_female, type="response")))
plot(density(resid(LM_Styg_female, type="pearson")))

pairs(emmeans(LM_Styg_female, ~Treatment))

effects_female <- data.frame(summary(LM_Styg_female)$coef)
effects_female$upper <- effects_female$Estimate + 1.96 * effects_female$Std..Error
effects_female$lower <- effects_female$Estimate - 1.96 * effects_female$Std..Error

########## Keeping only interaction terms ###########
effects_female_int <- effects_female[grep(":", rownames(effects_female)), ]
effects_female_int$Comp <- c("50/50","50/50","70s","70s","70l","70l")
effects_female_int$Dens <- c("10","50","10","50","10","50")
effects_female_int$sig <- ifelse(effects_female_int$Pr...t.. <= 0.05, "sig", "nonsig")


########Plotting stygia female interspecific effects plot (Figure 3B) #######
Styg_Female_Plot <- ggplot(effects_female_int,
                           aes(x = Estimate, y = Comp, xmin = lower, xmax = upper)) +
  theme_bw() +
  geom_errorbar(height = 0) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "black") +
  geom_point(aes(fill = sig, size = sig), shape = 21) +
  scale_fill_manual(values = c("sig"="black","nonsig"="white")) +
  scale_size_manual(values = c("sig"=2,"nonsig"=2)) +
  facet_grid(. ~ Dens) +
  xlab("Effect size") + ylab("Competition") +
  scale_y_discrete(limits = c("50/50","70s","70l")) +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_rect(fill = "#FFD700")
  )

Styg_Female_Plot


######### Stitching together the four plots#######

# Shared y-axis scale for all four plots
comp_scale <- scale_y_discrete(
  limits = c("50/50","70s","70l"),
  labels = c(
    "50/50" = "50:50:S:L",
    "70s"   = "70% S",
    "70l"   = "70% L"
  )
)
# Shared theme for larger fonts and no axis titles
big_theme <- theme(
  text = element_text(size = 12),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = "none",
  axis.text.x = element_text(size = 10, angle = 60, hjust = 1)
)

# Plot versions with matching themes and y axis scale
Styg_Male_clean    <- Styg_Male_Plot   + comp_scale + big_theme
Styg_Female_clean  <- Styg_Female_Plot + comp_scale + big_theme
Luc_Male_clean     <- Plot_Male.Luc    + comp_scale + big_theme
Luc_Female_clean   <- Plot_Female.luc  + comp_scale + big_theme

Combined_4panel <-
  (Styg_Male_clean + Styg_Female_clean) /
  (Luc_Male_clean  + Luc_Female_clean) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 12)
    )
  )
##### View the combined panel output#######
Combined_4panel


################## Interspecific adult emergence plots ####################
#########################################################################
######## Subsetting for survival based on species ##########
lucilia_Survivaliso <- subset(Surv, select = c(Treatment, Density, Surviving.luc, Lfail, Date, Total.Lucilia, prop_lucilia_survived))
stygia_Survivaliso <- subset(Surv, select = c(Treatment, Density, Surviving.Styg, Sfail , Date, Total.stygia,prop_stygia_survived))


lucilia_Survivaliso$Treatment <- factor(lucilia_Survivaliso$Treatment,levels=c("Fl","50/50","70l","70s"))
lucilia_Survivaliso_clean <- na.omit(lucilia_Survivaliso)

stygia_Survivaliso$Treatment <- factor(stygia_Survivaliso$Treatment,levels=c("Fs","50/50","70l","70s"))
stygia_Survivaliso_clean <- na.omit(stygia_Survivaliso)

Stygia.response <- cbind (stygia_Survivaliso_clean$Surviving.Styg, stygia_Survivaliso_clean$Sfail)
lucilia.response <- cbind (lucilia_Survivaliso_clean$Surviving.luc, lucilia_Survivaliso_clean$Lfail)


#### Creating the glm for lucilia#######
GLMM_Lucilia<- glm(lucilia.response ~ Density+Density:Treatment, family = quasibinomial, data=lucilia_Survivaliso_clean)
plot(density(resid(GLMM_Lucilia , type='response'))) #plot residuals
plot(density(resid(GLMM_Lucilia , type='pearson'))) #plot residuals
summary(GLMM_Lucilia)

#Post hoc tests
pairs(emmeans(GLMM_Lucilia, ~Treatment))


###### Extracting coefficients from the GLM#################### 
effects_total<-data.frame(summary(GLMM_Lucilia)$coef) #Put GLM effecs summary into new data frame
effects_total$upper <- effects_total$Estimate + 1.96*effects_total$Std..Error #Add upper Std error
effects_total$lower <- effects_total$Estimate - 1.96*effects_total$Std..Error #Add lower Std error


##### Filter the data to include only rows that contain interaction terms#######
effects_total_interactions <- effects_total[grep(":", rownames(effects_total)), ]
effects_total_interactions$Comp<-c("50/50","50/50","70l","70l","70s","70s")
effects_total_interactions$Dens <- rep(c("10","50"),3)
effects_total_interactions$significance <- ifelse(effects_total_interactions$Pr...t.. <= 0.05, "sig", "nonsig")
effects_total_interactions 

########Plotting lucilia interspecific adult emergence#######
Luc.emergeplot <- ggplot(data = effects_total_interactions, aes(x = Estimate, y = Comp, xmin = lower, xmax = upper)) +
  theme_bw() +
  geom_errorbar(height = 0) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey69", linewidth = 0.8) +  # Updated to 'linewidth'
  geom_point(aes(fill = significance, size = significance), shape = 21) +  # Corrected to use 'significance' instead of 'sig'
  xlab("Effect size") +
  ylab("") +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    strip.background = element_rect(fill = "#90EE90"), # Light green background
    plot.margin = margin(0.75, 0.75, 0.75, 0.75, "cm")
  ) +
  scale_fill_manual(values = c("sig" = "black", "nonsig" = "white")) +
  scale_y_discrete(name = "Competition", 
                   limits = c("50/50", "70s", "70l")) +  # Updated to reflect the correct competition categories
  scale_size_manual(values = c("sig" = 2, "nonsig" = 2), drop = FALSE) +
  facet_grid(. ~ Dens, labeller = label_wrap_gen(width = 15), scales = "free") +
  guides(fill = "none", size = "none", alpha = "none", shape = guide_legend(override.aes = list(size = 5))) +  # Corrected guides
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4))

Luc.emergeplot


######## Generating interspecific C. stygia adult emergence success ###########
###############################################################################
stygia_Survivaliso$Treatment <- droplevels(stygia_Survivaliso$Treatment)
stygia_Survivaliso$Treatment <- factor(stygia_Survivaliso$Treatment,levels=c("Fs","50/50","70l","70s"))


# dropping and viewing unused levels
droplevels(stygia_Survivaliso_clean$Treatment) 
levels(stygia_Survivaliso_clean$Treatment)

GLMM_StygiaT<- glm(Stygia.response ~ Density+Density:Treatment, family = quasibinomial, data=stygia_Survivaliso_clean)
plot(density(resid(GLMM_StygiaT , type='response'))) #plot residuals
plot(density(resid(GLMM_StygiaT , type='pearson'))) #plot residuals
summary(GLMM_StygiaT)
pairs(emmeans(GLMM_StygiaT , ~Treatment))


######### Extracting coefficients from stygia GLM##############
effects_total.Se<-data.frame(summary(GLMM_StygiaT )$coef) #Put GLM effects summary into new data frame
effects_total.Se$upper <- effects_total.Se$Estimate + 1.96*effects_total.Se$Std..Error #Add upper Std error
effects_total.Se$lower <- effects_total.Se$Estimate - 1.96*effects_total.Se$Std..Error #Add lower Std error


######## Extracting only interaction terms #######################
effects_total_interactions.SE <- effects_total.Se[grep(":", rownames(effects_total.Se)), ]
effects_total_interactions.SE $Comp<-c("50/50","50/50","70l","70l","70s","70s")
effects_total_interactions.SE$Dens<-c("10","50","10","50","10","50")
effects_total_interactions.SE$sig<-ifelse(effects_total_interactions.SE$Pr...t..<=0.05,"sig","nonsig")
effects_total_interactions.SE #looks right

########Plotting stygia interspecific adult emergence (Figure 4) #######
Stygia.EmergPlot <- ggplot(data = effects_total_interactions.SE, aes(x = Estimate, y = Comp, xmin = lower, xmax = upper)) +
  theme_bw() +
  geom_errorbar(height = 0) +
  geom_vline(xintercept = 0, linetype = "longdash", colour = "grey69", size = 0.8) +
  geom_point(aes(fill = sig, size = sig), shape = 21) +
  xlab("Effect size") +
  ylab("") +
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(angle = 0, vjust = 0.5, color = "black"),
    axis.text.x = element_text(size = 8, color = "black"),
    strip.background = element_rect(fill = "#FFD700"),
    plot.margin = margin(0.75, 0.75, 0.75, 0.75, "cm")
  ) +
  scale_fill_manual(values = c("sig" = "black", "nonsig" = "white")) +
  scale_y_discrete(name = "Competition", 
                   limits = c("50/50","70s","70l")) +  # Updated to reflect the correct competition categories
  scale_size_manual(values = c("sig" = 2, "nonsig" = 2), drop = FALSE) +
  facet_grid(. ~ Dens, labeller = label_wrap_gen(width = 15), scales = "free") +
  guides(fill = FALSE, size = FALSE, alpha = FALSE, shape = guide_legend(override.aes = list(size = 5))) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4))



# Shared y-axis scale
comp_scale <- scale_y_discrete(
  limits = c("50/50","70s","70l"),
  labels = c(
    "50/50" = "50:50:S:L",
    "70s"   = "70% S",
    "70l"   = "70% L"
  )
)
# Shared theme with larger fonts
big_theme <- theme(
  text = element_text(size = 15),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = "none"
)

# Apply to the two plots
Stygia_clean <- Stygia.EmergPlot + comp_scale + big_theme
Luc_clean    <- Luc.emergeplot + comp_scale + big_theme

####Combining plots#####
Combined_2panel <- 
  (Stygia_clean + Luc_clean) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 16),
      plot.margin = margin(t = 10, r = 10, b = 40, l = 10),
      axis.title.x = element_blank()
    )
  )
###### View the combined interspecific adult emergence output######
Combined_2panel
