# Load data and packages --------------------------------------------------
source(here::here("R/Library.R"))

#load data
df <- vroom::vroom(here::here("data-raw/csrp3_aav_results.csv"))

#create delta data frame
df_delta <- df %>%
    rename(
        weight_gast = gast_weight,
        weight_ta = ta_weight,
        csrp3_gast = csrp3_gast,
        csrp3_ta = csrp3_ta
    ) %>%
pivot_longer(
    cols = c(weight_gast, weight_ta, csrp3_gast, csrp3_ta, flag_gast, flag_ta),
    names_to = c(".value", "muscle"),
    names_pattern = "(weight|csrp3|flag)_(gast|ta)"
) %>%
    group_by(bl6, muscle) %>%
    mutate(
        delta_weight = weight[aav == "csrp3"] - weight[aav == "egfp"],
        delta_csrp3 = csrp3[aav == "csrp3"] - csrp3[aav == "egfp"],
        delta_flag = flag[aav == "csrp3"] - flag[aav == "egfp"]
    ) %>%
    ungroup() %>%
    filter(aav == "csrp3") %>%
    mutate(
        normalized_weight = delta_weight / weight
    )



# GASTOCNEMIUS MUSCLE WEIGHT ----------------------------------------------

#Linear mixed model
lmm_gast <- lmer(gast_weight ~ aav + (1 | bl6), data = df, REML = FALSE)
fixed_effects_gast <- anova(lmm_gast)

#calculate means for each aav
means_gast <- emmeans(lmm_gast, ~ aav)
print(means_gast)


# TIBIALIS ANTERIOR MUSCLE WEIGHT -----------------------------------------

#Linear mixed model
lmm_ta <- lmer(ta_weight ~ aav + (1 | bl6), data = df, REML = FALSE)
fixed_effects_ta <- anova(lmm_ta)

#calculate means for each aav
means_ta <- emmeans(lmm_ta, ~ aav)
print(means_ta)

# MUSCLE MASS FIGURE ------------------------------------------------------

#Calculate CSRP3-eGFP and 95%CI
contrast_gast <- contrast(means_gast, method = "pairwise", adjust = "none", infer = TRUE) %>%
    as.data.frame() %>%
    mutate(muscle = "gast")

contrast_ta <- contrast(means_ta, method = "pairwise", adjust = "none", infer = TRUE) %>%
    as.data.frame() %>%
    mutate(muscle = "ta")

#merge data frames
emm_muscle <- rbind(contrast_gast, contrast_ta)

#define p-values for figure
p_muscle <- tibble(
    x = c(1, 2),
    y = c(7.5, 6),
    label = c("p<0.001", "p=0.032")
)

#Figure
muscle_fig <- emm_muscle %>%
    ggplot(aes(x = muscle, y = estimate, fill = muscle)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    #geom_jitter(size=3, width=0, alpha = 0.5, stroke = 0)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values=c("#E61717", "#F6CDCD"))+
    scale_x_discrete(labels = c("gast" = "Gastrocnemius", "ta" = "Tibialis anterior")) +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.5, "lines"),
        legend.box.margin = ggplot2::margin(-10, 0, -10, 0),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = element_text(size = 8),
        text = element_text(size = 8),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(-1, 8)) +
    geom_text(data = p_muscle, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.75) +
    xlab("none") +
    ylab(expression(Delta~"Muscle mass (mg)")) +
    ggtitle("Change in muscle mass \n (CSRP3 - sham)")

#Normalize change in muscle mass to total muscle mass

#Get model-estimated EGFP means for normalization
egfp_gast <- means_gast %>%
    as.data.frame() %>%
    filter(aav == "egfp") %>%
    pull(emmean)

egfp_ta <- means_ta %>%
    as.data.frame() %>%
    filter(aav == "egfp") %>%
    pull(emmean)

#Normalize model-based contrast and CI by EGFP model mean
contrast_gast_norm <- contrast_gast %>%
    mutate(
        estimate = estimate / egfp_gast * 100,
        lower.CL = lower.CL / egfp_gast * 100,
        upper.CL = upper.CL / egfp_gast * 100
    )

contrast_ta_norm <- contrast_ta %>%
    mutate(
        estimate = estimate / egfp_ta * 100,
        lower.CL = lower.CL / egfp_ta * 100,
        upper.CL = upper.CL / egfp_ta * 100
    )

#merge data frames
emm_muscle_norm <- rbind(contrast_gast_norm, contrast_ta_norm)

#define p-values for figure
p_muscle <- tibble(
    x = c(1, 2),
    y = c(6, 13),
    label = c("p<0.001", "p=0.032")
)

#Figure
muscle_fig <- emm_muscle_norm %>%
    ggplot(aes(x = muscle, y = estimate, fill = muscle)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("#E61717", "#F6CDCD")) +
    scale_x_discrete(labels = c("gast" = "Gastrocnemius", "ta" = "Tibialis anterior")) +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.5, "lines"),
        legend.box.margin = ggplot2::margin(-10, 0, -10, 0),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = element_text(size = 8),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 6),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    #scale_y_continuous(limits = c(-0.02, 0.14)) +
    geom_text(data = p_muscle, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.75) +
    xlab("none") +
    ylab(expression(Delta~"Muscle mass (%)")) +
    ggtitle("Change in muscle mass \n (CSRP3 - sham)")

ggsave(plot = muscle_fig, here::here('figures/figure_5/aav_muscle.pdf'), height = 65, width = 47, units = "mm")

# Overexpression CSRP3 ----------------------------------------------------


##GASTROC##

#Linear mixed model of overexpression in gastroc
lmm_csrp3_gast <- lmer(csrp3_gast ~ aav + (1 | bl6), data = df, REML = FALSE)
effects_csrp3_gast <- anova(lmm_csrp3_gast)

#calculate means for each aav
means_csrp3_gast <- emmeans(lmm_csrp3_gast, ~ aav)
print(means_csrp3_gast)

##TIBIALIS ANTERIOR##

#Linear mixed model of overexpression in tibialis anterior
lmm_csrp3_ta <- lmer(csrp3_ta ~ aav + (1 | bl6), data = df, REML = FALSE)
effects_csrp3_ta <- anova(lmm_csrp3_ta)

#calculate means for each aav
means_csrp3_ta <- emmeans(lmm_csrp3_ta, ~ aav)
print(means_csrp3_ta)



# OVEREXPRESSION FIGURE ---------------------------------------------------

#Calculate CSRP3-eGFP and 95%CI
emm_gast <- contrast(means_csrp3_gast, method = "pairwise", adjust = "none", infer = TRUE) %>%
    as.data.frame() %>%
    mutate(muscle = "gast")

emm_ta <- contrast(means_csrp3_ta, method = "pairwise", adjust = "none", infer = TRUE) %>%
    as.data.frame() %>%
    mutate(muscle = "ta")

#merge data frames
emm_csrp3 <- rbind(emm_gast, emm_ta)

#define p-values for figure
p_csrp3 <- tibble(
    x = c(1, 2),
    y = c(6000000, 16000000),
    label = c("p=0.010", "p=0.014")
)

#Figure
csrp3_fig <- emm_csrp3 %>%
    ggplot(aes(x = muscle, y = estimate, fill = muscle)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    #geom_jitter(size=3, width=0, alpha = 0.5, stroke = 0)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values=c("#E61717", "#F6CDCD"))+
    scale_x_discrete(labels = c("gast" = "Gastrocnemius", "ta" = "Tibialis anterior")) +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.5, "lines"),
        legend.box.margin = ggplot2::margin(-10, 0, -10, 0),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = element_text(size = 8),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 6),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    #scale_y_continuous(limits = c(-2000000, 17000000)) +
    geom_text(data = p_csrp3, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.75) +
    xlab("none") +
    ylab(expression(Delta~"CSRP3 (a.u.)")) +
    ggtitle("Change in CSRP3 \n (CSRP3 - sham)")

ggsave(plot = csrp3_fig, here::here('figures/figure_5/aav_csrp3.pdf'), height = 65, width = 53, units = "mm")



# FLAG --------------------------------------------------------------------

##GASTROC##

#Linear mixed model of overexpression in gastroc
lmm_flag_gast <- lmer(flag_gast ~ aav + (1 | bl6), data = df, REML = FALSE)
effects_flag_gast <- anova(lmm_flag_gast)

#calculate means for each aav
means_flag_gast <- emmeans(lmm_flag_gast, ~ aav)
print(means_flag_gast)

##TIBIALIS ANTERIOR##

#Linear mixed model of overexpression in tibialis anterior
lmm_flag_ta <- lmer(flag_ta ~ aav + (1 | bl6), data = df, REML = FALSE)
effects_flag_ta <- anova(lmm_flag_ta)

#calculate means for each aav
means_flag_ta <- emmeans(lmm_flag_ta, ~ aav)
print(means_flag_ta)


# FLAG FIGURE -------------------------------------------------------------

#Calculate CSRP3-eGFP and 95%CI
emm_flag_gast <- contrast(means_flag_gast, method = "pairwise", adjust = "none", infer = TRUE) %>%
    as.data.frame() %>%
    mutate(muscle = "gast")

emm_flag_ta <- contrast(means_flag_ta, method = "pairwise", adjust = "none", infer = TRUE) %>%
    as.data.frame() %>%
    mutate(muscle = "ta")

#merge data frames
emm_flag <- rbind(emm_flag_gast, emm_flag_ta)

#define p-values for figure
p_flag <- tibble(
    x = c(1, 2),
    y = c(4500000, 22000000),
    label = c("p<0.001", "p=0.005")
)

#Figure
flag_fig <- emm_flag %>%
    ggplot(aes(x = muscle, y = estimate, fill = muscle)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    #geom_jitter(size=3, width=0, alpha = 0.5, stroke = 0)+
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values=c("#E61717", "#F6CDCD"))+
    scale_x_discrete(labels = c("gast" = "Gastrocnemius", "ta" = "Tibialis anterior")) +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.5, "lines"),
        legend.box.margin = ggplot2::margin(-10, 0, -10, 0),
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = element_text(size = 8),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 6),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    #scale_y_continuous(limits = c(-3000000, 23000000)) +
    geom_text(data = p_flag, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.75) +
    xlab("none") +
    ylab(expression(Delta~"rAAV:CSRP3 FLAG (a.u.)")) +
    ggtitle("Change in rAAV:CSRP3 FLAG \n (CSRP3 - sham)")

ggsave(plot = flag_fig, here::here('figures/figure_5/aav_flag.pdf'), height = 65, width = 53, units = "mm")


# REPLICATION GASTROCNEMIUS -----------------------------------------------

#load data
df_replicate <- vroom::vroom(here::here("data-raw/csrp3_aav_replication.csv"))

#create delta data frame
df_delta_replicate <- df_replicate %>%
    group_by(bl6) %>%
    mutate(
        delta_gast = gast[aav == "csrp3"] - gast[aav == "egfp"]
    ) %>%
    ungroup() %>%
    filter(aav == "csrp3")

#Linear mixed model
lmm_replicate <- lmer(gast ~ aav + (1 | bl6), data = df_replicate, REML = FALSE)
fixed_effects_replicate <- anova(lmm_replicate)

#define p-value for figure
p_replicate <- tibble(
    x = 4,
    y = 25,
    label = "p=0.017"
)

##FIGURE##
replicate_fig <- df_delta_replicate %>%
    ggplot(aes(x = week, y = delta_gast, color = aav)) +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.8, fatten = 1.5, color = "black", alpha = 0.7) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
                 width = 0.4, linewidth = 0.25, color = "black") +
    geom_quasirandom(aes(color = aav),
                     size = 2.5, alpha = 0.7, width = 0.25, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_color_manual(values = "#E61717") +
    scale_x_discrete(labels = "Gastrocnemius") +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    geom_text(data = p_replicate,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 2.75) +
    xlab("Gastrocnemius") +
    ylab(expression(Delta~"Muscle mass (mg)")) +
    ggtitle("Replication \n (CSRP3 - sham)")

ggsave(plot = replicate_fig, here::here('figures/figure_5/replicate_gast_aav.pdf'), height = 65, width = 33, units = "mm")

#bracket for replication figure
bracket_replication <- tibble(
    x = c(1, 1, 2),
    xend = c(2, 1, 2),
    y = c(169, 169, 169),
    yend = c(169, 160, 167)
)

#define p-value for figure
p_replicate <- tibble(
    x = 1.5,
    y = 171,
    label = "p=0.017"
)

##ABSOLUTE VALUES FIGURE##
replicate_fig <- df_replicate %>%
    mutate(aav = factor(aav, levels = c("egfp", "csrp3"))) %>%
    ggplot(aes(x = aav, y = gast, color = aav)) +
    geom_jitter(size = 2.5, alpha = 0.6, width = 0, stroke = 0) +
    scale_color_manual(
        name = NULL,
        values = c("csrp3" = "#E61717", "egfp" = "gray"))+
    geom_line(aes(group = bl6),
              color = "black",
              linewidth = 0.25,
              alpha = 0.7) +
    scale_x_discrete(labels = c("csrp3" = "CSRP3", "egfp" = "Sham")) +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(size = 8),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
    geom_segment(data = bracket_replication, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_replicate, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.5) +
    ylab("Gastrocnemius muscle mass (mg)") +
    xlab(NULL) +
    ggtitle("Replication")

ggsave(plot = replicate_fig, here::here('figures/figure_5/replicate_gast_aav.pdf'), height = 65, width = 33, units = "mm")

