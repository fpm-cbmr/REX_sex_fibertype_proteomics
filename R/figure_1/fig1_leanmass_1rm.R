source(here::here("R/Library.R"))
#load 1RM and leg lean mass data
lbm_1rm_df <- vroom::vroom(here::here("data-raw/leanmass_1rm_data.csv"))

#Make pre and post to factors in data_long so they appear in correct order
lbm_1rm_df$trial <- factor(lbm_1rm_df$trial, levels = c("pre", "post"))

#create delta dataframe
pre_df <- lbm_1rm_df %>%
    filter(trial == "pre")

post_df <- lbm_1rm_df %>%
    filter(trial == "post")

delta_df <- post_df %>%
    select(!lp_1rm) %>%
    select(!lbm_legs) %>%
    mutate(delta_1rm = post_df$lp_1rm - pre_df$lp_1rm) %>%
    mutate(delta_lbm = post_df$lbm_legs - pre_df$lbm_legs)

#run linear mixed model for 1RM for main effects and interaction
lm_1rm <- lmer(lp_1rm ~ trial * sex + (1 | id), data = lbm_1rm_df, REML = FALSE)
effects_1rm <- anova(lm_1rm)

#calculate means for each trial for each sex
means_1rm <- emmeans(lm_1rm, ~ trial | sex)
print(means_1rm)

#run linear model of trial for each sex
lm_1rm_sex <- contrast(means_1rm, method = "pairwise", by = "sex", adjust = "none")
lm_1rm_sex_df <- summary(lm_1rm_sex)

#run linear mixed model for lean mass for main effects and interaction
lm_leanmass <- lmer(lbm_legs ~ trial * sex + (1 | id), data = lbm_1rm_df, REML = FALSE)
effects_leanmass <- anova(lm_leanmass)

#calculate means for each trial for each sex
means_leanmass <- emmeans(lm_leanmass, ~ trial | sex)
print(means_leanmass)

#run linear model of trial for each sex
lm_leanmass_sex <- contrast(means_leanmass, method = "pairwise", by = "sex", adjust = "none")
lm_leanmass_sex_df <- summary(lm_leanmass_sex)

#bracket for leg lean mass figure
bracket_leanmass <- tibble(
    x = c(1, 1, 2),
    xend = c(2, 1, 2),
    y = c(2.25, 2.25, 2.25),
    yend = c(2.25, 1.75, 2.1)
)

#define p-value for figure
p_leanmass <- tibble(
    x = c(1, 2, 1.5),
    y = c(1.60, 1.9, 2.4),
    label = c("p<0.001", "p<0.001", "p=0.813")
)

#Box Plot of leg lean mass
leg_lean <- delta_df %>%
    ggplot(aes(x = sex, y = delta_lbm, fill = sex)) +
    geom_violin(trim = TRUE, width=1, linewidth = 0.5, alpha=0.5)+
    geom_boxplot(width=0.25, color="black", fill="white", alpha=0.5, outlier.size = 2.5, outlier.stroke = 0)+
    geom_jitter(size=2.5, width=0, aes(color = sex), alpha = 0.5, stroke=0)+
    scale_fill_manual(values=c("#000000", "#FF7518"))+
    scale_color_manual(values = c("#000000", "#FF7518")) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_x_discrete(labels=c("Females", "Males"))+
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )+
    scale_y_continuous(limits = c(-0.5, 2.5)) +
    geom_segment(data = bracket_leanmass, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_leanmass, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Change in leg lean mass (kg)")

##figure of leg lean mass##
leg_lean <- delta_df %>%
    ggplot(aes(x = sex, y = delta_lbm, fill = sex)) +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.8, fatten = 1.5, color = "black", alpha = 0.7) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
                 width = 0.4, linewidth = 0.25) +
    geom_quasirandom(aes(color = sex),
                     size = 2, alpha = 0.7, width = 0.25, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values = c("#000000", "#FF7518")) +
    scale_color_manual(values = c("#000000", "#FF7518")) +
    scale_x_discrete(labels = c("Females", "Males")) +
    scale_y_continuous(limits = c(-0.5, 2.5)) +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    geom_segment(data = bracket_leanmass,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_leanmass,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Change in leg lean mass (kg)")

ggsave(plot = leg_lean, here::here('figures/figure_1/leg_leanmass.pdf'), height = 35, width = 35, units = "mm")

##1RM##

#bracket for 1RM figure
bracket_1rm <- tibble(
    x = c(1, 1, 2),
    xend = c(2, 1, 2),
    y = c(200, 200, 200),
    yend = c(200, 175, 190)
)

#define p-value for figure
p_1rm <- tibble(
    x = c(1, 2, 1.5),
    y = c(160, 175, 215),
    label = c("p<0.001", "p<0.001", "p=0.203")
)


#Box Plot of 1RM leg press
rm <- delta_df %>%
    ggplot(aes(x = sex, y = delta_1rm, fill = sex)) +
    geom_violin(trim = TRUE, width=1, linewidth = 0.5, alpha=0.5)+
    geom_boxplot(width=0.25, color="black", fill="white", alpha=0.5, outlier.size = 2.5, outlier.stroke = 0)+
    geom_jitter(size=2.5, width=0, aes(color = sex), alpha = 0.5, stroke=0)+
    scale_fill_manual(values=c("#000000", "#FF7518"))+
    scale_color_manual(values = c("#000000", "#FF7518")) +
    geom_hline(yintercept=0, linetype="dashed") +
    scale_x_discrete(labels=c("Females", "Males"))+
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )+
    scale_y_continuous(limits = c(-25, 225)) +
    geom_segment(data = bracket_1rm, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_1rm, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Change in 1RM leg press (kg)")


##1RM leg press figure##
rm <- delta_df %>%
    ggplot(aes(x = sex, y = delta_1rm, fill = sex)) +
    stat_summary(fun = mean, geom = "crossbar",
                 width = 0.8, fatten = 1.5, color = "black", alpha = 0.7) +
    stat_summary(fun.data = mean_cl_normal, geom = "errorbar",
                 width = 0.4, linewidth = 0.25) +
    geom_quasirandom(aes(color = sex),
                     size = 2, alpha = 0.7, width = 0.25, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values = c("#000000", "#FF7518")) +
    scale_color_manual(values = c("#000000", "#FF7518")) +
    scale_x_discrete(labels = c("Females", "Males")) +
    scale_y_continuous(limits = c(-25, 225)) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    geom_segment(data = bracket_1rm,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_1rm,
              aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Change in 1RM leg press (kg)")


ggsave(plot = rm, here::here('figures/figure_1/legpress.pdf'), height = 35, width = 35, units = "mm")
