source(here::here("R/Library.R"))

#load data of fibers
df_ihc <- read_excel(here::here('data-raw/fibers_ihc_data.xlsx'))

# CREATE FIBER TYPE DATA FRAME ---------------------------------------------------
df_ihc <- df_ihc %>%
    pivot_longer(
        cols = matches("^(fibers|ihc_fibers|csa|dist|dist_ihc|area_dist)_mhc"),
        names_to = c(".value", "mhc"),
        names_pattern = "(.*)_mhc(.*)"
    ) %>%
    mutate(
        mhc = paste0("mhc", mhc))

#write.csv(df_ihc,
#here::here("data-raw/fibers_ihc_data.csv"))

# Make trial a factor
df_ihc$trial <- factor(df_ihc$trial, levels = c("pre", "post"))

#filter for fiber types
df_mhc1 <- df_ihc %>%
    filter(mhc == "mhc1")

df_mhc2 <- df_ihc %>%
    filter(mhc == "mhc2")


# SINGLE FIBER DISTRIBUTION -----------------------------------------------

##TYPE I FIBERS##

#run linear mixed model for main effects and interaction
lm_dist_mhc1 <- lmer(dist ~ sex * trial + (1 | id), data = df_mhc1)
dist_mhc1_effects <- anova(lm_dist_mhc1)

#Main effect of sex (P=0.0509)

# Calculate means for each trial within each sex
means_dist_mhc1 <- emmeans(lm_dist_mhc1, ~ trial | sex)
means_dist_mhc1_df <- summary(means_dist_mhc1)

#run linear mixed model of trial for each sex
lm_dist_mhc1_trial <- contrast(means_dist_mhc1, method = "pairwise", by = "sex", adjust = "none")
print(lm_dist_mhc1_trial)

#No change in fiber type in females (p=0.3867) or males (p=0.7871)

#run linear mixed model of sex for each trial
lm_dist_mhc1_sex <- contrast(means_dist_mhc1, method = "pairwise", by = "trial", adjust = "none")
print(lm_dist_mhc1_sex)

#greater proportion of type I fibers in females at pre (p=0.0349)

##TYPE II FIBERS##

#run linear mixed model for main effects and interaction
lm_dist_mhc2 <- lmer(dist ~ sex * trial + (1 | id), data = df_mhc2)
dist_mhc2_effects <- anova(lm_dist_mhc2)

#Main effect of sex (P=0.0509)

# Calculate means for each trial within each sex
means_dist_mhc2 <- emmeans(lm_dist_mhc2, ~ trial | sex)
means_dist_mhc2_df <- summary(means_dist_mhc2)

#run linear mixed model of trial for each sex
lm_dist_mhc2_trial <- contrast(means_dist_mhc2, method = "pairwise", by = "sex", adjust = "none")
print(lm_dist_mhc2_trial)

#No change in fiber type in females (p=0.3867) or males (p=0.7871)

#run linear mixed model of sex for each trial
lm_dist_mhc2_sex <- contrast(means_dist_mhc2, method = "pairwise", by = "trial", adjust = "none")
print(lm_dist_mhc2_sex)

#greater proportion of type II fibers in males at pre (p=0.0349)


# AREA PROPORTION ---------------------------------------------------------

##TYPE I FIBERS##

#run linear mixed model for main effects and interaction
lm_area_mhc1 <- lmer(area_dist ~ sex * trial + (1 | id), data = df_mhc1)
area_mhc1_effects <- anova(lm_area_mhc1)

#No main effect sex (P=0.1523)

# Calculate means for each trial within each sex
means_area_mhc1 <- emmeans(lm_area_mhc1, ~ trial | sex)
means_area_mhc1_df <- summary(means_area_mhc1)

#run linear mixed model of trial for each sex
lm_area_mhc1_trial <- contrast(means_area_mhc1, method = "pairwise", by = "sex", adjust = "none")
print(lm_area_mhc1_trial)

#No change in fiber type area distribution in females (p=0.4667) or males (p=0.7799)

#run linear mixed model of sex for each trial
lm_area_mhc1_sex <- contrast(means_area_mhc1, method = "pairwise", by = "trial", adjust = "none")
print(lm_area_mhc1_sex)

#No sex differences in area proportion of type I fibers at pre (p=0.1115) or post (p=2856)

##TYPE II FIBERS##

#run linear mixed model for main effects and interaction
lm_area_mhc2 <- lmer(area_dist ~ sex * trial + (1 | id), data = df_mhc2)
area_mhc2_effects <- anova(lm_area_mhc2)

#No main effect sex (P=0.1523)

# Calculate means for each trial within each sex
means_area_mhc2 <- emmeans(lm_area_mhc2, ~ trial | sex)
means_area_mhc2_df <- summary(means_area_mhc2)

#run linear mixed model of trial for each sex
lm_area_mhc2_trial <- contrast(means_area_mhc2, method = "pairwise", by = "sex", adjust = "none")
print(lm_area_mhc2_trial)

#No change in fiber type area distribution in females (p=0.4667) or males (p=0.7799)

#run linear mixed model of sex for each trial
lm_area_mhc2_sex <- contrast(means_area_mhc2, method = "pairwise", by = "trial", adjust = "none")
print(lm_area_mhc2_sex)

#No sex differences in area proportion of type I fibers at pre (p=0.1115) or post (p=0.2856)



# CREATE SUMMARY DATA FRAME -----------------------------------------------


##NUMBER OF FIBERS##
means_fibers_df <- df_ihc %>%
    group_by(sex, mhc) %>%
    summarise(
        mean = mean(fibers, na.rm = TRUE),
        sd   = sd(fibers, na.rm = TRUE),
        n    = n(),
        .groups = "drop"
    )

##FIBERTYPE DISTRIBUTION##
means_dist_df <- df_ihc %>%
    group_by(sex, mhc) %>%
    summarise(
        mean = mean(dist, na.rm = TRUE),
        sd   = sd(dist, na.rm = TRUE),
        n    = n(),
        .groups = "drop"
    )

##AREA PROPORTION##
means_area_df <- df_ihc %>%
    group_by(sex, mhc) %>%
    summarise(
        mean = mean(area_dist, na.rm = TRUE),
        sd   = sd(area_dist, na.rm = TRUE),
        n    = n(),
        .groups = "drop"
    )


# FIGURE OF FIBER TYPE DISTRIBUTION ---------------------------------------

##SINGLE FIBER DISTRIBUTION##

#Define lines for figure
lines_dist <- tibble(
    mhc = c("mhc1", "mhc2"),
    x = c(1, 1),
    xend = c(2, 2),
    y = c(100, 100),
    yend = c(100, 100)
)

#define brackets for figure
brackets_dist <- tibble(
    mhc = c("mhc1", "mhc1", "mhc1",
            "mhc1", "mhc1", "mhc1",
            "mhc2", "mhc2", "mhc2",
            "mhc2", "mhc2", "mhc2"),
    x = c(0.75, 0.75, 1.25,
          1.75, 1.75, 2.25,
          0.75, 0.75, 1.25,
          1.75, 1.75, 2.25),
    xend = c(1.25, 0.75, 1.25,
             2.25, 1.75, 2.25,
             1.25, 0.75, 1.25,
             2.25, 1.75, 2.25),
    y = c(80, 80, 80,
          85, 85, 85,
          90, 90, 90,
          90, 90, 90),
    yend = c(80, 75, 65,
             85, 80, 80,
             90, 85, 85,
             90, 86, 86))

#define asterix for figure
p_dist <- tibble(
    mhc = c("mhc1", "mhc1", "mhc1",
            "mhc2", "mhc2", "mhc2"),
    x = c(1, 2, 1.5,
          1, 2, 1.5),
    y = c(85, 90, 105,
          95, 95, 105),
    label = c("p=0.035", "p=0.184", "Main effect of training: p=0.671",
              "p=0.035", "p=0.184", "Main effect of training: p=0.671"))

#Figure
dist_fig <- df_ihc %>%
    ggplot(aes(x = trial, y = dist*100, fill = sex)) +
    stat_summary(fun = mean, geom = "bar",
                 position = position_dodge(width = 0.95),
                 width = 0.9, color = NA, alpha = 0.8) +
    geom_jitter(size = 2,
                aes(color = sex),
                alpha = 0.5,
                stroke = 0,
                position = position_jitterdodge(jitter.width = 0, dodge.width = 0.95)) +
    scale_fill_manual(values=c("female" = "#000000", "male" = "#FF7518"),
                      labels = c("female" = "Females", "male" = "Males"),
                      name = NULL)+
    scale_color_manual(values = c("female" = "#000000", "male" = "#FF7518"),
                       name = NULL, guide = "none") +
    scale_x_discrete(labels=c("pre" = "Pre", "post" = "Post"))+
    facet_wrap(~mhc, labeller = labeller(mhc = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "right",
        legend.key.size = unit(4, "mm"),
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )+
    coord_cartesian(ylim = c(10, 100)) +
    #geom_segment(data = lines_dist, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    #geom_segment(data = brackets_dist, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    #geom_text(data = p_dist, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Fiber type distribution (%)") +
    ggtitle("Distribution of single muscle fibers")

ggsave(plot = dist_fig, here::here('figures/supplementary_figure_1/fiber_distribution.pdf'), height = 60, width = 140, units = "mm")

##FIBER AREA DISTRIBUTION##

#Define lines for figure
lines_dist_area <- tibble(
    mhc = c("mhc1", "mhc2"),
    x = c(1, 1),
    xend = c(2, 2),
    y = c(102, 102),
    yend = c(102, 102)
)

#define brackets for figure
brackets_dist_area <- tibble(
    mhc = c("mhc1", "mhc1", "mhc1",
            "mhc1", "mhc1", "mhc1",
            "mhc2", "mhc2", "mhc2",
            "mhc2", "mhc2", "mhc2"),
    x = c(0.75, 0.75, 1.25,
          1.75, 1.75, 2.25,
          0.75, 0.75, 1.25,
          1.75, 1.75, 2.25),
    xend = c(1.25, 0.75, 1.25,
             2.25, 1.75, 2.25,
             1.25, 0.75, 1.25,
             2.25, 1.75, 2.25),
    y = c(95, 95, 95,
          90, 90, 90,
          85, 85, 85,
          85, 85, 85),
    yend = c(95, 92, 92,
             90, 85, 85,
             85, 80, 80,
             85, 80, 80))

#define p-values for figure
p_dist_area <- tibble(
    mhc = c("mhc1", "mhc1", "mhc1",
            "mhc2", "mhc2", "mhc2"),
    x = c(1, 2, 1.5,
          1, 2, 1.5),
    y = c(98, 95, 105,
          90, 90, 105),
    label = c("p=0.112", "p=0.286", "Main effect of training: p=0.739",
              "p=0.112", "p=0.286", "Main effect of training: p=0.739"))

#Figure
area_dist_fig <- df_ihc %>%
    ggplot(aes(x = trial, y = area_dist*100, fill = sex)) +
    stat_summary(fun = mean, geom = "bar",
                 position = position_dodge(width = 0.95),
                 width = 0.9, color = NA, alpha = 0.8) +
    geom_jitter(size = 2,
                aes(color = sex),
                alpha = 0.5,
                stroke = 0,
                position = position_jitterdodge(jitter.width = 0, dodge.width = 0.95)) +
    scale_fill_manual(values=c("#000000", "#FF7518"))+
    scale_color_manual(values = c("#000000", "#FF7518")) +
    scale_x_discrete(labels=c("pre" = "Pre", "post" = "Post"))+
    facet_wrap(~mhc, labeller = labeller(mhc = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
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
    coord_cartesian(ylim = c(10, 100)) +
    #geom_segment(data = lines_dist_area, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    #geom_segment(data = brackets_dist_area, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    #geom_text(data = p_dist_area, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Fiber area proportion (%)") +
    ggtitle("Area proportion of muscle fiber cross-sections")

ggsave(plot = area_dist_fig, here::here('figures/supplementary_figure_1/fiber_area_distribution.pdf'), height = 60, width = 120, units = "mm")
