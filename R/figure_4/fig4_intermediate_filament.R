source(here::here("R/Library.R"))

#load dataframes of protein log2fc with annotations
results_log2fc <- vroom::vroom(here::here("data/results_log2fc_keywords.csv"))

log2fc_df <- readRDS(here::here("data/data_log2fc_long_keywords.rds"))

#load long form data frame
df_long <- readRDS(here::here("data/data_long_keywords.rds"))

#create dataframe of mean log2fc of each protein
df_mean_log2fc <- log2fc_df %>%
    group_by(protein, group, fibertype, sex, gobp) %>%
    summarise(mean_log2fc = mean(log2fc, na.rm = TRUE)) %>%
    ungroup()

#create dataframe of mean abundance of each protein
df_mean <- df_long %>%
    group_by(protein, group, fibertype, sex, trial, gobp) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
    ungroup()

# INTERMEDIATE FILAMENT PROTEINS ------------------------------------------

#filter data frames for intermediate filament proteins
if_log2fc_df <-  df_mean_log2fc %>%
    filter(
        grepl("intermediate filament organization|intermediate filament cytoskeleton organizaton|intermediate filament-based process",
              gobp, ignore.case = TRUE) |
            protein %in% c("KRT16", "KRT6A", "KRT6C", "KRT71", "KRT6B", "KRT14", "KRT17", "KRT5", "DSP", "KRT9", "KRT1", "KRT2", "KRT10", "VIM",
                           "NES", "SYNC", "AGFG1", "DES"))

if_df <- df_mean %>%
    filter(
        grepl("intermediate filament organization|intermediate filament cytoskeleton organizaton|intermediate filament-based process",
              gobp, ignore.case = TRUE) |
            protein %in% c("KRT16", "KRT6A", "KRT6C", "KRT71", "KRT6B", "KRT14", "KRT17", "KRT5", "DSP", "KRT9", "KRT1", "KRT2", "KRT10", "VIM",
                           "NES", "SYNC", "AGFG1", "DES"))



#linear mixed model of intermediate filament proteins
lmm_if <- lmer(mean_expression ~ trial * fibertype * sex + (1 | protein), data = if_df, REML = FALSE)
effects_if <- anova(lmm_if)

#calculate means of each trial for each fiber type and sex
means_if <- emmeans(lmm_if, ~ trial | fibertype | sex)
print(means_if)

#run linear mixed model of trial for each fiber type and sex
lm_if_trial <- contrast(means_if, method = "pairwise", by = c("fibertype", "sex"), adjust = "none")
print(lm_if_trial)

#linear mixed model of intermediate filament log2fc
lmm_if_fc <- lmer(mean_log2fc ~ fibertype * sex + (1 | protein), data = if_log2fc_df, REML = FALSE)
effects_if_fc <- anova(lmm_if_fc)

#calculate means of each fiber type and sex
means_if_fc <- emmeans(lmm_if_fc, ~ fibertype | sex)
print(means_if_fc)

#run linear mixed model of sex for each fiber type
lm_if_fc_sex <- contrast(means_if_fc, method = "pairwise", by = "fibertype", adjust = "none")
print(lm_if_fc_sex)

#run linear mixed model of fibertype for each sex
lm_if_fc_trial <- contrast(means_if_fc, method = "pairwise", by = "sex", adjust = "none")
print(lm_if_fc_trial)

#define brackets for figure
brackets_if <- tibble(
    fibertype = c("mhc1", "mhc1", "mhc1",
                  "mhc2", "mhc2", "mhc2"),
    x = c(1, 1, 2,
          1, 1, 2),
    xend = c(2, 1, 2,
             2, 1, 2),
    y = c(2.2, 2.2, 2.2,
          2.5, 2.5, 2.5),
    yend = c(2.2, 2, 1.1,
            2.5, 2.3, 1.4)
)

#define p-values for figure
p_if <- tibble(
    fibertype = c("mhc1", "mhc1", "mhc1",
                  "mhc2", "mhc2", "mhc2"),
    x = c(1, 2, 1.5,
          1, 2, 1.5),
    y = c(1.9, 0.9, 2.4,
          2.2, 1.3, 2.7),
    label = c("p<0.001", "p=0.655", "p<0.001",
              "p<0.001", "p=0.015", "p<0.001")
)

#Box plot of mean intermediate filament log2fc
if_fig <- if_log2fc_df %>%
    ggplot(aes(x = sex, y = mean_log2fc, fill = sex)) +
    geom_violin(trim = TRUE, width=1, linewidth = 0.5, alpha=0.5)+
    geom_boxplot(width=0.25, color="black", fill="white", alpha=0.5, outlier.size = 2, outlier.stroke = 0)+
    geom_jitter(size=2, width=0, alpha = 0.5, stroke=0)+
    scale_fill_manual(values=c("#000000", "#FF7518"))+
    scale_color_manual(values = c("#000000", "#FF7518")) +
    scale_x_discrete(labels=c("Females", "Males"))+
    geom_hline(yintercept=0, linetype="dashed") +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = element_text(size = 7),
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )+
    facet_wrap(~fibertype, labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    scale_y_continuous(limits = c(-1, 3)) +
    geom_segment(data = brackets_if, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_if, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("Intermediate filament proteins ")

#Bar plot of mean intermediate filament log2fc
#Compute emmeans from linear mixed model
emm_if <- emmeans(lmm_if_fc, ~ fibertype | sex) %>%
    as.data.frame()

#define brackets for figure
brackets_if <- tibble(
    fibertype = c("mhc1", "mhc1", "mhc1",
                  "mhc2", "mhc2", "mhc2"),
    x = c(1, 1, 2,
          1, 1, 2),
    xend = c(2, 1, 2,
             2, 1, 2),
    y = c(1.7, 1.7, 1.7,
          1.7, 1.7, 1.7),
    yend = c(1.7, 1.3, 0.6,
             1.7, 1.5, 0.9)
)

#define p-values for figure
p_if <- tibble(
    fibertype = c("mhc1", "mhc1", "mhc1",
                  "mhc2", "mhc2", "mhc2"),
    x = c(1, 2, 1.5,
          1, 2, 1.5),
    y = c(1.2, 0.5, 1.8,
          1.4, 0.8, 1.8),
    label = c("p<0.001", "p=0.655", "p<0.001",
              "p<0.001", "p=0.015", "p<0.001")
)


##FIGURE##
if_fig <- ggplot(emm_if, aes(x = sex, y = emmean, fill = sex)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    scale_fill_manual(values=c("female" = "#000000", "male" = "#FF7518"))+
    scale_x_discrete(labels=c("Females", "Males"))+
    geom_hline(yintercept=0, linetype="dashed", linewidth = 0.25) +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill=NA, linewidth = 0.5),
        panel.grid.minor=element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = element_text(size = 7),
        text = element_text(size = 6),
        axis.text.x= element_text(color="black", size = 6),
        axis.text.y= element_text(color="black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    )+
    facet_wrap(~fibertype, labeller = labeller(fibertype = c("mhc1" = "Type I", "mhc2" = "Type II"))) +
    scale_y_continuous(limits = c(-0.5, 2)) +
    geom_segment(data = brackets_if, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_if, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("none") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("Intermediate filament proteins ")

ggsave(plot = if_fig, here::here('figures/figure_4/i_filament_log2fc.pdf'), height = 70, width = 60, units = "mm")




