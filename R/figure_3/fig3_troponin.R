source(here::here("R/Library.R"))

#load differential expression between sex with annotations
results_tnn <- vroom::vroom(here::here("data/results_sex_keywords.csv")) %>%
    filter(protein %in% c("TNNC2", "TNNI2", "TNNT3"))

#load data long
data_long_pre <- vroom::vroom(here::here("data/data_log2_long_normalized.csv")) %>%
    filter(trial == "pre")


#create dataframe of mean expression
data_mean_pre <- data_long_pre %>%
    group_by(gene, group, fibertype, sex) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
    ungroup()

#filter for troponin isoforms in data_long
tnn_long <- data_long_pre %>%
    select(!(1)) %>%
    filter(gene %in% c("TNNT3", "TNNC2", "TNNI2"))

# Make sex a factor with male as reference
tnn_long$sex <- factor(tnn_long$sex, levels = c("male", "female"))

#filter for each isoform
tnnc2_df <- tnn_long %>%
    filter(gene == "TNNC2")

tnni2_df <- tnn_long %>%
    filter(gene == "TNNI2")

tnnt3_df <- tnn_long %>%
    filter(gene == "TNNT3")

##TNNC2##
#Linear mixed model
lmm_tnnc2 <- lmer(expression ~ fibertype * sex + (1 | subject), data = tnnc2_df, REML = FALSE)
fixed_effects_tnnc2 <- anova(lmm_tnnc2)

#calculate means for each fibertype for each sex
means_tnnc2 <- emmeans(lmm_tnnc2, ~ sex | fibertype)
print(means_tnnc2)

#run linear model of sex for each fibertype
lm_tnnc2_sex <- contrast(means_tnnc2, method = "pairwise", by = "fibertype")
print(lm_tnnc2_sex)

#run linear model of fibertype for each sex
lm_tnnc2_fibertype <- contrast(means_tnnc2, method = "pairwise", by = "sex")
print(lm_tnnc2_fibertype)

##TNNI2##
#Linear mixed model
lmm_tnni2 <- lmer(expression ~ fibertype * sex + (1 | subject), data = tnni2_df, REML = FALSE)
fixed_effects_tnni2 <- anova(lmm_tnni2)

#calculate means for each fibertype for each sex
means_tnni2 <- emmeans(lmm_tnni2, ~ sex | fibertype)
print(means_tnni2)

#run linear model of sex for each fibertype
lm_tnni2_sex <- contrast(means_tnni2, method = "pairwise", by = "fibertype")
print(lm_tnni2_sex)

#run linear model of fibertype for each sex
lm_tnni2_fibertype <- contrast(means_tnni2, method = "pairwise", by = "sex")
print(lm_tnni2_fibertype)

##TNNT3##
#Linear mixed model
lmm_tnnt3 <- lmer(expression ~ fibertype * sex + (1 | subject), data = tnnt3_df, REML = FALSE)
fixed_effects_tnnt3 <- anova(lmm_tnnt3)

#calculate means for each fibertype for each sex
means_tnnt3 <- emmeans(lmm_tnnt3, ~ sex | fibertype)
print(means_tnnt3)

#run linear model of sex for each fibertype
lm_tnnt3_sex <- contrast(means_tnnt3, method = "pairwise", by = "fibertype")
print(lm_tnnt3_sex)

#run linear model of fibertype for each sex
lm_tnnt3_fibertype <- contrast(means_tnnt3, method = "pairwise", by = "sex")
print(lm_tnnt3_fibertype)

#define brackets for figure
brackets_tnn <- tibble(
    fibertype = c("mhc1", "mhc1", "mhc1",
              "mhc1", "mhc1", "mhc1",
              "mhc2", "mhc2", "mhc2"),
    x = c(0.75, 0.75, 1.25,
          1.75, 1.75, 2.25,
          2.75, 2.75, 3.25),
    xend = c(1.25, 0.75, 1.25,
             2.25, 1.75, 2.25,
             3.25, 2.75, 3.25),
    y = c(16, 16, 16,
          17, 17, 17,
          23, 23, 23),
    yend = c(16, 14.5, 15.5,
             17, 15, 16.5,
             23, 22.5, 22)
)

#define asterix for figure
stars_tnn <- tibble(
    fibertype = c("mhc1", "mhc1", "mhc2"),
    x = c(1, 2, 3),
    y = c(16.5, 17.5, 23.5),
    label = c("*", "*", "*")
)

#boxplot of troponin isoforms
tnn_plot <- tnn_long %>%
    ggplot(aes(x = gene, y = expression, fill = sex)) +
    geom_boxplot(width = 1, linewidth = 0.25, alpha = 0.5, outlier.size = 0, outlier.stroke = 0) +
    geom_jitter(position = position_dodge(width = 1),
                aes(color = sex), size = 1, alpha = 0.5, stroke = 0) +
    scale_fill_manual(values = c("#000000", "#FF7518")) +
    scale_color_manual(values = c("#000000", "#FF7518")) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    scale_y_continuous(limits = c(11, 27)) +
    geom_segment(data = brackets_tnn, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = stars_tnn, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 6) +
    xlab("Protein") +
    ylab("Abundance (a.u.)") +
    ggtitle("Fast troponin isoforms") +
    facet_wrap(
        ~fibertype,
        labeller = as_labeller(c(
            mhc1 = "Type I",
            mhc2 = "Type II"
        )))

ggsave(plot = tnn_plot, here::here('figures/figure_3/tnn_figure.pdf'), height = 55, width = 70, units = "mm")


#bar plot of troponin isoforms
emm_tnnc2 <- emmeans(lm_tnnc2_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    mutate(
        regulated = case_when(
            fibertype == "mhc1" ~ "no",
            fibertype == "mhc2" ~ "no")) %>%
    mutate(protein = "TNNC2")

emm_tnni2 <- emmeans(lm_tnni2_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    mutate(
        regulated = case_when(
            fibertype == "mhc1" ~ "female",
            fibertype == "mhc2" ~ "no")) %>%
    mutate(protein = "TNNI2")

emm_tnnt3 <- emmeans(lm_tnnt3_sex, ~ fibertype | fibertype) %>%
    as.data.frame() %>%
    mutate(
        regulated = case_when(
            fibertype == "mhc1" ~ "female",
            fibertype == "mhc2" ~ "male")) %>%
    mutate(protein = "TNNT3")

emm_tnn <- rbind(emm_tnnc2, emm_tnni2, emm_tnnt3)


##FIGURE##
tnn_plot <- emm_tnn %>%
    ggplot(aes(x = protein, y = estimate, fill = regulated)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values = c("female" = "#000000",
                                 "male" = "#FF7518",
                                 "no" = "gray"),
                      labels = c("male" = "Greater abundance in males", "female" = "Greater abundance in females"),
                      breaks = c("female", "male"),
                      name = NULL) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(2, "mm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    #coord_cartesian(ylim = c(12, 25)) +
    xlab("Protein") +
    ylab("Log2fold difference (male - female)") +
    ggtitle("Fast troponin isoforms") +
    facet_wrap(~fibertype,
               labeller = as_labeller(c(mhc1 = "Type I", mhc2 = "Type II")))




