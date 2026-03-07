source(here::here("R/Library.R"))

#load dataframes of protein log2fc with annotations
results_log2fc <- vroom::vroom(here::here("data/results_log2fc_keywords.csv"))

log2fc_df <- readRDS(here::here("data/data_log2fc_long_keywords.rds"))

#load long form data frame
df_long <- readRDS(here::here("data/data_long_keywords.rds"))

#create dataframe of mean log2fc of each protein
df_mean_log2fc <- log2fc_df %>%
    group_by(protein, group, fibertype, sex) %>%
    summarise(mean_log2fc = mean(log2fc, na.rm = TRUE)) %>%
    ungroup()

#create dataframe of mean abundance of each protein
df_mean <- df_long %>%
    group_by(protein, group, fibertype, sex, trial) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
    ungroup()


# MITOCHONDRIAL PROTEINS --------------------------------------------------

#Load mitocarta
mitocarta <- read_excel(here::here('data-raw/mitocarta.xls'))%>%
    dplyr::select('symbol', 'pathways') %>%
    dplyr::rename(protein=symbol)

#Add mitocarta to df
df_mean <- df_mean %>%
    merge(mitocarta, by="protein", all.x = T) %>%
    dplyr::rename(mito = pathways)

oxphos_df <- df_mean %>%
    dplyr::filter(grepl('OXPHOS', mito))

mitoribo_df <- df_mean %>%
    dplyr::filter(grepl('Mitochondrial ribosome', mito))

mito_translation_df <- df_mean %>%
    dplyr::filter(grepl('Translation', mito))

mito_df <- rbind(oxphos_df, mitoribo_df)

#Add mitocarta to log2fc
df_mean_log2fc <- df_mean_log2fc %>%
    merge(mitocarta, by="protein", all.x = T) %>%
    dplyr::rename(mito = pathways)

oxphos_log2fc <- df_mean_log2fc %>%
    dplyr::filter(grepl('OXPHOS', mito)) %>%
    dplyr::mutate(component = "OXPHOS")

mitoribo_log2fc <- df_mean_log2fc %>%
    dplyr::filter(grepl('Mitochondrial ribosome', mito)) %>%
    dplyr::mutate(component = "mitoribosome")

mito_translation_log2fc <- df_mean_log2fc %>%
    dplyr::filter(grepl('Translation', mito))

mito_log2fc <- rbind(oxphos_log2fc, mitoribo_log2fc)


##OXPHOS##
#linear mixed model
lmm_oxphos <- lmer(mean_expression ~ trial * fibertype * sex + (1 | protein), data = oxphos_df, REML = FALSE)
effects_oxphos <- anova(lmm_oxphos)

#calculate means of each trial for each fiber type and sex
means_oxphos <- emmeans(lmm_oxphos, ~ trial | fibertype | sex)
print(means_oxphos)

#run linear mixed model of trial for each fiber type and sex
lm_oxphos_trial <- contrast(means_oxphos, method = "pairwise", by = c("fibertype", "sex"), adjust = "none")
print(lm_oxphos_trial)

#linear mixed model of log2fc
lmm_oxphos_log2fc <- lmer(mean_log2fc ~ fibertype * sex + (1 | protein), data = oxphos_log2fc, REML = FALSE)
effects_oxphos_log2fc <- anova(lmm_oxphos_log2fc)

#calculate means of each fiber type and sex
means_oxphos_log2fc <- emmeans(lmm_oxphos_log2fc, ~ fibertype | sex)
print(means_oxphos_log2fc)

#run linear mixed model of sex for each fiber type
lm_oxphos_log2fc_sex <- contrast(means_oxphos_log2fc, method = "pairwise", by = "fibertype", adjust = "none")
print(lm_oxphos_log2fc_sex)

#run linear mixed model of fiber type for each sex
lm_oxphos_log2fc_fibertype <- contrast(means_oxphos_log2fc, method = "pairwise", by = "sex", adjust = "none")
print(lm_oxphos_log2fc_fibertype)

#define brackets for box plot
brackets_oxphos <- tibble(
    sex = c("female", "female", "female",
            "male", "male", "male"),
    x = c(1, 1, 2,
          1, 1, 2),
    xend = c(2, 1, 2,
             2, 1, 2),
    y = c(1.1, 1.1, 1.1,
          0.9, 0.9, 0.9),
    yend = c(1.1, 0.8, 1.0,
             0.9, 0.5, 0.8)
)

#define p-values for box plot
p_oxphos <- tibble(
    sex = c("female", "female", "female",
                  "male", "male", "male"),
    x = c(1, 2, 1.5,
          1, 2, 1.5),
    y = c(0.7, 0.9, 1.2,
          0.4, 0.7, 1.0),
    label = c("p=0.009", "p=0.032", "p=0.108",
              "p=0.976", "p<0.001", "p<0.001")
)

##OXPHOS BOX PLOT##
oxphos_fig <- oxphos_log2fc %>%
    ggplot(aes(x = fibertype, y = mean_log2fc,
               fill = fibertype,
               group = interaction(fibertype, sex))) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.5, alpha = 0.5,
                position = position_dodge(width = 1),
                aes(color = sex)) +
    geom_boxplot(width=0.25, position = position_dodge(width = 1), color="black", fill="white", alpha=0.5, outlier.size = 1.5, outlier.stroke = 0)+
    geom_jitter(position = position_dodge(width = 1),color = "black", size = 1.5, alpha = 0.25, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("female" = "black", "male" = "#FF7518"),
                       labels = c("female" = "Female", "male" = "Male")) +
    scale_fill_manual(values = c("mhc1" = "#440154FF", "mhc2" = "#67CC5CFF"),
                      labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    scale_x_discrete(labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    guides(fill = guide_legend(title = NULL)) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    facet_wrap(~sex, labeller = labeller(sex = c("female" = "Females", "male" = "Males"))) +
    scale_y_continuous(limits = c(-0.9, 1.2)) +
    geom_segment(data = brackets_oxphos, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_oxphos, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("OXPHOS")

##OXPHOS BAR PLOT##

#Compute emmeans from linear mixed model
emm_oxphos <- emmeans(lmm_oxphos_log2fc, ~ fibertype | sex) %>%
    as.data.frame()

#define brackets for bar plot
brackets_oxphos <- tibble(
    sex = c("female", "female", "female",
            "male", "male", "male"),
    x = c(1, 1, 2,
          1, 1, 2),
    xend = c(2, 1, 2,
             2, 1, 2),
    y = c(0.15, 0.15, 0.15,
          0.15, 0.15, 0.15),
    yend = c(0.15, 0.12, 0.09,
             0.15, 0.05, 0.13)
)

#define p-values for bar plot
p_oxphos <- tibble(
    sex = c("female", "female", "female",
            "male", "male", "male"),
    x = c(1, 2, 1.5,
          1, 2, 1.5),
    y = c(0.11, 0.08, 0.16,
          0.04, 0.12, 0.16),
    label = c("p=0.009", "p=0.032", "p=0.108",
              "p=0.976", "p<0.001", "p<0.001")
)

##OXPHOS BAR PLOT##
oxphos_fig <- ggplot(emm_oxphos, aes(x = fibertype, y = emmean, fill = fibertype)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values = c("mhc1" = "#440154FF", "mhc2" = "#67CC5CFF"),
                      labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    scale_x_discrete(labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    facet_wrap(~sex, labeller = labeller(sex = c("female" = "Females", "male" = "Males"))) +
    scale_y_continuous(limits = c(-0.05, 0.2)) +
    geom_segment(data = brackets_oxphos, aes(x = x, xend = xend, y = y, yend = yend),
                 size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_oxphos, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 2) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    xlab("") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("OXPHOS")

#ggsave(plot = oxphos_fig, here::here('figures/figure_4/oxphos_log2fc.pdf'), height = 70, width = 60, units = "mm")

##MITORIBOSOME##
#linear mixed model
lmm_mitoribo <- lmer(mean_expression ~ trial * fibertype * sex + (1 | protein), data = mitoribo_df, REML = FALSE)
effects_mitoribo <- anova(lmm_mitoribo)

#calculate means of each trial for each fiber type and sex
means_mitoribo <- emmeans(lmm_mitoribo, ~ trial | fibertype | sex)
print(means_mitoribo)

#run linear mixed model of trial for each fiber type and sex
lm_mitoribo_trial <- contrast(means_mitoribo, method = "pairwise", by = c("fibertype", "sex"), adjust = "none")
print(lm_mitoribo_trial)

#linear mixed model of log2fc
lmm_mitoribo_log2fc <- lmer(mean_log2fc ~ fibertype * sex + (1 | protein), data = mitoribo_log2fc, REML = FALSE)
effects_mitoribo_log2fc <- anova(lmm_mitoribo_log2fc)

#calculate means of each fiber type and sex
means_mitoribo_log2fc <- emmeans(lmm_mitoribo_log2fc, ~ fibertype | sex)
print(means_mitoribo_log2fc)

#run linear mixed model of sex for each fiber type
lm_mitoribo_log2fc_sex <- contrast(means_mitoribo_log2fc, method = "pairwise", by = "fibertype", adjust = "none")
print(lm_mitoribo_log2fc_sex)

#run linear mixed model of fiber type for each sex
lm_mitoribo_log2fc_fibertype <- contrast(means_mitoribo_log2fc, method = "pairwise", by = "sex", adjust = "none")
print(lm_mitoribo_log2fc_fibertype)

#define brackets for box plot
brackets_mitoribo <- tibble(
    sex = c("female", "female", "female",
            "male", "male", "male"),
    x = c(1, 1, 2,
          1, 1, 2),
    xend = c(2, 1, 2,
             2, 1, 2),
    y = c(1.1, 1.1, 1.1,
          1.1, 1.1, 1.1),
    yend = c(1.1, 1.0, 1.0,
             1.1, 0.8, 1.0)
)

#define p-values for box plot
p_mitoribo <- tibble(
    sex = c("female", "female", "female",
            "male", "male", "male"),
    x = c(1, 2, 1.5,
          1, 2, 1.5),
    y = c(0.9, 0.9, 1.2,
          0.7, 0.9, 1.2),
    label = c("p=0.007", "p=0.043", "p=0.861",
              "p=0.021", "p<0.001", "p=0.236")
)

##MITORIBOSOME BOX PLOT##
mitoribo_fig <- mitoribo_log2fc %>%
    ggplot(aes(x = fibertype, y = mean_log2fc,
               fill = fibertype,
               group = interaction(fibertype, sex))) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.5, alpha = 0.5,
                position = position_dodge(width = 1),
                aes(color = sex)) +
    geom_boxplot(width=0.25, position = position_dodge(width = 1), color="black", fill="white", alpha=0.5, outlier.size = 1.5, outlier.stroke = 0)+
    geom_jitter(position = position_dodge(width = 1),color = "black", size = 1.5, alpha = 0.25, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("female" = "black", "male" = "#FF7518"),
                       labels = c("female" = "Female", "male" = "Male")) +
    scale_fill_manual(values = c("mhc1" = "#440154FF", "mhc2" = "#67CC5CFF"),
                      labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    scale_x_discrete(labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    guides(fill = guide_legend(title = NULL)) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    facet_wrap(~sex, labeller = labeller(sex = c("female" = "Females", "male" = "Males"))) +
    scale_y_continuous(limits = c(-0.8, 1.2)) +
    geom_segment(data = brackets_mitoribo, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_mitoribo, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("Mitochondrial ribosome")


##MITORIBOSOME BAR PLOT##

#Compute emmeans from linear mixed model
emm_mitoribo <- emmeans(lmm_mitoribo_log2fc, ~ fibertype | sex) %>%
    as.data.frame()

#define brackets for bar plot
brackets_mitoribo <- tibble(
    sex = c("female", "female", "female",
            "male", "male", "male"),
    x = c(1, 1, 2,
          1, 1, 2),
    xend = c(2, 1, 2,
             2, 1, 2),
    y = c(0.21, 0.21, 0.21,
          0.26, 0.26, 0.26),
    yend = c(0.21, 0.18, 0.19,
             0.26, 0.19, 0.24)
)

#define p-values for bar plot
p_mitoribo <- tibble(
    sex = c("female", "female", "female",
            "male", "male", "male"),
    x = c(1, 2, 1.5,
          1, 2, 1.5),
    y = c(0.17, 0.18, 0.22,
          0.18, 0.23, 0.27),
    label = c("p=0.007", "p=0.043", "p=0.861",
              "p=0.021", "p<0.001", "p=0.236")
)

#Mitoribosome figure
mitoribo_fig <- ggplot(emm_mitoribo, aes(x = fibertype, y = emmean, fill = fibertype)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values = c("mhc1" = "#440154FF", "mhc2" = "#67CC5CFF"),
                      labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    scale_x_discrete(labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    facet_wrap(~sex, labeller = labeller(sex = c("female" = "Females", "male" = "Males"))) +
    scale_y_continuous(limits = c(-0.05, 0.30)) +
    geom_segment(data = brackets_mitoribo, aes(x = x, xend = xend, y = y, yend = yend),
                 size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_mitoribo, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 2) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.title.y = element_text(size = 7),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    xlab("") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("Mitochondrial ribosome")

ggsave(plot = mitoribo_fig, here::here('figures/figure_4/mitoribosome_log2fc.pdf'), height = 70, width = 60, units = "mm")

# MITOCHONDRIAL TRANSLATION -----------------------------------------------

#linear mixed model
lmm_mito_translation <- lmer(mean_expression ~ trial * fibertype * sex + (1 | protein), data = mito_translation_df, REML = FALSE)
effects_mito_translation <- anova(lmm_mito_translation)

#calculate means of each trial for each fiber type and sex
means_mito_translation <- emmeans(lmm_mito_translation, ~ trial | fibertype | sex)
print(means_mito_translation)

#run linear mixed model of trial for each fiber type and sex
lm_mito_translation_trial <- contrast(means_mito_translation, method = "pairwise", by = c("fibertype", "sex"), adjust = "none")
print(lm_mito_translation_trial)

#linear mixed model of log2fc
lmm_mito_translation_log2fc <- lmer(mean_log2fc ~ fibertype * sex + (1 | protein), data = mito_translation_log2fc, REML = FALSE)
effects_mito_translation_log2fc <- anova(lmm_mito_translation_log2fc)

#calculate means of each fiber type and sex
means_mito_translation_log2fc <- emmeans(lmm_mito_translation_log2fc, ~ fibertype | sex)
print(means_mito_translation_log2fc)

#run linear mixed model of sex for each fiber type
lm_mito_translation_log2fc_sex <- contrast(means_mito_translation_log2fc, method = "pairwise", by = "fibertype", adjust = "none")
print(lm_mito_translation_log2fc_sex)

#run linear mixed model of fiber type for each sex
lm_mito_translation_log2fc_fibertype <- contrast(means_mito_translation_log2fc, method = "pairwise", by = "sex", adjust = "none")
print(lm_mito_translation_log2fc_fibertype)

##MITOCHONDRIAL TRANSLATION FIGURE##
mitotranslation_fig <- mito_translation_log2fc %>%
    ggplot(aes(x = group, y = mean_log2fc,
               fill = fibertype,
               group = interaction(fibertype, sex))) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.5, alpha = 0.5,
                position = position_dodge(width = 1),
                aes(color = sex)) +
    geom_boxplot(width=0.25, position = position_dodge(width = 1), color="black", fill="white", alpha=0.5, outlier.size = 1.5, outlier.stroke = 0)+
    geom_jitter(position = position_dodge(width = 1),color = "black", size = 1.5, alpha = 0.25, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("female" = "black", "male" = "#FF7518"),
                       labels = c("female" = "Female", "male" = "Male")) +
    scale_fill_manual(values = c("mhc1" = "#440154FF", "mhc2" = "#67CC5CFF"),
                      labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    scale_x_discrete(labels = c("female_mhc1" = "Female \n type I", "female_mhc2" = "Female \n type II",
                                "male_mhc1" = "Male \n type I", "male_mhc2" = "Male \n type II")) +
    guides(fill = guide_legend(title = NULL)) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5)
    ) +
    #scale_y_continuous(limits = c(-1, 1)) +
    #geom_text(data = p_mito_complex, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("Mitochondrial translation")


# MITOCHONDRIAL ENCODED PROTEINS ------------------------------------------
mt_encoded <- mito_df %>%
    filter(protein %in% c("MT-ND1", "MT-ND2", "MT-ND4", "MT-ND5", "MT-CYB", "MT-CO1", "MT-CO2", "MT_CO3", "MT-ATP6", "MT-ATP8"))

mt_encoded_log2fc <- mito_log2fc %>%
    filter(protein %in% c("MT-ND1", "MT-ND2", "MT-ND4", "MT-ND5", "MT-CYB", "MT-CO1", "MT-CO2", "MT_CO3", "MT-ATP6", "MT-ATP8"))

results_mt_encoded <- results_log2fc %>%
    filter(protein %in% c("MT-ND1", "MT-ND2", "MT-ND4", "MT-ND5", "MT-CYB", "MT-CO1", "MT-CO2", "MT_CO3", "MT-ATP6", "MT-ATP8"))

#linear mixed model
lmm_mt_encoded <- lmer(mean_expression ~ trial * fibertype * sex + (1 | protein), data = mt_encoded, REML = FALSE)
effects_mt_encoded <- anova(lmm_mt_encoded)

#No main effect of trial (p=0.386) or interactions

#calculate means of each trial for each fiber type and sex
means_mt_encoded <- emmeans(lmm_mt_encoded, ~ trial | fibertype | sex)
print(means_mt_encoded)

#run linear mixed model of trial for each fiber type and sex
lm_mt_encoded_trial <- contrast(means_mt_encoded, method = "pairwise", by = c("fibertype", "sex"), adjust = "none")
print(lm_mt_encoded_trial)


##MITOCHONDRIAL ENCODED PROTEINS FIGURE##
mito_encoded_fig <- mt_encoded_log2fc %>%
    ggplot(aes(x = group, y = mean_log2fc,
               fill = fibertype,
               group = interaction(fibertype, sex))) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.5, alpha = 0.5,
                position = position_dodge(width = 1),
                aes(color = sex)) +
    geom_boxplot(width=0.25, position = position_dodge(width = 1), color="black", fill="white", alpha=0.5, outlier.size = 1.5, outlier.stroke = 0)+
    geom_jitter(position = position_dodge(width = 1),color = "black", size = 1.5, alpha = 0.25, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("female" = "black", "male" = "#FF7518"),
                       labels = c("female" = "Female", "male" = "Male")) +
    scale_fill_manual(values = c("mhc1" = "#440154FF", "mhc2" = "#67CC5CFF"),
                      labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
    scale_x_discrete(labels = c("female_mhc1" = "Female \n type I", "female_mhc2" = "Female \n type II",
                                "male_mhc1" = "Male \n type I", "male_mhc2" = "Male \n type II")) +
    guides(fill = guide_legend(title = NULL)) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        legend.key.size = unit(0.3, "cm"),
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 6, face = "bold", hjust = 0.5)
    ) +
    #scale_y_continuous(limits = c(-1, 1)) +
    #geom_text(data = p_mito_complex, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("Mitochondrial encoded proteins")

