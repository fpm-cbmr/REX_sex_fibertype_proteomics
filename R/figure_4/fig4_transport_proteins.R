source(here::here("R/Library.R"))

#load dataframes of protein log2fc with annotations
results_log2fc <- vroom::vroom(here::here("data/results_log2fc_keywords.csv"))

log2fc_df <- readRDS(here::here("data/data_log2fc_long_keywords.RDS"))

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


# ONE CARBON-COMPOUND AND GAS TRANSPORT -----------------------------------

transport_log2fc_df <-  df_mean_log2fc %>%
    filter(grepl("one-compound carbon transport|gas transport", gobp, ignore.case = TRUE) |
               protein %in% c("AQP1","HBA2","CA4","CA2","SLC4A1","HBB","HBD","HBG1","HBG2","RHAG","MB","BPGM"))

transport_df <- df_mean %>%
    filter(grepl("one-compound carbon transport|gas transport", gobp, ignore.case = TRUE) |
               protein %in% c("AQP1","HBA2","CA4","CA2","SLC4A1","HBB","HBD","HBG1","HBG2","RHAG","MB","BPGM"))


#linear mixed model of transport proteins
lmm_transport <- lmer(mean_expression ~ trial * fibertype * sex + (1 | protein), data = transport_df, REML = FALSE)
effects_transport <- anova(lmm_transport)

#calculate means of each trial for each fiber type and sex
means_transport <- emmeans(lmm_transport, ~ trial | fibertype | sex)
print(means_transport)

#run linear mixed model of trial for each fiber type and sex
lm_transport_trial <- contrast(means_transport, method = "pairwise", by = c("fibertype", "sex"), adjust = "none")
lm_transport_trial_df <- summary(lm_transport_trial)

#linear mixed model of transport proteins log2fc
lmm_transport_fc <- lmer(mean_log2fc ~ fibertype * sex + (1 | protein), data = transport_log2fc_df, REML = FALSE)
effects_transport_fc <- anova(lmm_transport_fc)

#calculate means of each fiber type and sex
means_transport_fc <- emmeans(lmm_transport_fc, ~ fibertype | sex)
print(means_transport_fc)

#run linear mixed model of sex for each fiber type
lm_transport_fc_sex <- contrast(means_transport_fc, method = "pairwise", by = "fibertype", adjust = "none")
lm_transport_fc_sex_df <- summary(lm_transport_fc_sex)

#run linear mixed model of fiber type for each sex
lm_transport_fc_fibertype <- contrast(means_transport_fc, method = "pairwise", by = "sex", adjust = "none")
lm_transport_fc_fibertype_df <- summary(lm_transport_fc_fibertype)

# define lines for figure
lines_transport <- tibble(
    sex   = c("female", "male"),
    x     = c(1, 1),
    xend  = c(2, 2),
    y     = c(0.5, 0.5),
    yend  = c(0.5, 0.5)
)

# define brackets for figure
brackets_transport <- tibble(
    sex   = c("female", "female", "female",
              "male", "male", "male"),
    x     = c(1, 1, 2,
              1, 1, 2),
    xend  = c(2, 1, 2,
              2, 1, 2),
    y     = c(0.3, 0.3, 0.3,
              0.05, 0.05, 0.05),
    yend  = c(0.3, 0.2, -0.1,
              0.05, -0.05, -0.1)
)

# define p-values for figure
p_transport <- tibble(
    sex   = c("female", "female",
              "male", "male"),
    x     = c(1.5, 1.5,
              1.5, 1.5),
    y     = c(0.55, 0.35,
              0.55, 0.10),
    label = c("p=0.001", "p<0.001",
              "p=0.001", "p=0.064")
)

#Box plot
transport_fig <- transport_log2fc_df %>%
    ggplot(aes(x = fibertype, y = mean_log2fc,
               fill = fibertype,
               group = interaction(fibertype, sex))) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.25, alpha = 0.5,
                position = position_dodge(width = 1),
                aes(color = sex)) +
    geom_boxplot(width=0.25, linewidth = 0.25, position = position_dodge(width = 1), color="black", fill="white", alpha=0.5, outlier.size = 0, outlier.stroke = 0)+
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
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
    #scale_y_continuous(limits = c(-0.8, 0.6)) +
    geom_segment(data = brackets_transport, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_segment(data = lines_transport, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_transport, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold change (post - pre)") +
    ggtitle("Gas transport")

ggsave(plot = transport_fig, here::here('figures/figure_4/transport_log2fc.pdf'), height = 60, width = 60, units = "mm")

#Bar plot
#Compute emmeans from linear mixed model
emm_transport <- emmeans(lmm_transport_fc, ~ fibertype | sex) %>%
    as.data.frame()

# define lines for figure
lines_transport <- tibble(
    sex   = c("female", "male"),
    x     = c(1, 1),
    xend  = c(2, 2),
    y     = c(0.3, 0.3),
    yend  = c(0.3, 0.3)
)

# define brackets for figure
brackets_transport <- tibble(
    sex   = c("female", "female", "female",
              "male", "male", "male"),
    x     = c(1, 1, 2,
              1, 1, 2),
    xend  = c(2, 1, 2,
              2, 1, 2),
    y     = c(0.10, 0.10, 0.10,
              0.10, 0.10, 0.10),
    yend  = c(0.10, 0.05, 0.05,
              0.10, 0.05, 0.05)
)

# define p-values for figure
p_transport <- tibble(
    sex   = c("female", "female",
              "male", "male"),
    x     = c(1.5, 1.5,
              1.5, 1.5),
    y     = c(0.35, 0.15,
              0.35, 0.15),
    label = c("p<0.001", "p<0.001",
              "p<0.001", "p=0.064")
)

##FIGURE##
transport_fig <- ggplot(emm_transport, aes(x = fibertype, y = emmean, fill = fibertype)) +
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
    scale_y_continuous(limits = c(-0.60, 0.40)) +
    geom_segment(data = brackets_transport, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_segment(data = lines_transport, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_transport, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
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
    ggtitle("One-carbon compound \n and gas transport")




