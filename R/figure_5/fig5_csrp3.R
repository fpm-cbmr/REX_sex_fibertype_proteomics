source(here::here("R/Library.R"))


# TRIAL x FIBERTYPE INTERACTION -------------------------------------------


#load log2fold change data
data_log2fc <- vroom::vroom(here::here("data/data_log2fc.csv")) %>%
    rename_with(~ "protein", colnames(.)[1]) %>%
    column_to_rownames("protein")

#load metadata, create groups and filter to match log2fc data
metadata_log2fc <- vroom::vroom(here::here("data-raw/Metadata_REX_pooled_fibers.csv")) %>%
    filter(!subject == "rex09") %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE) %>%
    dplyr::filter(trial == "post") %>%
    dplyr::select(!fibers) %>%
    dplyr::select(!CSA) %>%
    column_to_rownames("sample_id")

#filter for proteins only regulated in type II fibers
mhc2_log2fc <- data_log2fc %>%
    rownames_to_column("protein") %>%
    filter(protein %in% c("ALDH1B1", "S100A6", "UCHL1", "CASP3", "CSRP3", "CES2", "PSME2", "CAMK2D", "GSN", "GBE1", "COMTD1", "ATP1A1",
                          "S100A4", "SYNC", "CBR1", "CALU", "TUBB2B", "ANXA2", "ME2", "METAP1", "NAMPT", "NPEPL1", "CAPN2", "ACY1", "GLUD1",
                          "MSN", "HMOX2", "GPD2", "CAPZA2", "CA2", "CACNB1", "SOD3", "TINAGL1", "MYLK2")) %>%
    column_to_rownames("protein")


#create summarized experiment
se_log2fc_mhc2 <- SummarizedExperiment(
    assays = list(counts = as.matrix(mhc2_log2fc)),
    colData = DataFrame(metadata_log2fc)
)

#create design matrix for interaction between log2fold change in type I and type II fibers
design_log2fc <- model.matrix(~0 + se_log2fc_mhc2$fibertype)
colnames(design_log2fc) <- c("typeI", "typeII")

#define comparison using contrast
contrast_log2fc <- makeContrasts(typeII - typeI, levels = design_log2fc)

#correlation between samples from same subject
correlation_log2fc <- duplicateCorrelation(assay(se_log2fc_mhc2), design_log2fc, block = se_log2fc_mhc2$subject)

#Use Bayes statistics to assess trial x fibertype interaction between type I and type II fibers
fit_log2fc <- eBayes(lmFit(assay(se_log2fc_mhc2), design_log2fc, block = se_log2fc_mhc2$subject, correlation = correlation_log2fc$consensus.correlation))
ebayes_log2fc <- eBayes(contrasts.fit(fit_log2fc, contrast_log2fc))

#extract results
results_mhc2_regulated <- topTable(ebayes_log2fc, coef = 1, number = Inf, sort.by = "logFC")

#filter for proteins that are more regulated in type II than type I fibers
regulated_mhc2 <- results_mhc2_regulated %>%
    filter(adj.P.Val < 0.05)



# CSRP3 -------------------------------------------------------------------


#load data long form with keywords
df_long <- readRDS(here::here("data/data_long_keywords.rds")) %>%
    select(!2) #random numbers

#filter for csrp3
csrp3_df <- df_long %>%
    filter(protein == "CSRP3")

#linear mixed model
lmm_csrp3 <- lmer(expression ~ trial * fibertype * sex + (1 | subject), data = csrp3_df, REML = FALSE)
fixed_effects_csrp3 <- anova(lmm_csrp3)

#calculate means of each trial for each fiber type
means_csrp3 <- emmeans(lmm_csrp3, ~ trial | fibertype)
print(means_csrp3)

#linear mixed model of trial for each fiber type
lm_csrp3_trial <- contrast(means_csrp3, method = "pairwise", by = "fibertype", adjust = "none")
print(lm_csrp3_trial)

#Load log2fc long form with keywords
df_log2fc <- readRDS(here::here("data/data_log2fc_long_keywords.rds"))

#filter for CSRP3, which  is a positive regulator of myogenesis and plays a role in sensing mechanical stretch
csrp3 <- df_log2fc %>%
    filter(protein == "CSRP3")

#define brackets for figure
brackets_csrp3 <- tibble(
    x = c(1, 1, 2),
    xend = c(2, 1, 2),
    y = c(1.7, 1.7, 1.7),
    yend = c(1.7, 1.5, 1.6)
)

#define fdr for figure
fdr_csrp3 <- tibble(
    x = 1.5,
    y = 1.8,
    label = "p=0.005"
)

# box plot of CSRP3
csrp3_fig <- csrp3 %>%
    ggplot(aes(x = fibertype, y = log2fc)) +
    geom_violin(aes(fill = fibertype), trim = TRUE, width = 1, linewidth = 0.5, alpha = 0.5, show.legend = FALSE) +
    geom_boxplot(width = 0.25, color = "black", fill = "white", alpha = 0.5) +
    geom_jitter(aes(color = sex), size = 3, width = 0, alpha = 0.5, stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(
        values = c("mhc1" = "#440154FF", "mhc2" = "#67CC5CFF"),
        labels = c("mhc1" = "Type I", "mhc2" = "Type II")
    ) +
    scale_color_manual(
        name = NULL,
        values = c("male" = "#FF7518", "female" = "#000000"),
        labels = c("male" = "Male", "female" = "Female")
    ) +
    scale_x_discrete(labels = c("Type I", "Type II")) +
    ggplot2::theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "top",
        legend.key.size = unit(0.5, "lines"),
        legend.box.margin = ggplot2::margin(-10, 0, -10, 0),
        axis.title.x = ggplot2::element_blank(),
        text = element_text(size = 8),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8)
    ) +
    scale_y_continuous(limits = c(-0.5, 1.9)) +
    geom_segment(data = brackets_csrp3, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = fdr_csrp3, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.75) +
    xlab("none") +
    ylab("CSRP3 log2fold change (post - pre)")


#Bar plot of CSRP3

#Compute log2fc from linear mixed model
emm_csrp3 <- emmeans(lm_csrp3_trial, ~ fibertype) %>%
    as.data.frame()

#define brackets for figure
brackets_csrp3 <- tibble(
    x = c(1, 1, 2),
    xend = c(2, 1, 2),
    y = c(1.0, 1.0, 1.0),
    yend = c(1.0, 0.6, 0.9)
)

#define p for figure
p_csrp3 <- tibble(
    x = 1.5,
    y = 1.1,
    label = "p=0.005"
)

##CSRP3 FIGURE##
csrp3_fig <- emm_csrp3 %>%
    ggplot(aes(x = fibertype, y = estimate, fill = fibertype)) +
    geom_col(position = position_dodge(width = 0.9),
             width = 0.9, color = NA, alpha = 0.7) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, position = position_dodge(width = 0.9),
                  linewidth = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(
        name = NULL,
        values = c("mhc1" = "#440154FF", "mhc2" = "#67CC5CFF"),
        labels = c("mhc1" = "Type I", "mhc2" = "Type II")
    ) +
    scale_x_discrete(labels = c("Type I", "Type II")) +
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
        text = element_text(size = 8),
        axis.text.x = element_text(color = "black", size = 8),
        axis.text.y = element_text(color = "black", size = 8),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
    ) +
    #scale_y_continuous(limits = c(-0.2, 1.2)) +
    geom_segment(data = brackets_csrp3, aes(x = x, xend = xend, y = y, yend = yend), size = 0.25, inherit.aes = FALSE) +
    geom_text(data = p_csrp3, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2.75) +
    xlab("none") +
    ylab("CSRP3 log2fold change (post - pre)") +
    ggtitle("Change in CSRP3 in humans")


ggsave(plot = csrp3_fig, here::here('figures/figure_5/csrp3.pdf'), height = 60, width = 60, units = "mm")
