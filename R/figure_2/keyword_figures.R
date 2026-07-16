source(here::here("R/Library.R"))

#load pre data
df_pre <- readRDS(here::here("data/data_long_keywords.rds")) %>%
    filter(trial == "pre")

#Add keywords and GO terms
annotations <- vroom::vroom(here::here("data/keywords.csv")) %>%
    dplyr::rename_with(snakecase::to_snake_case) %>%
    dplyr::select(c("gene_names", "keywords", "gene_ontology_biological_process", "gene_ontology_cellular_component", "gene_ontology_molecular_function")) %>%
    dplyr::rename(gobp = gene_ontology_biological_process,
                  gocc = gene_ontology_cellular_component,
                  gomf = gene_ontology_molecular_function,
                  protein = gene_names) %>%
    dplyr::mutate(protein = gsub("\\ .*","", protein)) %>%
    dplyr::mutate(protein = make.names(protein, unique=TRUE), protein)

#create dataframe of mean expression of each protein pre
df_pre_mean <- df_pre %>%
    group_by(protein, group, fibertype, sex) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
    ungroup() %>%
    merge(annotations, by="protein", all.x = TRUE)

#Load log2fold difference between fiber types
log2fd_m <- vroom::vroom(here::here("results/differential_expression_male_fibertypes.csv")) %>%
    mutate(sex = "male")
log2fd_f <- vroom::vroom(here::here("results/differential_expression_female_fibertypes.csv")) %>%
    mutate(sex = "female")

#create keyword data frame for log2fd for fibertypes
log2fd_df <- rbind(log2fd_m, log2fd_f) %>%
    merge(annotations, by="protein", all.x = TRUE)

#write.csv(log2fd_df,
#here::here("data/results_fibertypes_keywords.csv"))

# CYTOSOLIC RIBOSOME ------------------------------------------

# Make fibertype a factor with type II as reference
df_pre_mean$fibertype <- factor(df_pre_mean$fibertype, levels = c("mhc2", "mhc1"))

#filter for cytosolic ribosome
ribo_proteins <- df_pre_mean %>%
    filter(
        grepl("cytosolic ribosome", gocc, ignore.case = TRUE))

#linear mixed model of difference in intermediate filament protein abundance
lmm_ribo <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = ribo_proteins, REML = FALSE)
effects_ribo <- anova(lmm_ribo)

#calculate means for each fibertype for each sex
means_ribo <- emmeans(lmm_ribo, ~ sex | fibertype)
print(means_ribo)

#run linear model of sex for each fibertype
lm_ribo_sex <- contrast(means_ribo, method = "pairwise", by = "fibertype")
lm_ribo_sex_df <- summary(lm_ribo_sex)


#run linear model of fibertype for each sex
lm_ribo_fibertype <- contrast(means_ribo, method = "pairwise", by = "sex")
lm_ribo_fibertype_df <- summary(lm_ribo_fibertype)

# MITOCHONDRIA ------------------------------------------------------------

#filter for mitochondria
mito_proteins <- df_pre_mean %>%
    filter(
        grepl("mitochondrion", gocc, ignore.case = TRUE))

#linear mixed model of difference in intermediate filament protein abundance
lmm_mito <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = mito_proteins, REML = FALSE)
effects_mito <- anova(lmm_mito)

#calculate means for each fibertype for each sex
means_mito <- emmeans(lmm_mito, ~ sex | fibertype)
print(means_mito)

#run linear model of sex for each fibertype
lm_mito_sex <- contrast(means_mito, method = "pairwise", by = "fibertype")
lm_mito_sex_df <- summary(lm_mito_sex)

#run linear model of fibertype for each sex
lm_mito_fibertype <- contrast(means_mito, method = "pairwise", by = "sex")
lm_mito_fibertype_df <- summary(lm_mito_fibertype)


# MITOCHONDRIA AND RIBOSOME FIGURES -----------------------------------------------------------------

#filter for cc
ribo_log2fd <- log2fd_df %>%
    filter(
        grepl("cytosolic ribosome", gocc, ignore.case = TRUE)) %>%
    mutate(cc = "ribosome")

mito_log2fd <- log2fd_df %>%
    filter(
        grepl("mitochondrion", gocc, ignore.case = TRUE)) %>%
    mutate(cc = "mitochondria")

#combine data frames for figure
log2fd_cc <- rbind(ribo_log2fd, mito_log2fd)


#define text for figure
p_cc <- tibble(
    sex = c("female", "female",
            "male", "male"),
    x = c(1, 2,
          1, 2),
    y = c(3.5, 1.0,
          5.0, 1.25),
    label = c("p<0.001", "p<0.001",
              "p<0.001", "p<0.001")
)

##MITOCHONDRIA AND RIBOSOME FIGURE##
cc_fig <- log2fd_cc %>%
    ggplot(aes(x = cc, y = logFC, fill = cc)) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.25, alpha = 0.5) +
    geom_boxplot(width = 0.25, linewidth = 0.25, color = "black", fill = "white", alpha = 0.5, outlier.size = 0, outlier.stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("mitochondria" = "Mitochondria", "ribosome" = "Ribosome")) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    facet_wrap(~sex, labeller = labeller(sex = c("female" = "Females", "male" = "Males"))) +
    #scale_y_continuous(limits = c(-1, 0.5)) +
    geom_text(data = p_cc, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold difference (type II - type I)") +
    ggtitle("Mitochondria and ribosomes")

ggsave(plot = cc_fig, here::here('figures/figure_2/mito_ribo_figure.pdf'), height = 60, width = 90, units = "mm")

# GLYCOLYSIS --------------------------------------------------------------
#filter for glycolysis
glyco_proteins <- df_pre_mean %>%
    filter(
        grepl("glycolysis", gobp, ignore.case = TRUE))

#linear mixed model of difference in intermediate filament protein abundance
lmm_glyco <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = glyco_proteins, REML = FALSE)
effects_glyco <- anova(lmm_glyco)

#calculate means for each fibertype for each sex
means_glyco <- emmeans(lmm_glyco, ~ sex | fibertype)
print(means_glyco)

#run linear model of sex for each fibertype
lm_glyco_sex <- contrast(means_glyco, method = "pairwise", by = "fibertype")
lm_glyco_sex_df <- summary(lm_glyco_sex)

#run linear model of fibertype for each sex
lm_glyco_fibertype <- contrast(means_glyco, method = "pairwise", by = "sex")
lm_glyco_fibertype_df <- summary(lm_glyco_fibertype)


# FATTY ACID METABOLISM ---------------------------------------------------

#filter for fatty acid metabolism
fa_proteins <- df_pre_mean %>%
    filter(
        grepl("fatty acid metabolic process", gobp, ignore.case = TRUE))

#linear mixed model of difference in intermediate filament protein abundance
lmm_fa <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = fa_proteins, REML = FALSE)
effects_fa <- anova(lmm_fa)

#calculate means for each fibertype for each sex
means_fa <- emmeans(lmm_fa, ~ sex | fibertype)
print(means_fa)

#run linear model of sex for each fibertype
lm_fa_sex <- contrast(means_fa, method = "pairwise", by = "fibertype")
lm_fa_sex_df <- summary(lm_fa_sex)

#run linear model of fibertype for each sex
lm_fa_fibertype <- contrast(means_fa, method = "pairwise", by = "sex")
lm_fa_fibertype_df <- summary(lm_fa_fibertype)


#filter for bp
glyco_log2fd <- log2fd_df %>%
    filter(
        grepl("Glycolysis", gobp, ignore.case = TRUE)) %>%
    mutate(bp = "Glycolysis")

fa_log2fd <- log2fd_df %>%
    filter(
        grepl("fatty acid metabolic process", gobp, ignore.case = TRUE)) %>%
    mutate(bp = "Fatty acid metabolism")

#combine data frames for figure
log2fd_bp <- rbind(glyco_log2fd, fa_log2fd)

#Compute emmeans from linear mixed model
emm_glyco <- emmeans(lm_glyco_fibertype, ~ sex | sex) %>%
    as.data.frame() %>%
    mutate(bp = "Glycolysis") %>%
    mutate(
        fibertype = case_when(
            estimate < 0 ~ "type1",
            estimate > 0 ~ "type2"))


emm_fa <- emmeans(lm_fa_fibertype, ~ sex | sex) %>%
    as.data.frame() %>%
    mutate(bp = "Fatty acid metabolism") %>%
    mutate(
        fibertype = case_when(
            estimate < 0 ~ "type1",
            estimate > 0 ~ "type2"))

emm_bp <- rbind(emm_glyco, emm_fa)


#define text for figure
p_bp <- tibble(
    sex = c("female", "female",
            "male", "male"),
    x = c(1, 2,
          1, 2),
    y = c(0.8, 1.3,
          0.4, 1.6),
    label = c("p<0.001", "p<0.001",
              "p<0.001", "p<0.001")
)

##FIGURE OF GLYCOLYSIS AND FATTY ACID METABOLISM##
bp_fig <- log2fd_bp %>%
    ggplot(aes(x = bp, y = logFC, fill = bp)) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.25, alpha = 0.5) +
    geom_boxplot(width = 0.25, linewidth = 0.25, color = "black", fill = "white", alpha = 0.5, outlier.size = 0, outlier.stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("Fatty acid metabolism" = "Fatty acid metabolism", "Glycolysis" = "Glycolysis")) +
    theme_bw() +
    theme(
        panel.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_blank(),
        legend.position = "none",
        text = element_text(size = 6),
        axis.text.x = element_text(color = "black", size = 6),
        axis.text.y = element_text(color = "black", size = 6),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 8),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5)
    ) +
    facet_wrap(~sex, labeller = labeller(sex = c("female" = "Females", "male" = "Males"))) +
    #scale_y_continuous(limits = c(-1, 0.5)) +
    geom_text(data = p_bp, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold difference (type II - type I)") +
    ggtitle("Glycolysis and fatty acid metabolism")

ggsave(plot = bp_fig, here::here('figures/figure_2/glyco_fa_figure.pdf'), height = 60, width = 90, units = "mm")

