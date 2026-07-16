# Load data and packages --------------------------------------------------
source(here::here("R/Library.R"))

#load data
long_df <- vroom::vroom(here::here("data/data_log2_long_normalized.csv")) %>%
    rename(protein = gene)

#Load log2fold difference between sexes
log2fd_sex <- vroom::vroom(here::here("data/results_sex_keywords.csv"))

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

#Add annotations to keyword dataframe
df_keywords <- long_df %>%
    merge(annotations, by="protein", all.x = TRUE)

#saveRDS(df_keywords,
#here::here("data/data_long_keywords.rds"))

#Filter data for pre
df_pre <- df_keywords %>%
    filter(trial == "pre")

#create dataframe of mean expression of each protein pre
df_pre_mean <- long_df %>%
    filter(trial == "pre") %>%
    group_by(protein, group, fibertype, sex) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE)) %>%
    ungroup() %>%
    merge(annotations, by="protein", all.x = TRUE)

# INTERMEDIATE FILAMENT ---------------------------------------------------

#filter for intermediate filament proteins
if_proteins <- df_pre_mean %>%
    filter(
        grepl("intermediate filament organization|intermediate filament cytoskeleton organizaton|intermediate filament-based process", gobp, ignore.case = TRUE) |
            grepl("KRT", protein, ignore.case = TRUE) |
            protein == "NES"
    )

# Make sex a factor with male as reference
if_proteins$sex <- factor(if_proteins$sex, levels = c("male", "female"))

#linear mixed model of difference in intermediate filament protein abundance
lmm_if <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = if_proteins, REML = FALSE)
effects_if <- anova(lmm_if)

#calculate means for each fibertype for each sex
means_if <- emmeans(lmm_if, ~ sex | fibertype)
print(means_if)

#run linear model of sex for each fibertype
lm_if_sex <- contrast(means_if, method = "pairwise", by = "fibertype")
print(lm_if_sex)

#run linear model of fibertype for each sex
lm_if_fibertype <- contrast(means_if, method = "pairwise", by = "sex")
print(lm_if_fibertype)

#filter for intermediate filament protein log2fold difference between sexes
if_log2fd <- log2fd_sex %>%
    filter(
        grepl("intermediate filament organization|intermediate filament cytoskeleton organizaton|intermediate filament-based process", gobp, ignore.case = TRUE) |
            grepl("KRT", protein, ignore.case = TRUE) |
            protein == "NES"
    )

#define text for figure
p_if <- tibble(
    fibertype = c("mhc1", "mhc2"),
    x = c(1, 2),
    y = c(0.4, 0.4),
    label = c("p=0.019", "p=0.006")
)

#log2fold difference figure (box plot)
if_fig <- if_log2fd %>%
    ggplot(aes(x = fibertype, y = logFC, fill = fibertype)) +
    geom_violin(trim = TRUE, width = 1, linewidth = 0.25, alpha = 0.5) +
    geom_boxplot(width = 0.25, linewidth = 0.25, color = "black", fill = "white", alpha = 0.5, outlier.size = 0, outlier.stroke = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.25) +
    scale_fill_manual(values=c("#440154FF", "#67CC5CFF"))+
    scale_x_discrete(labels = c("mhc1" = "Type I", "mhc2" = "Type II")) +
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
    scale_y_continuous(limits = c(-1, 0.5)) +
    geom_text(data = p_if, aes(x = x, y = y, label = label), inherit.aes = FALSE, size = 2) +
    xlab("") +
    ylab("Log2fold difference (male - female)") +
    ggtitle("Intermediate filament proteins")

ggsave(plot = if_fig, here::here('figures/figure_3/if_pre_figure.pdf'), height = 45, width = 60, units = "mm")
