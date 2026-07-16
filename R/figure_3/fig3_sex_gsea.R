source(here::here("R/Library.R"))

##SEX DIFFERENCES PRE##

#open gsea_bp files of sex differences

mhc1_sex_bp <- vroom::vroom(here::here("results/gsea_bp_mhc1_sex.csv")) %>%
    mutate(fibertype = "mhc1") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "male",
        enrichmentScore < 0 ~ "female"))

mhc2_sex_bp <- vroom::vroom(here::here("results/gsea_bp_mhc2_sex.csv")) %>%
    mutate(fibertype = "mhc2") %>%
    mutate(regulated = case_when(
        enrichmentScore > 0 ~ "male",
        enrichmentScore < 0 ~ "female"))

#merge BP data frames
sex_bp <- rbind(mhc1_sex_bp, mhc2_sex_bp)
sex_bp$regulated <- factor(sex_bp$regulated, levels = c("female", "male"))

#write.csv(sex_bp,
#here::here("results/gsea_bp_sex.csv"))

#Enrichment plot of BP for sex differences
gsea <- sex_bp %>%
    dplyr::filter(Description %in% c(
        "intermediate filament organization",
        "mitochondrial gene expression",
        "mitochondrial translation",
        "energy derivation by oxidation of organic compounds",
        "long-chain fatty acid transport",
        "lipoprotein metabolic process",
        "skeletal muscle contraction",
        "mRNA destabilization",
        "extracellular matrix organization",
        "tissue remodeling",
        "ribose phosphate metabolic process",
        #"cell adhesion",
        "proton transmembrane transport",
        "ATP metabolic process",
        "regulation of fatty acid metabolic process",
        "response to external stimulus",
        "translation",
        "negative regulation of peptidase activity",
        "regulation of cellular catabolic process"
    )) |>
    ggplot(
        aes(
            x = fibertype,
            y = Description,
            color = qvalue,
            size = abs(NES)
        )
    ) +
    geom_point() +
    facet_wrap(
        ~regulated,
        labeller = labeller(
            regulated = c(
                "female" = "Enriched in females",
                "male" = "Enriched in males"
            )
        )
    ) +
    scale_x_discrete(
        breaks = c("mhc1", "mhc2"),
        labels = c("Type I", "Type II")
    ) +
    theme_bw() +
    theme(
        text = ggplot2::element_text(size = 7),
        strip.text = ggplot2::element_text(size = 6),
        legend.key.size = ggplot2::unit(2, units = "mm"),

        legend.position = "right",
        legend.location = "plot",
        legend.justification.right = "top",
        legend.box = "vertical",
        legend.box.just = "top",

        legend.box.margin = margin(0, -10, 0, -10),
        plot.margin = margin(1, 5, 1, 1)
    ) +
    xlab("") +
    ylab("") +
    labs(
        color = "FDR",
        size = "NES"
    )
ggsave(plot = gsea, here::here('figures/figure_3/sex_gsea.pdf'), height = 80, width = 110, units = "mm")


# LINEAR MIXED MODEL OF ENRICHED TERMS ------------------------------------

#load pre data
df_pre <- vroom::vroom(here::here("data/data_long_keywords.csv")) %>%
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

#create dataframes for enriched terms
df_translation <- df_pre_mean %>%
    filter(grepl("translation",
                 gobp, ignore.case = TRUE) |
               protein %in% c("YBX1","GIGYF2","MRPL38","TNKS1BP1","MRPL4","MTFMT","PELO","RPS17","CAPRIN1",
                          "RPL3","MRPL43","GAPDH","RBM3","NDUFA7","CELF1","RPS29","ROCK1","RPL38",
                          "RPS7","MRPL19","MRPS31","RPS28","KHDRBS1","RPS13","GADD45GIP1","RPS15A",
                          "MRPL53","PNPT1","MRPL3","EEF1A2","RPS26","SARNP","LRRC47","EIF2S2",
                          "MRPS25","MRPS14","RPS14","RPL22","EEF1A1","MTPN","TRMT10C","PABPC4",
                          "RPS19","MRPS7","EEF1G","RPS11","TACO1","MRPL12","YARS1","LSM14A","RPS10",
                          "RPL18A","RPS20","TSFM","PDF","ACO1","RPL24","RPS27","MRPS17","SARS1",
                          "YBX3","AIMP2","TUFM","MRPS30","AARSD1","SYNCRIP","RPL12","SHMT2","PTCD3",
                          "MTOR","EIF3G","RPS25","MRPL37","EIF4A1","FECH","MRPL40","HNRNPR","EEF1D",
                          "MRPL30","EIF3F","MRPS23","SSB","MRPL46","MRPL48","CIRBP","EIF3A","MRPS28",
                          "RPL23","RPLP2","LTN1","LSM14B","EIF3L","RPS3A","SERBP1","EIF3J","PKM",
                          "MRPS5","ETF1","RWDD1","NARS1","CARS1","CALR","EIF5A","MTIF2","LARP4",
                          "GATC","MTRF1L","RPL10","MCTS1","EEF1B2","EEF2","PDE12","LARS2","RBM8A",
                          "MRPL13","TRAP1","EIF3M","YARS2","WARS1","EIF4G1","EIF4G2","TPR","RCC1L",
                          "EIF2B5","RPL3L","RPL8","NCL","EIF4A2","RPS9","NARS2","EEF2K","LRPPRC",
                          "C1QBP","MRPL50","ABCE1","IARS1","MAGOH","CSNK2A1","RPL27A","MSI2","MRPL47",
                          "RPSA","RPS23","FARSA","DHPS"))

df_mrna <- df_pre_mean %>%
    filter(grepl("mRNA destabilization",
                 gobp, ignore.case = TRUE) |
               protein %in% c("YBX1", "GIGYF2", "TNKS1BP1", "CAPRIN1", "CELF1", "ROCK1", "PNPT1"))

df_fa <- df_pre_mean %>%
    filter(grepl("regulation of fatty acid metabolic process | long-chain fatty acid transport | lipoprotein metabolic process",
                 gobp, ignore.case = TRUE) |
               protein %in% c("APOA4","PDK1","PLIN5","ANXA1","FABP5","ADIPOQ","APOC1","FABP3","APOC3",
                       "ERLIN2","SNCA","PDK4","ACSL3","APOE","CD36","PLIN2","SLC2A1","FABP4",
                       "ITGAV","ATG7","HHATL","APOB","APOA2","NMT1","APOA1","PIGK","APOD"))

#linear mixed model of difference in translation
lmm_translation <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = df_translation, REML = FALSE)
effects_translation <- anova(lmm_translation)

#linear mixed model of difference in mRNA destabilization
lmm_mrna <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = df_mrna, REML = FALSE)
effects_mrna <- anova(lmm_mrna)

#linear mixed model of difference in fatty acid metabolism
lmm_fa <- lmer(mean_expression ~ fibertype * sex + (1 | protein), data = df_fa, REML = FALSE)
effects_fa <- anova(lmm_fa)
