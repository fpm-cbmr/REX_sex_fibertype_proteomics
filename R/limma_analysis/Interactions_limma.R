source(here::here("R/Library.R"))

#Load results of regulated proteins (main effect)
regulated_proteins <- vroom::vroom(here::here("results/trial_main_effect.csv")) %>%
    select(!1) %>%
    filter(q < 0.05)

#load log2fold change data
data_log2fc <- vroom::vroom(here::here("data/data_log2fc.csv")) %>%
    rename_with(~ "protein", colnames(.)[1]) %>%
    filter(protein %in% regulated_proteins$protein) %>%
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

#create summarized experiment
se_log2fc <- SummarizedExperiment(
    assays = list(counts = as.matrix(data_log2fc)),
    colData = DataFrame(metadata_log2fc)
)

# TRIAL x SEX INTERACTION -------------------------------------------------

#create design matrix for interaction between log2fold change in males and females
design_sex <- model.matrix(~0 + se_log2fc$sex)
colnames(design_sex) <- c("female", "male")

#define comparison using contrast
contrast_sex <- makeContrasts(male - female, levels = design_sex)

#correlation between samples from same subject
correlation_sex <- duplicateCorrelation(assay(se_log2fc), design_sex, block = se_log2fc$subject)

#Use Bayes statistics to assess differential expression interaction between males and females
fit_sex <- eBayes(lmFit(assay(se_log2fc), design_sex, block = se_log2fc$subject, correlation = correlation_sex$consensus.correlation))
ebayes_sex <- eBayes(contrasts.fit(fit_sex, contrast_sex))

#extract results
results_sex <- topTable(ebayes_sex, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))

# TRIAL x FIBERTYPE INTERACTION -------------------------------------------

#create design matrix for interaction between log2fold change in type I and type II fibers
design_fibertype <- model.matrix(~0 + se_log2fc$fibertype)
colnames(design_fibertype) <- c("typeI", "typeII")

#define comparison using contrast
contrast_fibertype <- makeContrasts(typeII - typeI, levels = design_fibertype)

#correlation between samples from same subject
correlation_fibertype <- duplicateCorrelation(assay(se_log2fc), design_fibertype, block = se_log2fc$subject)

#Use Bayes statistics to assess differential expression interaction between type I and type II fibers
fit_fibertype <- eBayes(lmFit(assay(se_log2fc), design_fibertype, block = se_log2fc$subject, correlation = correlation_fibertype$consensus.correlation))
ebayes_fibertype <- eBayes(contrasts.fit(fit_fibertype, contrast_fibertype))

#extract results
results_fibertype <- topTable(ebayes_fibertype, coef = 1, number = Inf, sort.by = "logFC") %>%
    dplyr::mutate(protein = row.names(.),
                  q = qvalue(.$P.Value)$qvalues,
                  "-log10p" = -log10(.$P.Value),
                  regulated_q = ifelse(q < 0.05, "+", "")
    )%>%
    arrange(desc(logFC))






