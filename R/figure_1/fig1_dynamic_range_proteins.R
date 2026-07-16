# Load data and packages --------------------------------------------------
source(here::here("R/Library.R"))

#load data
data_long <- vroom::vroom(
    here::here("data/data_log2_long_normalized.csv")
)

#load raw data
data_raw <- vroom::vroom(
    here::here("data-raw/library_based_data_raw.csv")
)

# load Metadata:

metadata <- vroom::vroom(
    here::here("data-raw/Metadata_REX_pooled_fibers.csv")
)


# Setup of dataframe:

data_wrangled <- data_raw |>
    dplyr::select(!Protein.Group) |> #! means "don't"
    dplyr::select(!First.Protein.Description) |>
    dplyr::select(!Protein.Names) |>
    dplyr::select(!Protein.Ids) |>
    dplyr::filter(!duplicated(Genes)) |>
    dplyr::filter(!is.na(Genes) == TRUE) |>
    tibble::column_to_rownames("Genes") |>
    t() |>
    as.data.frame() |>
    tibble::rownames_to_column("column_names") |>
    dplyr::mutate(new_column_names = stringr::str_extract(column_names, "(?<=-)(.*?)(?=_)"))

metadata <- metadata |>
    dplyr::mutate(position = toupper(position)) # toupper changes to upper case

data_wrangled <- data_wrangled |>
    dplyr::rename("position" = new_column_names) |>
    dplyr::group_by(
        metadata |>
            dplyr::select(c(position,
                            sample_id))
    ) |>
    tibble::column_to_rownames("sample_id") |>
    dplyr::select(!position) |>
    dplyr::select(!column_names) |>
    t() |>
    as.data.frame()

#Make sample_id row names in metadata
metadata <- metadata %>%
    tibble::column_to_rownames("sample_id")

#create groups in metadata
metadata <- metadata %>%
    tidyr::unite(col = "group",
                 c(sex, fibertype), sep ="_", remove = FALSE, na.rm = FALSE)

#need to filter for proteins identified in 70% of samples in at least one group
data_filtered <- data_wrangled %>%
    PhosR::selectGrps(grps=metadata$group, percent = 0.7, n = 1)

#setup raw data in long format
data_long <- data_filtered %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    group_by(metadata %>%
                 rownames_to_column("sample_id") %>%
                 dplyr::select(!fibers) %>%
                 dplyr::select(!CSA)) %>%  #merge with metadata and make sample_id column for common identifyer
    pivot_longer(where(is.numeric), names_to = "gene", values_to = "expression")


# Group by protein
data_mean <- data_long %>%
    group_by(gene) %>%
    summarise(expression = sum(expression, na.rm = TRUE)) %>%
    arrange(desc(expression)) %>%
    mutate(rank = rank(-expression, ties.method = "first")) %>%
    mutate(percent = expression / sum(expression, na.rm = TRUE) * 100) %>%
    mutate(log2_percent = log2(percent)) %>%
    mutate(label = case_when(
        gene %in% c("PECAM1", "CD34", "MYH7", "MYH2") ~ gene,
        TRUE ~ ""
    )) %>%
    mutate(cell = case_when(
        label == "CD34" ~ "Satellite cells and FAPs",
        label == "PECAM1" ~ "Endothelial cells",
        TRUE ~ NA_character_  # Use NA for proteins without a cell association
    ))

# Define the desired percentages
desired_percentages <- c(10, 1, 0.1, 0.01, 0.001)

# Convert percentages to log2 scale
log2_breaks <- log2(desired_percentages)

# Create plot
dynamic_range_fig <- data_mean %>%
    ggplot(aes(x = rank, y = log2_percent)) +
    geom_point(alpha = 0.5, color = "black", size = 0.5) +
    scale_y_continuous(
        breaks = log2_breaks,  # Set y-axis breaks in log2 scale
        labels = desired_percentages  # Use desired percentages as labels
    ) +
    theme_classic() +
    theme(
        legend.position = "top",
        legend.key.size = unit(0.4, "lines"),
        text = element_text(size = 8),
        strip.text = element_text(size = 8),
        legend.title = element_blank()
    ) +
    geom_hline(yintercept = log2_breaks, linetype = "dashed", color = "gray") +
    geom_label_repel(
        data = data_mean %>%
            filter(label != ""),  # Ensure only labeled proteins are used for plotting
        aes(label = label, fill = cell),  # Use fill only for proteins with cell associations
        color = "black",
        size = 2,
        label.padding = 0.1,
        min.segment.length = 0.1,
        segment.size = 0.2,
        force = 20,
        max.overlaps = Inf
    ) +
    scale_fill_manual(
        values = c(
            "Satellite cells and FAPs" = "#F9E076",
            "Endothelial cells" = "#8b6c5c"
        ), # Only include actual cell types
        labels = c(
            "Satellite cells and FAPs" = "Satellite cells and FAPs",
            "Endothelial cells" = "Endothelial cells"
        ),
        breaks = c("Satellite cells and FAPs", "Endothelial cells"),
        na.translate = FALSE  # Exclude NA from the legend
    ) +
    xlab("Protein rank") +
    ylab("% total intensity (log2)")

ggsave(plot = dynamic_range_fig, here::here('figures/figure_1/figure1_dynamic_range.pdf'), height = 60, width = 60, units = "mm")
