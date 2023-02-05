suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(dplyr))

# Set paths and constants
set.seed(1234)
sample_cell_frac = 0.1

plate_to_focus <- "localhost220513100001_KK22-05-198_FactinAdjusted"

# input_data_dir <- file.path("..", "..", "..", "3.process-cfret-features", "data")
input_data_dir <- file.path("..", "UMAP", "temp")
file_suffix <- "_sc_norm_fs_cellprofiler.csv.gz_dropna.csv.gz"
cp_file <- file.path(
    input_data_dir, paste0(plate_to_focus, file_suffix)
)

output_figure_dir <- "figures"
cp_heatmap_file_noext <- file.path(output_figure_dir, "cp_complex_heatmap")

# Load data
cp_df <- readr::read_csv(
    cp_file,
    col_types = readr::cols(
        .default="d",
        "Metadata_WellRow" = "c",
        "Metadata_WellCol" = "c",
        "Metadata_heart_number" = "c",
        "Metadata_treatment" = "c",
        "Metadata_dose" = "c",
        "Metadata_ImageNumber" = "c",
        "Metadata_Plate" = "c",
        "Metadata_Well" = "c"
    )
) %>%  # Add dose numeric column
    dplyr::mutate(
        Metadata_dose_value = as.numeric(stringr::str_remove(Metadata_dose, "uM"))
    )

# Add cell count column
cell_count_df <- cp_df %>%
    dplyr::group_by(Metadata_Well) %>%
    dplyr::count()

cp_df <- cp_df %>%
    dplyr::left_join(cell_count_df, by = "Metadata_Well") %>%
    dplyr::rename(Metadata_cell_count_per_well = n)

print(dim(cp_df))
head(cp_df, 3)

# Setup heatmap colors
treatment_colors = c(
    "drug_x" = "#785EF0",
    "DMSO" = "#DC267F"
)

# From colorbrewer2.org
spectral_palette <- c(
    "#fff7ec",
    "#fee8c8",
    "#fdd49e",
    "#fdbb84",
    "#fc8d59",
    "#ef6548",
    "#d7301f",
    "#b30000",
    "#7f0000"
)

spectral_limits <- c(0, max(cp_df$Metadata_dose_value))

spectral_breaks <- seq(
    spectral_limits[1],
    spectral_limits[2],
    length.out = length(spectral_palette)
) 

dose_col <- circlize::colorRamp2(spectral_breaks, spectral_palette)

# Randomly sample cells
# Split metadata and feature data
cp_feature_df <- cp_df %>%
    dplyr::sample_frac(sample_cell_frac)

cp_metadata_df <- cp_feature_df %>%
    dplyr::select(tidyr::starts_with("Metadata"))
cp_meta_cols <- colnames(cp_metadata_df)

cp_feature_df <- cp_feature_df %>% dplyr::select(-!!cp_meta_cols)

# Calculate correlation matrix from feature data
cp_cor_matrix <- t(cp_feature_df) %>% cor()

print(dim(cp_cor_matrix))
head(cp_cor_matrix, 3)

# Generate heatmap
ht <- Heatmap(
    cp_cor_matrix,
    name = "Pearson\nCorrelation",
    column_dend_side = "top",
    
    clustering_method_columns = "average",
    clustering_method_rows = "average",
    
    top_annotation = HeatmapAnnotation(
        Treatment = cp_metadata_df$Metadata_treatment,
        CellCount = anno_barplot(
            cp_metadata_df$Metadata_cell_count_per_well,
            height = unit(1, "cm")
        ),
        Dose = cp_metadata_df$Metadata_dose_value,
        WellColumn = cp_metadata_df$Metadata_WellCol,
        
        col = list(
            Treatment = treatment_colors,
            Dose = dose_col
        )

    )
)

draw(ht)

# Save heatmap to file
pdf(paste0(cp_heatmap_file_noext, ".pdf"))
draw(ht)
dev.off()

png(paste0(cp_heatmap_file_noext, ".png"), width = 6.5, height = 6, units = "in", res = 500)
draw(ht)
dev.off()
