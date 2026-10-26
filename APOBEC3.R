library(tidyverse)
setwd("C:/Users/lenovo/Desktop/Dossiers/Projets de Recherche/Mpox/Seq Mpox Projet rim/")
# 1. Charger le fichier TSV exporté depuis Nextclade
file_path <- "Mpox2025.tsv"
df <- read_tsv(file_path, show_col_types = FALSE)

# Déterminer le nom de la colonne contenant les substitutions
mut_col <- if ("substitutions" %in% colnames(df)) "substitutions" else "mutations"

# 2. Fonction d'analyse de la signature APOBEC3 sur une chaîne de mutations
analyze_apobec3 <- function(sub_string) {
  if (is.na(sub_string) || nchar(trimws(sub_string)) == 0) {
    return(tibble(
      total_mutations = 0,
      apobec3_mutations = 0,
      other_mutations = 0,
      pct_apobec3 = 0.0
    ))
  }
  
  # Séparation des mutations (ex: "C1234T,G5678A,A9101G")
  mut_list <- str_split(sub_string, ",")[[1]] %>% str_trim()
  mut_list <- mut_list[nchar(mut_list) >= 3]
  
  total <- length(mut_list)
  if (total == 0) {
    return(tibble(
      total_mutations = 0,
      apobec3_mutations = 0,
      other_mutations = 0,
      pct_apobec3 = 0.0
    ))
  }
  
  # Extraction de la base initiale (première lettre) et mutée (dernière lettre)
  ref_bases <- str_sub(mut_list, 1, 1) %>% toupper()
  alt_bases <- str_sub(mut_list, -1, -1) %>% toupper()
  
  # Identification des mutations C->T et G->A
  is_apobec3 <- (ref_bases == "C" & alt_bases == "T") | (ref_bases == "G" & alt_bases == "A")
  
  apobec3_count <- sum(is_apobec3)
  other_count <- total - apobec3_count
  pct_apobec3 <- round((apobec3_count / total) * 100, 2)
  
  return(tibble(
    total_mutations = total,
    apobec3_mutations = apobec3_count,
    other_mutations = other_count,
    pct_apobec3 = pct_apobec3
  ))
}

# 3. Application du calcul sur l'ensemble des échantillons
results <- df %>%
  select(seqName, all_of(mut_col)) %>%
  rename(sub_data = all_of(mut_col)) %>%
  mutate(metrics = map(sub_data, analyze_apobec3)) %>%
  unnest(metrics) %>%
  select(-sub_data)

# 4. Affichage d'un aperçu des résultats dans la console
print(head(results, 15))

# 5. Export des résultats au format CSV
write_csv(results, "apobec3_summary_Mpox2025.csv")
cat("\n[OK] Analyse terminée. Fichier 'apobec3_summary_Mpox2025.csv' généré avec succès.\n")




##############APOBEC3 selon le temps##################"
library(tidyverse)
library(lubridate)

# 1. Load both CSV files
results <- read_csv("apobec3_summary_Mpox2025.csv")
meta <- read_csv2("Mpox_Date_2025.csv") # read_csv2 handles semicolon-separated values (;)

# Clean potential leading/trailing whitespace in column headers
colnames(meta) <- str_trim(colnames(meta))

# 2. Merge metadata and parse collection dates (DD/MM/YYYY format)
df_temporal <- results %>%
  bind_cols(meta %>% select(Date, clade, Host, Villes)) %>%
  mutate(
    date = dmy(Date), # Automatically converts DD/MM/YYYY format
    year_decimal = decimal_date(date)
  ) %>%
  filter(!is.na(date))

# 3. Generate plot with a pure white background (Figure 3)
fig3 <- ggplot(df_temporal, aes(x = date, y = apobec3_mutations)) +
  geom_point(aes(color = pct_apobec3, size = total_mutations), alpha = 0.8) +
  geom_smooth(method = "lm", color = "#D55E00", fill = "#D55E00", alpha = 0.2, se = TRUE) +
  scale_color_viridis_c(option = "viridis", name = "APOBEC3 (%)") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(
    title = "Figure 3: Temporal Accumulation of APOBEC3-Mediated Mutations",
    subtitle = "Linear increase in C->T / G->A mutations over time (Côte d'Ivoire, 2025)",
    x = "Collection Date",
    y = "Number of APOBEC3 Mutations (C->T & G->A)",
    size = "Total Mutations"
  ) +
  theme_classic(base_size = 12) + # Clean classic theme (white background, no gray grid)
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA),
    plot.title        = element_text(face = "bold", size = 13),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    legend.position   = "right",
    axis.line         = element_line(color = "black")
  )

# 4. Display and export image with explicit white background
print(fig3)
ggsave("Temporal_Accumulation.png", fig3, width = 9, height = 6, dpi = 300, bg = "white")


################Par Ville#######################


# ==============================================================================
# R PIPELINE: Temporal and Geographic Analysis of APOBEC3 Mutations (Mpox 2025)
# Publication-ready layout: Clean white background & English labels
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, # Data manipulation (dplyr, readr) and plotting (ggplot2)
  lubridate, # Date handling
  patchwork, # Plot combination
  ggpubr     # Statistical annotations
)

# ------------------------------------------------------------------------------
# 2. Import and Clean Data
# ------------------------------------------------------------------------------
cat("\n[1/5] Importing CSV files...\n")

# Load APOBEC3 summary results (comma-separated)
results <- read_csv("apobec3_summary_Mpox2025.csv", show_col_types = FALSE)

# Load metadata (semicolon-separated)
meta <- read_csv2("Mpox_Date_2025.csv", show_col_types = FALSE)

# Clean column headers
colnames(meta) <- str_trim(colnames(meta))

# Merge datasets and format dates/locations
df_data <- results %>%
  bind_cols(meta %>% select(Date, clade, Host, Villes)) %>%
  mutate(
    # Clean string spaces
    Villes = str_trim(Villes),
    clade  = str_trim(clade),
    
    # Parse collection dates (DD/MM/YYYY format)
    date = dmy(Date),
    year_decimal = decimal_date(date),
    days_since_start = as.numeric(date - min(date, na.rm = TRUE))
  ) %>%
  filter(!is.na(date))

# Standardize regional groupings (e.g. Abidjan District)
df_data <- df_data %>%
  mutate(
    District_Region = case_when(
      str_detect(Villes, "(?i)Abidjan|Abobo|PORT BOUET") ~ "Abidjan District",
      TRUE ~ Villes
    )
  )

cat("-> Number of analyzed sequences:", nrow(df_data), "\n")

# ------------------------------------------------------------------------------
# 3. Statistical Analysis: Linear Regression
# ------------------------------------------------------------------------------
cat("\n[2/5] Calculating temporal linear regression...\n")

lm_model <- lm(apobec3_mutations ~ days_since_start, data = df_data)
lm_summary <- summary(lm_model)

r2_val <- lm_summary$r.squared
p_val  <- lm_summary$coefficients[2, 4]
slope  <- lm_summary$coefficients[2, 1]

cat(sprintf("   - R²: %.4f\n", r2_val))
cat(sprintf("   - p-value: %.4e\n", p_val))
cat(sprintf("   - Accumulation rate: %.4f APOBEC3 mutations / day\n", slope))

# ------------------------------------------------------------------------------
# 4. Generate Figures (Pure White Background / English Labels)
# ------------------------------------------------------------------------------
cat("\n[3/5] Generating graphical representations...\n")

# --- PANEL A: Temporal Accumulation by Sampling Location ---
p_temporal <- ggplot(df_data, aes(x = date, y = apobec3_mutations)) +
  geom_smooth(method = "lm", color = "#222222", fill = "gray85", alpha = 0.4, se = TRUE, linetype = "dashed") +
  geom_point(aes(color = District_Region, size = total_mutations), alpha = 0.85) +
  scale_size_continuous(range = c(3, 7), name = "Total Mutations") +
  scale_color_brewer(palette = "Set1", name = "Sampling Location") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %Y") +
  labs(
    title = "A. Temporal accumulation of APOBEC3-mediated mutations",
    subtitle = sprintf("R² = %.2f | p = %.4f | Rate = +%.3f mut/day", r2_val, p_val, slope),
    x = "Collection Date",
    y = "Number of APOBEC3 Mutations (C->T / G->A)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA),
    plot.title        = element_text(face = "bold", size = 13),
    axis.text.x       = element_text(angle = 45, hjust = 1),
    legend.position   = "right",
    axis.line         = element_line(color = "black")
  )

# --- PANEL B: Percentage of APOBEC3 Mutations by Sampling Site ---
p_box <- ggplot(df_data, aes(x = reorder(Villes, pct_apobec3, FUN = median), y = pct_apobec3)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.2, fill = "#0072B2", color = "black") +
  geom_jitter(aes(color = clade), width = 0.15, size = 3, alpha = 0.8) +
  coord_flip() +
  scale_color_viridis_d(option = "plasma", name = "Clade") +
  labs(
    title = "B. Proportion of APOBEC3 activity (% C->T / G->A) by collection site",
    x = "Sampling Location / City",
    y = "APOBEC3 Mutations (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key        = element_rect(fill = "white", color = NA),
    plot.title        = element_text(face = "bold", size = 13),
    axis.line         = element_line(color = "black")
  )

# ------------------------------------------------------------------------------
# 5. Combine and Export Figures
# ------------------------------------------------------------------------------
cat("\n[4/5] Assembling final figure...\n")

final_figure <- p_temporal / p_box + 
  plot_annotation(
    title = "Genomic Surveillance of Mpox in Côte d'Ivoire (2025)",
    subtitle = "Spatiotemporal analysis of APOBEC3 mutational signature and viral evolution",
    caption = "Source: Institut Pasteur de Côte d'Ivoire (IPCI)",
    theme = theme(
      plot.background = element_rect(fill = "white", color = NA),
      plot.title      = element_text(face = "bold", size = 16),
      plot.subtitle   = element_text(size = 12, color = "gray30")
    )
  )

# Save High-Resolution Figures with explicit white background
ggsave("Figure_Mpox_APOBEC3_Temporal_Geo_2025.png", final_figure, width = 11, height = 10, dpi = 300, bg = "white")
ggsave("Figure_Mpox_APOBEC3_Temporal_Geo_2025.pdf", final_figure, width = 11, height = 10, bg = "white")

cat("\n[5/5] Exporting merged dataset...\n")
write_csv(df_data, "Mpox_2025_APOBEC3_Merged_Analysis.csv")

cat("\nPipeline finished successfully! Output files:\n")
cat(" - Figure_Mpox_APOBEC3_Temporal_Geo_2025.png\n")
cat(" - Figure_Mpox_APOBEC3_Temporal_Geo_2025.pdf\n")
cat(" - Mpox_2025_APOBEC3_Merged_Analysis.csv\n\n")