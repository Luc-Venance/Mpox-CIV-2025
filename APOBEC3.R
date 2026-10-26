
# 1. Charger les packages ------------------------------------------------------

packages <- c("tidyverse", "lubridate", "patchwork", "viridis")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

library(tidyverse)
library(lubridate)
library(patchwork)
library(viridis)

# 2. Définir le répertoire de travail -----------------------------------------

setwd(
  "C:/Users/lenovo/Desktop/Dossiers/Projets de Recherche/Mpox/Seq Mpox Projet rim/"
)

# 3. Importer les fichiers CSV ------------------------------------------------

results <- read_csv(
  "apobec3_summary_Mpox2025.csv",
  show_col_types = FALSE
)

meta <- read_csv2(
  "Mpox_Date_2025.csv",
  show_col_types = FALSE
)

# Nettoyer les noms de colonnes
names(results) <- str_trim(names(results))
names(meta) <- str_trim(names(meta))

# Vérifier la correspondance des lignes
if (nrow(results) != nrow(meta)) {
  stop(
    "Erreur : les deux fichiers n'ont pas le même nombre de lignes. ",
    "Vérifiez la correspondance des séquences avant de continuer."
  )
}

# Vérifier la présence des colonnes requises
required_results <- c(
  "apobec3_mutations",
  "pct_apobec3",
  "total_mutations"
)

required_meta <- c("Date", "clade", "Host", "Villes")

if (!all(required_results %in% names(results))) {
  stop(
    "Colonnes manquantes dans le fichier APOBEC3 : ",
    paste(setdiff(required_results, names(results)), collapse = ", ")
  )
}

if (!all(required_meta %in% names(meta))) {
  stop(
    "Colonnes manquantes dans le fichier de métadonnées : ",
    paste(setdiff(required_meta, names(meta)), collapse = ", ")
  )
}

# 4. Fusionner et nettoyer les données ----------------------------------------

df_data <- bind_cols(
  results,
  meta %>% select(Date, clade, Host, Villes)
) %>%
  mutate(
    date = dmy(str_trim(as.character(Date))),
    Villes = str_trim(as.character(Villes)),
    clade = str_trim(as.character(clade)),
    apobec3_mutations = as.numeric(apobec3_mutations),
    pct_apobec3 = as.numeric(pct_apobec3),
    total_mutations = as.numeric(total_mutations)
  ) %>%
  filter(
    !is.na(date),
    !is.na(apobec3_mutations),
    !is.na(pct_apobec3),
    !is.na(total_mutations),
    !is.na(Villes),
    Villes != ""
  ) %>%
  mutate(
    Site = case_when(
      str_detect(
        Villes,
        regex("Abidjan|Abobo|Port[- ]Bouet", ignore_case = TRUE)
      ) ~ "Abidjan District",
      TRUE ~ Villes
    )
  )

if (nrow(df_data) < 3 || n_distinct(df_data$date) < 2) {
  stop("Données insuffisantes pour réaliser la régression temporelle.")
}

cat("Nombre de séquences analysées :", nrow(df_data), "\n")

# 5. Calculer la régression linéaire ------------------------------------------

modele <- lm(
  apobec3_mutations ~ as.numeric(date),
  data = df_data
)

resume_modele <- summary(modele)

r2 <- resume_modele$r.squared
p_value <- coef(resume_modele)[2, 4]
pente_jour <- unname(coef(modele)[2])

subtitle_text <- sprintf(
  "Linear regression: R² = %.3f | p = %.4g | slope = %.3f mutations/day",
  r2,
  p_value,
  pente_jour
)

cat("R² :", round(r2, 3), "\n")
cat("p-value :", signif(p_value, 4), "\n")
cat(
  "Pente :", round(pente_jour, 3),
  "mutations par jour\n"
)

# 6. Paramètres communs aux deux panneaux -------------------------------------

theme_commun <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12
    ),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_blank(), # SUPPRIMER le titre de l'axe Y
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    panel.background = element_rect(
      fill = "white",
      color = NA
    )
  )

# 7. PANEL A : évolution temporelle selon la ville ----------------------------

p_location <- ggplot(
  df_data,
  aes(
    x = date,
    y = apobec3_mutations
  )
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "#D55E00",
    fill = "#D55E00",
    alpha = 0.20,
    se = TRUE
  ) +
  geom_point(
    aes(
      color = Site,
      size = total_mutations
    ),
    alpha = 0.85
  ) +
  scale_color_brewer(
    palette = "Set1",
    name = "Sampling location"
  ) +
  scale_size_continuous(
    range = c(2, 7),
    name = "Total mutations"
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %Y",
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "A. Temporal trend by sampling location",
    subtitle = subtitle_text,
    x = "Collection date",
    y = NULL
  ) +
  theme_commun

# 8. PANEL B : évolution temporelle selon le pourcentage APOBEC3 --------------

p_percentage <- ggplot(
  df_data,
  aes(
    x = date,
    y = apobec3_mutations
  )
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    color = "#D55E00",
    fill = "#D55E00",
    alpha = 0.20,
    se = TRUE
  ) +
  geom_point(
    aes(
      color = pct_apobec3,
      size = total_mutations
    ),
    alpha = 0.85
  ) +
  scale_color_viridis_c(
    option = "viridis",
    name = "APOBEC3 (%)"
  ) +
  scale_size_continuous(
    range = c(2, 7),
    name = "Total mutations"
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b %Y",
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  labs(
    title = "B. Temporal trend by APOBEC3 percentage",
    subtitle = subtitle_text,
    x = "Collection date",
    y = NULL
  ) +
  theme_commun

# 9. Combiner les panneaux et regrouper les légendes --------------------------

figure_final <- (p_location | p_percentage) +
  plot_layout(
    guides = "collect",
    widths = c(1, 1)
  ) +
  plot_annotation(
    title = paste(
      "Temporal accumulation of APOBEC3-associated mutations",
      "in MPXV, Côte d'Ivoire, 2025"
    ),
    caption = paste(
      "Panel A: point colour indicates sampling location.",
      "Panel B: point colour indicates APOBEC3 percentage.",
      "Point size indicates total mutations."
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15
      ),
      plot.caption = element_text(size = 9, color = "grey40"),
      plot.background = element_rect(
        fill = "white",
        color = NA
      )
    )
  ) &
  theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left"
  )

# 10. Afficher la figure ------------------------------------------------------

print(figure_final)

# 11. Exporter la figure en PNG et PDF ----------------------------------------

ggsave(
  filename = "Figure_Mpox_APOBEC3_Temporal_Final_2025.png",
  plot = figure_final,
  width = 16,
  height = 7.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "Figure_Mpox_APOBEC3_Temporal_Final_2025.pdf",
  plot = figure_final,
  width = 16,
  height = 7.5,
  units = "in",
  bg = "white"
)

# 12. Exporter les données analysées ------------------------------------------

write_csv(
  df_data,
  "Mpox_2025_APOBEC3_Temporal_Analysis.csv"
)

cat("\nAnalyse terminée avec succès.\n")
cat("Fichiers générés :\n")
cat("- Figure_Mpox_APOBEC3_Temporal_Final_2025.png\n")
cat("- Figure_Mpox_APOBEC3_Temporal_Final_2025.pdf\n")
cat("- Mpox_2025_APOBEC3_Temporal_Analysis.csv\n")

