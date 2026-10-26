# =====================================================================
# Arbre phylogénétique MPXV 2025 en Côte d'Ivoire  (version 3)
#  - colonnes Clade / Country / Host RESSERRÉES
#  - Country : tous les pays, une couleur par pays
#  - Host    : origine de l'espèce hôte
# Entrées : Mpox_Tree_file.treefile  +  annotations.txt
# Sortie  : tree_mpxv.png
# =====================================================================

# install.packages(c("ape", "phytools", "ggplot2", "dplyr", "ggnewscale"))
# if (!requireNamespace("BiocManager")) install.packages("BiocManager")
# BiocManager::install("ggtree")

library(ape)
library(phytools)
library(ggtree)
library(ggplot2)
library(dplyr)
library(ggnewscale)
setwd("D:/Seq Mpox Projet rim/files/")

# ---- 1. PARAMÈTRES MODIFIABLES --------------------------------------
tree_file  <- "Mpox_Tree_file.treefile"
annot_file <- "annotations.txt"
out_file   <- "tree_mpxv.png"

titre         <- "Phylogenetic tree of MPXV sequences 2025 in CIV"
boot_min      <- 70
taille_tiplab <- 3.5
taille_point  <- 3

# ---- Positions des colonnes (multiples de la profondeur de l'arbre) --
# Plus ces valeurs sont petites, plus les colonnes sont serrées.
# Il faut garder : pos_labels < pos_clade < pos_pays < pos_hote
pos_labels <- 1.35   # début des noms de séquences
pos_clade  <- 2.35   # colonne Clade
pos_pays   <- 3.00   # colonne Country
pos_hote   <- 3.45   # colonne Host
larg_col   <- 0.07   # largeur de chaque colonne
xlim_max   <- 3.90   # limite droite de la figure

couleurs_clade <- c(
  "Ia"        = "#E41A1C",
  "Ib"        = "#377EB8",
  "IIa"       = "#4DAF4A",
  "IIb A"     = "#984EA3",
  "IIb A.1"   = "#F781BF",
  "IIb A.2.2" = "#A65628"
)

couleurs_pays <- c(
  "Cote d'Ivoire" = "#E69F00",
  "Congo"         = "#D55E00",
  "Netherlands"   = "#0072B2",
  "USA"           = "#56B4E9",
  "Ghana"         = "#009E73",
  "Liberia"       = "#CC79A7",
  "Guinea"        = "#F0E442"
)

couleurs_hote <- c(
  "Human"    = "#1B9E77",
  "monkey"   = "#7570B3",
  "squirrel" = "#E7298A",
  "Other"    = "#999999"
)

# Clades à fixer à la main (lignes "a_verifier" + écarts avec la figure)
clade_manuel <- c(
  "hMpxV_COG_LNSP-CG062_2025"      = "Ia",
  "hMpxV_COG_LNSP-CG064_2025"      = "Ia",
  "hMpxV_USA_IL-PHNPO-0100_2025"   = "Ib",
  "hMpxV_NLD_NH-AUMC-0006_2025"    = "Ib",
  "hMpxV_GHA_NPHRL-BIW-22-24_2022" = "IIa",
  "hMpxV_LBR_NPHIL-007_2024"       = "IIa",
  "hMpxV_LBR_NPHIL-037_2024"       = "IIa",
  "hMpxV_GIN_CRV-LFHVG-G1027_2024" = "IIa",
  "MpxV_CIV_HIOH-03_2022"          = "IIa",
  "MpxV_CIV_HIOH-12_2023"          = "IIa",
  "MpxV_CIV_HIOH-31_2023"          = "IIa",
  "hMpxV_CIV_MP338_2025"           = "IIb A.1",
  # écarts constatés entre la figure d'origine et le fichier d'annotations
  "hMpxV_CIV_MP285_2025-06-20"     = "IIb A.1",
  "hMpxV_CIV_MP317_2025-08-05"     = "IIb A.1",
  "hMpxV_CIV_MP339_2025-08-20"     = "IIb A"
)

# Ordre des séquences (de haut en bas) identique à la figure d'origine.
# Mettre ordre_tips <- NULL pour utiliser simplement ladderize().
ordre_tips <- c(
  paste0("hMpxV_CIV_", c(
    "MP276_2025-06-20","MP363_2025-08-27","MP224_2025-04-25","MP279_2025-06-20",
    "MP285_2025-06-20","MP219_2025-04-29","MP235_2025-05-09","MP260_2025-06-03",
    "MP284_2025-06-21","MP249_2025-05-22","MP278_2025-06-20","MP180_2025-03-23",
    "MP338_2025","MP341_2025-08-10","MP317_2025-08-05","MP339_2025-08-20",
    "MP357_2025-08-26","MP356_2025-08-26","MP321_2025-07-11","MP345_2025-08-12",
    "MP355_2025-08-23","MP354_2025-08-23","MP286_2025-06-18","MP371_2025-09-08",
    "MP350_2025-08-22","MP320_2025-07-11","MP302_2025-06-30")),
  "MpxV_CIV_HIOH-31_2023","MpxV_CIV_HIOH-03_2022","MpxV_CIV_HIOH-12_2023",
  "hMpxV_LBR_NPHIL-037_2024","hMpxV_LBR_NPHIL-007_2024",
  "hMpxV_CIV_NBM-003_1971","hMpxV_GIN_CRV-LFHVG-G1027_2024",
  "hMpxV_CIV_MP280_2025-06-23","hMpxV_GHA_NPHRL-BIW-22-24_2022",
  "hMpxV_COG_LNSP-CG064_2025","hMpxV_COG_LNSP-CG062_2025",
  "hMpxV_USA_IL-PHNPO-0100_2025","hMpxV_NLD_NH-AUMC-0006_2025"
)

# ---- 2. Lecture des données -----------------------------------------
tree <- read.tree(tree_file)
ann  <- read.delim(annot_file, sep = "\t", stringsAsFactors = FALSE,
                   check.names = FALSE, quote = "")

tree <- midpoint.root(tree)

# Supports : gère "100" ou "98.5/100" (on garde la dernière valeur)
nl <- tree$node.label
nl[nl == ""] <- NA
tree$node.label <- sub("^.*/", "", nl)

# Ordre des tips
if (!is.null(ordre_tips) && all(tree$tip.label %in% ordre_tips)) {
  tree <- tryCatch(rotateConstr(tree, ordre_tips),
                   error = function(e) ladderize(tree, right = FALSE))
} else {
  tree <- ladderize(tree, right = FALSE)
}

# ---- 3. Préparation des annotations ---------------------------------
ann <- ann %>%
  mutate(
    clade_plot = ifelse(clade_full == "a_verifier", NA, clade_full),
    pays_plot  = gsub("_", " ", country),
    pays_plot  = ifelse(pays_plot == "Cote d Ivoire", "Cote d'Ivoire", pays_plot),
    hote_plot  = host
  )

idx <- match(names(clade_manuel), ann$taxon)
ann$clade_plot[idx[!is.na(idx)]] <- clade_manuel[!is.na(idx)]

manq <- setdiff(tree$tip.label, ann$taxon)
if (length(manq) > 0)
  warning("Tips absents des annotations : ", paste(manq, collapse = ", "))
if (length(setdiff(unique(ann$pays_plot), names(couleurs_pays))) > 0)
  warning("Pays sans couleur : ",
          paste(setdiff(unique(ann$pays_plot), names(couleurs_pays)), collapse = ", "))
if (length(setdiff(unique(ann$hote_plot), names(couleurs_hote))) > 0)
  warning("Hôtes sans couleur : ",
          paste(setdiff(unique(ann$hote_plot), names(couleurs_hote)), collapse = ", "))

# ---- 4. Arbre de base + positions -----------------------------------
p0   <- ggtree(tree, size = 0.6, ladderize = FALSE) %<+% ann
d    <- p0$data
xmax <- max(d$x, na.rm = TRUE)
tips <- d %>% filter(isTip) %>% arrange(y)

x_labels <- xmax * pos_labels
x_clade  <- xmax * pos_clade
x_pays   <- xmax * pos_pays
x_hote   <- xmax * pos_hote
largeur  <- xmax * larg_col

# ---- 5. Données des colonnes ----------------------------------------
runs  <- rle(as.character(tips$clade_plot))
fin   <- cumsum(runs$lengths)
debut <- fin - runs$lengths + 1
bar_clade <- data.frame(
  clade = runs$values,
  ymin  = tips$y[debut] - 0.4,
  ymax  = tips$y[fin]   + 0.4
) %>% filter(!is.na(clade))

bar_pays <- tips %>%
  transmute(pays = pays_plot, ymin = y - 0.5, ymax = y + 0.5) %>%
  filter(!is.na(pays))

bar_hote <- tips %>%
  transmute(hote = hote_plot, ymin = y - 0.5, ymax = y + 0.5) %>%
  filter(!is.na(hote))

# ---- 6. Figure -------------------------------------------------------
y_top <- max(tips$y) + 2.5

p <- p0 +
  geom_tippoint(aes(color = clade_plot), size = taille_point) +
  scale_color_manual(values = couleurs_clade, name = "Clade / Lineage",
                     breaks = names(couleurs_clade), na.value = "grey50") +
  geom_tiplab(align = TRUE, linetype = "dotted", linesize = 0.3,
              size = taille_tiplab, offset = x_labels - xmax) +
  geom_text2(aes(subset = !isTip & !is.na(label) &
                   suppressWarnings(as.numeric(label)) >= boot_min,
                 label = label),
             size = 2.3, colour = "grey30", hjust = 1.1, vjust = -0.4) +
  
  # Colonne Clade
  new_scale_fill() +
  geom_rect(data = bar_clade, inherit.aes = FALSE,
            aes(xmin = x_clade, xmax = x_clade + largeur,
                ymin = ymin, ymax = ymax, fill = clade),
            show.legend = FALSE) +
  scale_fill_manual(values = couleurs_clade) +
  geom_text(data = bar_clade, inherit.aes = FALSE,
            aes(x = x_clade + largeur * 1.5, y = (ymin + ymax) / 2, label = clade),
            hjust = 0, size = 4) +
  
  # Colonne Country
  new_scale_fill() +
  geom_rect(data = bar_pays, inherit.aes = FALSE,
            aes(xmin = x_pays, xmax = x_pays + largeur,
                ymin = ymin, ymax = ymax, fill = pays)) +
  scale_fill_manual(values = couleurs_pays, name = "Country",
                    breaks = names(couleurs_pays), drop = TRUE) +
  
  # Colonne Host
  new_scale_fill() +
  geom_rect(data = bar_hote, inherit.aes = FALSE,
            aes(xmin = x_hote, xmax = x_hote + largeur,
                ymin = ymin, ymax = ymax, fill = hote)) +
  scale_fill_manual(values = couleurs_hote, name = "Host",
                    breaks = names(couleurs_hote), drop = TRUE) +
  
  # Titres de colonnes
  annotate("text", x = x_clade, y = y_top, label = "Clade",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate("text", x = x_pays, y = y_top, label = "Country",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate("text", x = x_hote, y = y_top, label = "Host",
           hjust = 0, fontface = "bold", size = 4.5) +
  
  geom_treescale(x = 0, y = -1.5, fontsize = 4, linesize = 0.8) +
  xlim(NA, xmax * xlim_max) +
  coord_cartesian(clip = "off") +
  labs(title   = titre,
       caption = paste0("Midpoint rooted. Scale: substitutions per site.\n",
                        "Bootstrap >= ", boot_min, "%.")) +
  theme(
    legend.position = "left",
    legend.title    = element_text(face = "bold", size = 13),
    legend.text     = element_text(size = 12),
    plot.title      = element_text(face = "bold", size = 18),
    plot.caption    = element_text(hjust = 0, colour = "grey30"),
    plot.margin     = margin(10, 20, 10, 10)
  )

print(p)

# ---- 7. Export -------------------------------------------------------
ggsave(out_file, p, width = 16, height = 10, dpi = 300, bg = "white")