# ============================================================
# FINAL MANUSCRIPT FIGURES
# Six-Dataset Biomedical Classification Benchmark
# ============================================================

# Run this script from the root of the repository:
# biomedical-kpca-benchmark/

required_packages <- c(
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

library(ggplot2)
library(grid)

PROJECT_DIR <- getwd()

RESULT_DIR <- file.path(
  PROJECT_DIR,
  "results",
  "benchmark_results_final"
)

final_results <- readRDS(
  file.path(
    RESULT_DIR,
    "complete_benchmark_results_6datasets_FINAL.rds"
  )
)

outer_results   <- final_results$outer_results
summary_results <- final_results$summary_results

stopifnot(
  nrow(outer_results) == 720,
  nrow(summary_results) == 48
)
# ============================================================
# FIGURE 1
# Outer-test ROC-AUC across methods
# ============================================================

required_cols <- c(
  "Dataset",
  "Representation",
  "Classifier",
  "ROC_AUC"
)

stopifnot(all(required_cols %in% names(outer_results)))

# Dataset order
outer_results$Dataset <- factor(
  outer_results$Dataset,
  levels = c(
    "Colon",
    "Gravier",
    "Leukemia",
    "Pima",
    "Prostate",
    "Wisconsin"
  )
)

# Representation labels
outer_results$Representation_Label <- factor(
  outer_results$Representation,
  levels = c(
    "raw",
    "pca",
    "rbf_kpca",
    "local_kpca"
  ),
  labels = c(
    "Raw",
    "PCA",
    "RBF-KPCA",
    "Local-KPCA"
  )
)

# Classifier labels
outer_results$Classifier_Label <- factor(
  outer_results$Classifier,
  levels = c(
    "linear",
    "radial"
  ),
  labels = c(
    "Linear",
    "Radial"
  )
)

# Method labels
outer_results$Method <- paste(
  outer_results$Representation_Label,
  outer_results$Classifier_Label,
  sep = " + "
)

method_order <- c(
  "Raw + Linear",
  "Raw + Radial",
  "PCA + Linear",
  "PCA + Radial",
  "RBF-KPCA + Linear",
  "RBF-KPCA + Radial",
  "Local-KPCA + Linear",
  "Local-KPCA + Radial"
)

outer_results$Method <- factor(
  outer_results$Method,
  levels = rev(method_order)
)

# Representation colors
representation_colors <- c(
  "Raw"        = "#00BFC4",
  "PCA"        = "#7CAE00",
  "RBF-KPCA"   = "#C77CFF",
  "Local-KPCA" = "#F8766D"
)

# Mean ROC-AUC
mean_results <- aggregate(
  ROC_AUC ~
    Dataset +
    Method +
    Classifier_Label +
    Representation_Label,
  data = outer_results,
  FUN = mean
)

# Create Figure 1
ggfig1 <- ggplot(
  outer_results,
  aes(
    x = ROC_AUC,
    y = Method,
    color = Representation_Label,
    shape = Classifier_Label
  )
) +
  
  geom_point(
    alpha = 0.45,
    size = 2.0,
    position = position_jitter(
      width = 0,
      height = 0.045
    )
  ) +
  
  geom_point(
    data = mean_results,
    aes(
      x = ROC_AUC,
      y = Method,
      shape = Classifier_Label
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 3.6
  ) +
  
  facet_wrap(
    ~ Dataset,
    ncol = 2,
    scales = "free_x"
  ) +
  
  scale_color_manual(
    values = representation_colors,
    breaks = c(
      "Raw",
      "PCA",
      "RBF-KPCA",
      "Local-KPCA"
    ),
    labels = c(
      "Raw",
      "PCA",
      "RBF-KPCA",
      "Locally scaled KPCA"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Linear" = 16,
      "Radial" = 17
    )
  ) +
  
  labs(
    x = "Outer-test ROC-AUC",
    y = NULL,
    color = "Representation",
    shape = "Classifier"
  ) +
  
  guides(
    shape = guide_legend(
      order = 1,
      nrow = 1
    ),
    color = guide_legend(
      order = 2,
      nrow = 1
    )
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    strip.background = element_rect(
      fill = "grey88",
      color = "grey40",
      linewidth = 0.4
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 6)
    ),
    
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    
    panel.grid.major = element_line(
      linewidth = 0.35,
      color = "grey88"
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey35",
      fill = NA,
      linewidth = 0.5
    ),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    
    legend.title = element_text(
      size = 10,
      face = "plain"
    ),
    
    legend.text = element_text(size = 9),
    
    legend.key.width = grid::unit(
      0.55,
      "cm"
    ),
    
    legend.spacing.x = grid::unit(
      0.10,
      "cm"
    ),
    
    legend.margin = margin(
      t = 4,
      r = 0,
      b = 0,
      l = 0
    ),
    
    panel.spacing = grid::unit(
      0.25,
      "cm"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 6,
      b = 5,
      l = 5
    )
  )

# Display Figure 1
ggfig1

# Save Figure 1
ggsave(
  file.path(
    RESULT_DIR,
    "Figure1_Outer_Test_ROC_AUC_6Datasets.pdf"
  ),
  plot = ggfig1,
  width = 11,
  height = 7.4,
  units = "in",
  limitsize = FALSE
)

ggsave(
  file.path(
    RESULT_DIR,
    "Figure1_Outer_Test_ROC_AUC_6Datasets.png"
  ),
  plot = ggfig1,
  width = 11,
  height = 7.4,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)
# ============================================================
# FIGURE 2
# Fold-matched change in ROC-AUC relative to raw predictors
# ============================================================

raw_auc <- outer_results[
  outer_results$Representation == "raw",
  c(
    "Dataset",
    "Outer_Repeat",
    "Outer_Fold",
    "Classifier",
    "ROC_AUC"
  )
]

names(raw_auc)[
  names(raw_auc) == "ROC_AUC"
] <- "ROC_AUC_Raw"

transformed_auc <- outer_results[
  outer_results$Representation %in%
    c(
      "pca",
      "rbf_kpca",
      "local_kpca"
    ),
  c(
    "Dataset",
    "Outer_Repeat",
    "Outer_Fold",
    "Representation",
    "Classifier",
    "ROC_AUC"
  )
]

names(transformed_auc)[
  names(transformed_auc) == "ROC_AUC"
] <- "ROC_AUC_Transformed"

matched_auc <- merge(
  transformed_auc,
  raw_auc,
  by = c(
    "Dataset",
    "Outer_Repeat",
    "Outer_Fold",
    "Classifier"
  )
)

matched_auc$Delta_ROC_AUC <-
  matched_auc$ROC_AUC_Transformed -
  matched_auc$ROC_AUC_Raw

representation_names <- c(
  "pca"        = "PCA",
  "rbf_kpca"   = "RBF-KPCA",
  "local_kpca" = "Local-KPCA"
)

classifier_names <- c(
  "linear" = "Linear",
  "radial" = "Radial"
)

matched_auc$Contrast <- paste0(
  representation_names[
    matched_auc$Representation
  ],
  " vs Raw: ",
  classifier_names[
    matched_auc$Classifier
  ]
)

contrast_order <- c(
  "PCA vs Raw: Linear",
  "PCA vs Raw: Radial",
  "RBF-KPCA vs Raw: Linear",
  "RBF-KPCA vs Raw: Radial",
  "Local-KPCA vs Raw: Linear",
  "Local-KPCA vs Raw: Radial"
)

matched_auc$Contrast <- factor(
  matched_auc$Contrast,
  levels = rev(contrast_order)
)

matched_auc$Dataset <- factor(
  matched_auc$Dataset,
  levels = c(
    "Colon",
    "Gravier",
    "Leukemia",
    "Pima",
    "Prostate",
    "Wisconsin"
  )
)

stopifnot(nrow(matched_auc) == 540)

mean_delta <- aggregate(
  Delta_ROC_AUC ~ Dataset + Contrast,
  data = matched_auc,
  FUN = mean
)

ggfig2 <- ggplot(
  matched_auc,
  aes(
    x = Delta_ROC_AUC,
    y = Contrast
  )
) +
  
  geom_vline(
    xintercept = 0,
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_point(
    color = "grey45",
    alpha = 0.55,
    size = 2.0,
    position = position_jitter(
      width = 0,
      height = 0.045
    )
  ) +
  
  geom_point(
    data = mean_delta,
    aes(
      x = Delta_ROC_AUC,
      y = Contrast
    ),
    inherit.aes = FALSE,
    shape = 18,
    color = "black",
    size = 3.7
  ) +
  
  facet_wrap(
    ~ Dataset,
    ncol = 2,
    scales = "free_x"
  ) +
  
  labs(
    title = "Change in ROC-AUC relative to raw predictors",
    x = expression(Delta * "ROC-AUC"),
    y = NULL
  ) +
  
  theme_bw(
    base_size = 11
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0,
      margin = margin(
        b = 7
      )
    ),
    
    strip.background = element_rect(
      fill = "grey82",
      color = "grey40",
      linewidth = 0.5
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.title.x = element_text(
      size = 11,
      margin = margin(
        t = 6
      )
    ),
    
    axis.text.x = element_text(
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey40",
      fill = NA,
      linewidth = 0.5
    ),
    
    panel.spacing = grid::unit(
      0.25,
      "cm"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 6,
      b = 5,
      l = 5
    )
  )

ggfig2

ggsave(
  file.path(
    RESULT_DIR,
    "Figure2_Matched_Delta_ROC_AUC_6Datasets.pdf"
  ),
  plot = ggfig2,
  width = 11,
  height = 7.8,
  units = "in",
  limitsize = FALSE
)

ggsave(
  file.path(
    RESULT_DIR,
    "Figure2_Matched_Delta_ROC_AUC_6Datasets.png"
  ),
  plot = ggfig2,
  width = 11,
  height = 7.8,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)
# ============================================================
# FIGURE 3
# Predictive performance versus computational cost
# ============================================================

required_summary_cols <- c(
  "Dataset",
  "Representation",
  "Classifier",
  "ROC_AUC_Mean",
  "Mean_Elapsed_Seconds"
)

stopifnot(
  all(required_summary_cols %in% names(summary_results))
)

summary_results$Dataset <- factor(
  summary_results$Dataset,
  levels = c(
    "Colon",
    "Gravier",
    "Leukemia",
    "Pima",
    "Prostate",
    "Wisconsin"
  )
)

summary_results$Representation_Label <- factor(
  summary_results$Representation,
  levels = c(
    "local_kpca",
    "pca",
    "raw",
    "rbf_kpca"
  ),
  labels = c(
    "Local-KPCA",
    "PCA",
    "Raw",
    "RBF-KPCA"
  )
)

summary_results$Classifier_Label <- factor(
  summary_results$Classifier,
  levels = c(
    "linear",
    "radial"
  ),
  labels = c(
    "Linear",
    "Radial"
  )
)

ggfig3 <- ggplot(
  summary_results,
  aes(
    x = Mean_Elapsed_Seconds,
    y = ROC_AUC_Mean,
    color = Representation_Label,
    shape = Classifier_Label
  )
) +
  
  geom_point(
    size = 4.6,
    alpha = 0.95
  ) +
  
  facet_wrap(
    ~ Dataset,
    ncol = 2,
    scales = "free_y"
  ) +
  
  scale_x_log10() +
  
  scale_color_manual(
    values = representation_colors,
    breaks = c(
      "Local-KPCA",
      "PCA",
      "Raw",
      "RBF-KPCA"
    )
  ) +
  
  scale_shape_manual(
    values = c(
      "Linear" = 16,
      "Radial" = 17
    )
  ) +
  
  labs(
    title = "Predictive performance versus computational cost",
    x = "Mean computation time per outer evaluation (seconds, log scale)",
    y = "Mean ROC-AUC",
    color = "Representation",
    shape = "Classifier"
  ) +
  
  guides(
    shape = guide_legend(
      order = 1,
      nrow = 1
    ),
    color = guide_legend(
      order = 2,
      nrow = 1
    )
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 13,
      hjust = 0,
      margin = margin(b = 8)
    ),
    
    strip.background = element_rect(
      fill = "grey82",
      color = "grey40",
      linewidth = 0.5
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.title.x = element_text(
      size = 11,
      margin = margin(t = 6)
    ),
    
    axis.title.y = element_text(
      size = 11,
      margin = margin(r = 6)
    ),
    
    axis.text.x = element_text(
      size = 9
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    panel.grid.major = element_line(
      color = "grey88",
      linewidth = 0.35
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey40",
      fill = NA,
      linewidth = 0.5
    ),
    
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "center",
    
    legend.title = element_text(
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    legend.key.width = grid::unit(
      0.55,
      "cm"
    ),
    
    legend.spacing.x = grid::unit(
      0.10,
      "cm"
    ),
    
    panel.spacing = grid::unit(
      0.25,
      "cm"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 6,
      b = 5,
      l = 5
    )
  )

ggfig3

ggsave(
  file.path(
    RESULT_DIR,
    "Figure3_Performance_vs_Cost_6Datasets.pdf"
  ),
  plot = ggfig3,
  width = 11,
  height = 7.8,
  units = "in",
  limitsize = FALSE
)

ggsave(
  file.path(
    RESULT_DIR,
    "Figure3_Performance_vs_Cost_6Datasets.png"
  ),
  plot = ggfig3,
  width = 11,
  height = 7.8,
  units = "in",
  dpi = 600,
  limitsize = FALSE
)
