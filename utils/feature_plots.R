if (!require("pacman"))
  install.packages("pacman")
pacman::p_load("Seurat", "Signac", "tidyverse", "ggplot2"
)

# functions for feature plots
create_expr_feature_plot = function(marker_gene) {
  # feature plot
  plot = FeaturePlot(
    object = coembed.scrna,
    features = marker_gene,
    min.cutoff = min(coembed.scrna@assays$RNA@data[marker_gene, ]),
    max.cutoff = max(coembed.scrna@assays$RNA@data[marker_gene, ]),
    pt.size = 4,
    raster = TRUE,
    order = TRUE
  ) +
    xlim(-15, 15) +
    ylim(-10, 10) +
    scale_color_gradient2(
      low = "#edf8b1",
      mid = "#7fcdbb",
      high = "#225ea8",
      midpoint = mean(c(
        min(coembed.scrna@assays$RNA@data[marker_gene, ]),
        max(coembed.scrna@assays$RNA@data[marker_gene, ])
      ))
    ) +
    theme(
      legend.position = 'bottom',
      text = element_text(size = 15),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoAxes()
  return(print(plot))
}

create_g4_feature_plot = function(marker_gene) {
  plot = FeaturePlot(
    object = coembed.g4,
    features = marker_gene,
    min.cutoff = min(coembed.g4@assays$GA@data[marker_gene, ]),
    max.cutoff = max(coembed.g4@assays$GA@data[marker_gene, ]),
    pt.size = 4,
    raster = TRUE,
    order = TRUE
  ) +
    xlim(-15, 15) +
    ylim(-10, 10) +
    scale_color_gradient2(
      low = "#fee5d9",
      mid = "#fcae91",
      high = "#cb181d",
      midpoint = mean(c(
        min(coembed.g4@assays$GA@data[marker_gene, ]),
        max(coembed.g4@assays$GA@data[marker_gene, ])
      ))
    ) +
    theme(
      legend.position = 'bottom',
      text = element_text(size = 15),
      plot.title = element_text(size = 20),
      axis.text.x = element_text(size = 25, color = "black"),
      axis.text.y = element_text(size = 25, color = "black")
    ) +
    NoAxes()
  return(print(plot))
}

check_featureplot_warnings = function(gene) {
  if (!gene %in% rownames(coembed.g4@assays$GA@data))
  {
    return(FALSE)
  }
  tt <- tryCatch(
    create_g4_feature_plot(gene),
    error = function(e)
      e,
    warning = function(w)
      w
  )
  tt2 <- tryCatch(
    create_g4_feature_plot(gene),
    error = function(e)
      e,
    warning = function(w)
      w
  )
  
  if (is(tt, "warning")) {
    return(FALSE)
  } else {
    return(TRUE)
  }
  
}

get_best_marker = function(cell_type) {
  print(cell_type)
  marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
    arrange(avg_log2FC) %>% top_n(n = 1, wt = avg_log2FC) %>% pull(gene)
  
  if (marker %in% rownames(coembed.g4@assays$GA@data) &
      check_featureplot_warnings(marker)) {
    return(marker)
  } else {
    marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
      arrange(avg_log2FC) %>% top_n(n = 40, wt = avg_log2FC) %>% pull(gene)
    for (i in marker) {
      if (i %in% rownames(coembed.g4@assays$GA@data) &
          check_featureplot_warnings(i)) {
        return(i)
      }
    }
  }
}

get_underexpr_marker = function(cell_type) {
  print(cell_type)
  marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
    arrange(avg_log2FC) %>% top_n(n = -1, wt = avg_log2FC) %>% pull(gene)
  
  if (marker %in% rownames(coembed.g4@assays$GA@data) &
      check_featureplot_warnings(marker)) {
    return(marker)
  } else {
    marker = markers %>% dplyr::filter(str_detect(cluster, cell_type)) %>%
      arrange(avg_log2FC) %>% top_n(n = -40, wt = avg_log2FC) %>% pull(gene)
    for (i in marker) {
      if (i %in% rownames(coembed.g4@assays$GA@data) &
          check_featureplot_warnings(i)) {
        return(i)
      }
    }
  }
}