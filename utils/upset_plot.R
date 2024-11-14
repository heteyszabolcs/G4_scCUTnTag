# packages
if (!require("pacman"))
  install.packages("pacman")
pacman::p_load(
  "tidyverse",
  "glue",
  "ggplot2",
  "data.table",
  "GenomicFeatures",
  "GenomicRanges",
  "cowplot"
)

# customized upset plot function
plot_freq_intersect <-
  function(dat,
           .by = "Group",
           .levels = NA,
           .split = "Category",
           .color = 1:10,
           top_n = 10) {
    # Args:
    #   dat: gene_id (unique)
    #   .by: row feature
    #   .levels: dot plot row names
    #   .split: bar split feature
    
    if (is.na(.levels[1])) {
      .levels <- seq_len(as.character(nchar(dat[1, .by])))
    }
    
    # limit less frequent groups
    top_groups <- names(tail(sort(table(dat[, .by])), top_n))
    dat <- dat[dat[[.by]] %in% top_groups, ]
    
    dat$Group <-
      factor(dat[, .by], levels = names(sort(table(dat[, .by]), decreasing = T)))
    dat$Type <- dat[, .split]
    
    # barplot
    dat_g1 <- dplyr::count(dat, Group, Type, sort = TRUE)
    
    frac_tbl <- table(dat[, "Group"], dat[, "Type"])
    frac_tbl <- frac_tbl / rowSums(frac_tbl)
    
    dat_g1_text <- data.frame(Group = unique(dat_g1$Group),
                              #Fraction = paste0(round(frac_tbl[, 1], 2) * 100, "%"),
                              Type = colnames(frac_tbl)[1])
    
    Group_n <- dplyr::count(dat, Group, sort = TRUE)
    dat_g1_text$n  <-
      Group_n[match(dat_g1_text$Group, Group_n$Group), "n"]
    
    g1 <- ggplot(dat_g1, aes(x = Group, y = n, fill = Type)) +
      geom_bar(stat = "identity", width = 0.6) +
      geom_text(
        data = dat_g1_text,
        aes(x = Group, y = n, label = n),
        vjust = -0.25,
        hjust = 0.5,
        size = 5
      ) +
      xlab("") + ylab("Number of regions") +
      theme_classic() +
      scale_fill_manual(name = .split, values = c("#636363", .color)) +
      theme(
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(c(0, 1, 1, 1))
      )
    
    # dot plot
    dat_tile <- NULL
    for (m in seq_along(levels(dat$Group))) {
      tmp <- cbind(x = m,
                   y = .levels,
                   color = ifelse(strsplit(levels(dat$Group)[m], "")[[1]] == 1, 1, 0))
      dat_tile <- rbind(dat_tile, tmp)
    }
    
    dat_tile <- as.data.frame(dat_tile)
    dat_tile$x <-
      factor(dat_tile$x, levels = order(as.numeric(dat_tile$x)))
    dat_tile$y <- factor(dat_tile$y, levels = rev(.levels))
    
    g2 <-
      ggplot(dat_tile, aes(
        x = x,
        y = y,
        color = factor(color),
        group = x
      )) +
      geom_point(size = 4) +
      scale_color_manual(values = c("grey90", .color)) +
      xlab("") + ylab("") +
      ggpubr::theme_pubclean() +
      theme(
        panel.grid = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(c(1, 1, 0, 1))
      )
    
    # add lines
    dat_seg <- NULL
    for (i in unique(dat_tile$x)) {
      tmp <- dat_tile[dat_tile$x == i, ]
      tmp <- tmp[tmp$color == 1, ]
      if (nrow(tmp) > 1) {
        dat_seg <- rbind(dat_seg,
                         data.frame(
                           x = tmp$x[-nrow(tmp)],
                           y = tmp$y[-nrow(tmp)],
                           xend = tmp$x[-1],
                           yend = tmp$y[-1]
                         ))
      }
    }
    
    if (!is.null(dat_seg)) {
      levels(dat_seg$x) <- levels(dat_seg$xend) <- levels(dat_tile$x)
      levels(dat_seg$y) <- levels(dat_seg$yend) <- levels(dat_tile$y)
      dat_seg$color <- 1
      
      g2 <-
        g2 + geom_segment(data = dat_seg,
                          aes(
                            x = x,
                            y = y,
                            xend = xend,
                            yend = yend
                          ),
                          lty = 2)
    }
    
    cowplot::plot_grid(
      g1,
      g2,
      ncol = 1,
      align = "v",
      rel_heights = c(1, 0.4)
    )
  }

overlap = function(peak_set1, peak_set2) {
  ol = peak_set1[queryHits(findOverlaps(
    peak_set1,
    peak_set2,
    maxgap = -1L,
    minoverlap = 1,
    type = c("any")
  )), ]
  
  ol = as_tibble(ol)
  ol = ol %>% mutate(id = paste0(seqnames, "_", start, "_", end)) %>% pull(id)
  return(ol)
}

# create ID column
add_id = function(peak_set = peak0) {
  peak_set = as_tibble(peak_set)
  ids = peak_set %>% mutate(id = paste0(seqnames, "_", start, "_", end)) %>% pull(id)
  return(ids)
}

# create input data frame for upset plot
upset_input = function(peak_set_1, peak_set_2, bulk_set) {
  peak_all =
    GenomicRanges::reduce(c(bulk_set, peak_set_1, peak_set_2))
  peak_set_t = as_tibble(peak_all)
  peak_set_t = peak_set_t %>% mutate(id = paste0(seqnames, "_", start, "_", end))
  
  # create ids for upset plot
  groups = paste0(as.numeric(
    peak_set_t$id %in% overlap(peak_set1 = peak_all, peak_set2 = peak_set_1)
  ),
  paste0(
    as.numeric(
      peak_set_t$id %in% overlap(peak_set1 = peak_all, peak_set2 = peak_set_2)
    ),
    paste0(as.numeric(
      peak_set_t$id %in% overlap(peak_set1 = peak_all, peak_set2 = bulk_set)
    ))
  ))
  
  input = tibble(id = peak_set_t$id,
                 category = "",
                 group = groups)
  input = as.data.frame(input)
  
  return(input)
}
