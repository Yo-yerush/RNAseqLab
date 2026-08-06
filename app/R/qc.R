qc_numeric_matrix <- function(x, id_col = "gene_id") {
  if (is.null(x)) return(NULL)
  if (is.data.frame(x) && id_col %in% names(x)) {
    ids <- as.character(x[[id_col]])
    x <- x[, setdiff(names(x), id_col), drop = FALSE]
  } else {
    ids <- rownames(x)
  }
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  if (!is.null(ids) && length(ids) == nrow(x)) rownames(x) <- ids
  x
}

calculate_deseq2_qc <- function(dds, result, raw_counts, coldata,
                                initial_features = nrow(raw_counts),
                                min_count = NA_real_, input_label = "Count data") {
  raw_counts <- qc_numeric_matrix(raw_counts)
  fitted_counts <- as.matrix(DESeq2::counts(dds, normalized = FALSE))
  cd <- as.data.frame(coldata, stringsAsFactors = FALSE)
  sample_ids <- colnames(raw_counts)
  if (is.null(sample_ids)) sample_ids <- rownames(cd)
  cd_ids <- if ("sample_id" %in% names(cd)) as.character(cd$sample_id) else rownames(cd)
  condition <- if ("condition" %in% names(cd)) as.character(cd$condition[match(sample_ids, cd_ids)]) else rep("", length(sample_ids))
  label <- if ("sample_label" %in% names(cd)) as.character(cd$sample_label[match(sample_ids, cd_ids)]) else sample_ids

  size_factors <- tryCatch(DESeq2::sizeFactors(dds), error = function(e) NULL)
  normalization_label <- "DESeq2 size factor"
  if (is.null(size_factors) || !length(size_factors)) {
    normalization_factors <- tryCatch(DESeq2::normalizationFactors(dds), error = function(e) NULL)
    if (!is.null(normalization_factors) && ncol(normalization_factors) == length(sample_ids)) {
      size_factors <- apply(normalization_factors, 2, function(z) {
        z <- z[is.finite(z) & z > 0]
        if (length(z)) exp(mean(log(z))) else NA_real_
      })
      normalization_label <- "Median effective normalization factor"
    } else {
      size_factors <- rep(NA_real_, length(sample_ids))
      normalization_label <- "Normalization factor unavailable"
    }
  }
  if (!is.null(names(size_factors)) && all(sample_ids %in% names(size_factors))) {
    size_factors <- size_factors[sample_ids]
  } else if (length(size_factors) != length(sample_ids)) {
    size_factors <- rep(NA_real_, length(sample_ids))
  }

  sample_metrics <- data.frame(
    sample_id = sample_ids,
    sample_label = ifelse(is.na(label) | !nzchar(label), sample_ids, label),
    condition = condition,
    library_size = colSums(raw_counts, na.rm = TRUE),
    detected_features = colSums(raw_counts > 0, na.rm = TRUE),
    zero_fraction = colMeans(raw_counts == 0, na.rm = TRUE),
    size_factor = as.numeric(size_factors),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  cooks <- tryCatch(SummarizedExperiment::assays(dds)[["cooks"]], error = function(e) NULL)
  cooks_cutoff <- NA_real_
  if (!is.null(cooks)) {
    design_matrix <- tryCatch(stats::model.matrix(DESeq2::design(dds), data = as.data.frame(SummarizedExperiment::colData(dds))), error = function(e) NULL)
    p <- if (is.null(design_matrix)) NA_integer_ else ncol(design_matrix)
    m <- ncol(dds)
    if (is.finite(p) && m > p) cooks_cutoff <- stats::qf(0.99, p, m - p)
    sample_metrics$max_cooks <- apply(cooks, 2, function(z) suppressWarnings(max(z[is.finite(z)], na.rm = TRUE)))
    sample_metrics$max_cooks[!is.finite(sample_metrics$max_cooks)] <- NA_real_
    sample_metrics$genes_above_cooks_cutoff <- if (is.finite(cooks_cutoff)) colSums(cooks > cooks_cutoff, na.rm = TRUE) else NA_integer_
  } else {
    sample_metrics$max_cooks <- NA_real_
    sample_metrics$genes_above_cooks_cutoff <- NA_integer_
  }

  row_data <- as.data.frame(SummarizedExperiment::rowData(dds))
  dispersion <- data.frame(
    gene_id = rownames(row_data),
    base_mean = rowMeans(DESeq2::counts(dds, normalized = TRUE), na.rm = TRUE),
    gene_estimate = if ("dispGeneEst" %in% names(row_data)) row_data$dispGeneEst else NA_real_,
    fitted = if ("dispFit" %in% names(row_data)) row_data$dispFit else NA_real_,
    final = DESeq2::dispersions(dds),
    stringsAsFactors = FALSE
  )

  result_df <- as.data.frame(result)
  p_col <- if ("pvalue" %in% names(result_df)) result_df$pvalue else rep(NA_real_, nrow(result_df))
  padj_col <- if ("padj" %in% names(result_df)) result_df$padj else rep(NA_real_, nrow(result_df))
  filtering <- data.frame(
    stage = c("Imported features", "Retained after minimum-count filter", "Features with a DE p-value", "Features with adjusted p-value", "Independent-filtered features"),
    features = c(initial_features, nrow(fitted_counts), sum(is.finite(p_col)), sum(is.finite(padj_col)), sum(is.finite(p_col) & !is.finite(padj_col))),
    stringsAsFactors = FALSE
  )

  condition_counts <- table(condition)
  warnings <- character()
  if (any(condition_counts < 2)) warnings <- c(warnings, "At least one condition has fewer than two samples.")
  if (any(condition_counts < 3)) warnings <- c(warnings, "Three or more biological replicates per condition are recommended.")
  lib_ratio <- sample_metrics$library_size / stats::median(sample_metrics$library_size[sample_metrics$library_size > 0], na.rm = TRUE)
  if (any(lib_ratio < 0.5 | lib_ratio > 2, na.rm = TRUE)) warnings <- c(warnings, "At least one library size is less than half or more than twice the median.")
  if (any(sample_metrics$zero_fraction > 0.9, na.rm = TRUE)) warnings <- c(warnings, "At least one sample contains zero counts for more than 90% of imported features.")

  list(
    input_label = input_label,
    normalization_label = normalization_label,
    sample_metrics = sample_metrics,
    filtering = filtering,
    dispersion = dispersion,
    cooks_cutoff = cooks_cutoff,
    warnings = unique(warnings),
    overview = list(
      samples = ncol(raw_counts), conditions = length(unique(condition[!is.na(condition) & nzchar(condition)])),
      imported_features = initial_features, retained_features = nrow(fitted_counts), min_count = min_count
    )
  )
}

qc_relationship_matrices <- function(vst_counts) {
  mat <- qc_numeric_matrix(vst_counts)
  if (is.null(mat) || ncol(mat) < 2) stop("At least two samples are required for sample-relationship plots.")
  finite_rows <- apply(mat, 1, function(z) all(is.finite(z)))
  mat <- mat[finite_rows, , drop = FALSE]
  if (!nrow(mat)) stop("No finite VST expression rows are available.")
  list(
    correlation = stats::cor(mat, use = "pairwise.complete.obs", method = "pearson"),
    distance = as.matrix(stats::dist(t(mat))),
    hclust = stats::hclust(stats::dist(t(mat)))
  )
}

qc_matrix_long <- function(mat) {
  d <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(d) <- c("sample_x", "sample_y", "value")
  d
}

make_qc_heatmap_plot <- function(mat, title, low = "#2166AC", high = "#B2182B",
                                 midpoint = NULL, plot_theme = "minimal", font_family = "serif") {
  d <- qc_matrix_long(mat)
  p <- ggplot2::ggplot(d, ggplot2::aes(sample_x, sample_y, fill = value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::coord_equal() +
    plot_theme_choice(plot_theme, base_size = 11, font_family = font_family) +
    ggplot2::theme(axis.title = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = title, fill = NULL)
  if (is.null(midpoint)) {
    p + ggplot2::scale_fill_gradient(low = "#F7FBFF", high = high)
  } else {
    p + ggplot2::scale_fill_gradient2(low = low, mid = "white", high = high, midpoint = midpoint)
  }
}

make_qc_dendrogram_plot <- function(hc, plot_theme = "minimal", font_family = "serif") {
  n <- length(hc$order)
  leaf_x <- numeric(n)
  leaf_x[hc$order] <- seq_len(n)
  node_x <- numeric(n - 1L)
  segments <- vector("list", n - 1L)
  node_position <- function(v) if (v < 0) leaf_x[-v] else node_x[v]
  node_height <- function(v) if (v < 0) 0 else hc$height[v]
  for (i in seq_len(n - 1L)) {
    children <- hc$merge[i, ]
    x1 <- node_position(children[1]); x2 <- node_position(children[2]); y <- hc$height[i]
    node_x[i] <- mean(c(x1, x2))
    segments[[i]] <- data.frame(
      x = c(x1, x2, x1), y = c(node_height(children[1]), node_height(children[2]), y),
      xend = c(x1, x2, x2), yend = c(y, y, y)
    )
  }
  seg <- do.call(rbind, segments)
  labels <- data.frame(x = seq_len(n), y = 0, label = hc$labels[hc$order])
  ggplot2::ggplot(seg) +
    ggplot2::geom_segment(ggplot2::aes(x, y, xend = xend, yend = yend), linewidth = 0.55) +
    ggplot2::geom_text(data = labels, ggplot2::aes(x, y, label = label), angle = 45, hjust = 1, vjust = 1, inherit.aes = FALSE) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.03))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.05))) +
    plot_theme_choice(plot_theme, base_size = 11, font_family = font_family) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank()) +
    ggplot2::labs(title = "Sample dendrogram", y = "Euclidean distance")
}

make_qc_dispersion_plot <- function(qc, plot_theme = "minimal", font_family = "serif") {
  d <- qc$dispersion
  ggplot2::ggplot(d, ggplot2::aes(base_mean, final)) +
    ggplot2::geom_point(color = "gray55", alpha = 0.35, size = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = fitted), color = "#B2182B", alpha = 0.65, size = 0.55, na.rm = TRUE) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
    plot_theme_choice(plot_theme, base_size = 11, font_family = font_family) +
    ggplot2::labs(title = "DESeq2 dispersion estimates", x = "Mean normalized count", y = "Dispersion")
}

make_qc_sample_bar_plot <- function(qc, value, title, y_label, color_palette = "default",
                                    plot_theme = "minimal", font_family = "serif") {
  d <- qc$sample_metrics
  d$value <- d[[value]]
  cols <- pca_palette_values(length(unique(d$condition)), color_palette)
  p <- ggplot2::ggplot(d, ggplot2::aes(stats::reorder(sample_label, value), value, fill = condition)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::coord_flip() +
    plot_theme_choice(plot_theme, base_size = 11, font_family = font_family) +
    ggplot2::labs(title = title, x = NULL, y = y_label, fill = NULL)
  if (!is.null(cols)) p <- p + ggplot2::scale_fill_manual(values = cols)
  p
}

calculate_transcript_input_qc <- function(transcript_data) {
  if (is.null(transcript_data$counts) || is.null(transcript_data$tx2gene)) return(NULL)
  counts <- as.matrix(transcript_data$counts); storage.mode(counts) <- "numeric"
  map <- as.data.frame(transcript_data$tx2gene, stringsAsFactors = FALSE)
  names(map)[1:2] <- c("transcript_id", "gene_id")
  map <- map[!is.na(map$transcript_id) & nzchar(map$transcript_id) & !is.na(map$gene_id) & nzchar(map$gene_id), , drop = FALSE]
  mapped <- intersect(rownames(counts), map$transcript_id)
  map <- map[match(mapped, map$transcript_id), , drop = FALSE]
  tx_per_gene <- as.data.frame(table(map$gene_id), stringsAsFactors = FALSE)
  names(tx_per_gene) <- c("gene_id", "transcripts")
  cd <- as.data.frame(transcript_data$coldata, stringsAsFactors = FALSE)
  sample_ids <- colnames(counts)
  cd_ids <- if ("sample_id" %in% names(cd)) as.character(cd$sample_id) else rownames(cd)
  sample_table <- data.frame(
    sample_id = sample_ids,
    condition = if ("condition" %in% names(cd)) as.character(cd$condition[match(sample_ids, cd_ids)]) else "",
    total_transcript_count = colSums(counts, na.rm = TRUE),
    detected_transcripts = colSums(counts > 0, na.rm = TRUE),
    zero_fraction = colMeans(counts == 0, na.rm = TRUE),
    non_finite_values = colSums(!is.finite(counts)),
    stringsAsFactors = FALSE
  )
  imported_transcripts <- as.integer(transcript_data$imported_transcripts %||% nrow(counts))
  mapped_transcripts <- as.integer(transcript_data$mapped_transcripts %||% length(mapped))
  non_finite_abundance <- if (is.null(transcript_data$abundance)) NA_integer_ else sum(!is.finite(as.matrix(transcript_data$abundance)))
  summary <- data.frame(
    metric = c("Imported transcripts", "Mapped transcripts", "Mapped genes", "Genes with one transcript", "Genes with multiple transcripts", "Non-finite count values", "Non-finite abundance values"),
    value = c(imported_transcripts, mapped_transcripts, length(unique(map$gene_id)), sum(tx_per_gene$transcripts == 1), sum(tx_per_gene$transcripts >= 2), sum(!is.finite(counts)), non_finite_abundance),
    stringsAsFactors = FALSE
  )
  list(summary = summary, sample_metrics = sample_table, transcripts_per_gene = tx_per_gene,
       mapping_coverage = if (imported_transcripts) mapped_transcripts / imported_transcripts else NA_real_)
}

make_transcript_qc_plot <- function(tx_qc, plot_theme = "minimal", font_family = "serif") {
  d <- tx_qc$transcripts_per_gene
  max_display <- max(2, stats::quantile(d$transcripts, 0.99, na.rm = TRUE))
  d$display <- pmin(d$transcripts, max_display)
  ggplot2::ggplot(d, ggplot2::aes(display)) +
    ggplot2::geom_histogram(binwidth = 1, boundary = 0.5, fill = "#4C78A8", color = "white") +
    plot_theme_choice(plot_theme, base_size = 11, font_family = font_family) +
    ggplot2::labs(title = "Transcripts per mapped gene", x = "Transcripts per gene", y = "Genes")
}
