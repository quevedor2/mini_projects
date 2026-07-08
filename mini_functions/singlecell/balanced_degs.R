### Set of functions to get a balanced number of DEGs at a p_adjusted value of 
#. a set amount when comparing between two groups (e.g. WT vs KO) per celltype.
# It downsamples to Y number of cells per group so that each group is equally represented
# per celltype (e.g. 100 WT and 100 KO), then does the differential testing 
# X number of times at downsample Y.  It repeats this process for all values of 
# step_Z up to the max number of cells per cluster.

# After calclulating the mean/sd per step, it fits a LOESS curve to the data and uses 
# this fit to predict the number of DEGs for a given number of cells per cluster
# which can be used for plotting and making equivalence vlaues.  It can also plot
# the saturation curve of number of DEGs across steps.

#### Example usage:
# groupids <- c('monocle3_clusters', 'manual_anno')
# DefaultAssay(seu) <- 'SCT'
# seu <- PrepSCTFindMarkers(object = seu)
# 
# curve_data <- map_dfr(groupids, function(gid) {
#   build_saturation_curve(seu, groupid = gid, cond_i = c("KO", "WT"),
#                          step = 100, n_reps = 10)
# })
# 
# curve_summary <- get_curve_summary(curve_data)
# fits <- get_fits(curve_summary)
# resl <- get_predicted_at_ref(fits, curve_summary, ref_n=500)
# curve_fitted <- get_curve_fitted(fits)
# 
# predicted_at_ref <- resl$predicted
# ref_n=500
# 
# 
# pdf("~/xfer/fits.pdf")
# plot_saturation_curves(curve_summary, curve_fitted)
# dev.off()
# 

get_markers_at_n <- function(seu, groupid, cluster_i, cond_i, n_sub, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  Idents(seu) <- groupid
  pool_1 <- Cells(seu)[Idents(seu) == cluster_i & seu$condition == cond_i[1]]
  pool_2 <- Cells(seu)[Idents(seu) == cluster_i & seu$condition == cond_i[2]]
  if (length(pool_1) < n_sub || length(pool_2) < n_sub) return(NULL)
  
  cells_keep <- c(sample(pool_1, n_sub, replace = FALSE),
                  sample(pool_2, n_sub, replace = FALSE))
  seu_balanced <- subset(seu, cells = cells_keep)
  
  tryCatch({
    FindMarkers(
      object = seu_balanced, assay = "SCT",
      ident.1 = cond_i[1], ident.2 = cond_i[2], group.by = 'condition',
      subset.ident = cluster_i, recorrect_umi = FALSE
    ) %>%
      tibble::rownames_to_column("symbol")
  }, error = function(e) NULL)
}

build_saturation_curve <- function(seu, groupid, cond_i, step = 20, n_reps = 10, 
                                   max_n = NULL, padj_cutoff = 0.05) {
  Idents(seu) <- groupid
  tab <- table(Idents(seu), seu$condition)
  clusters <- rownames(tab)[apply(tab[, cond_i, drop = FALSE], 1, function(x) all(x >= step))]
  
  map_dfr(clusters, function(cluster_i) {
    message("\nCluster ", cluster_i, ": ", match(cluster_i, clusters), " / ", length(clusters))
    cluster_max <- min(tab[cluster_i, cond_i])
    top_n <- if (!is.null(max_n)) min(cluster_max, max_n) else cluster_max
    sizes <- seq(step, top_n, by = step)
    if (length(sizes) == 0) return(NULL)
    
    map_dfr(sizes, function(n) {
      message("Step: ", n, " / ", max(sizes))
      degs_per_rep <- map_dbl(1:n_reps, function(r) {
        m <- get_markers_at_n(seu, groupid, cluster_i, cond_i, n_sub = n, seed = r)
        if (is.null(m)) return(NA)
        sum(m$p_val_adj < padj_cutoff, na.rm = TRUE)
      })
      data.frame(groupid = groupid, cluster = cluster_i, n_cells = n,
                 rep = 1:n_reps, n_degs = degs_per_rep)
    })
  })
}

fit_saturation_curve <- function(df) {
  if (nrow(df) < 3) return(NULL)
  tryCatch({
    # nls(mean_degs ~ SSasymp(n_cells, Asym, R0, lrc), data = df)
    # lm(mean_degs ~ n_cells, data = df)
    loess(mean_degs ~ n_cells, data = df, span = 0.9)
  }, error = function(e) NULL)
}

get_curve_summary <- function(curve_data){
  curve_summary <- curve_data %>%
    group_by(groupid, cluster, n_cells) %>%
    summarize(mean_degs = mean(n_degs, na.rm = TRUE),
              sd_degs = sd(n_degs, na.rm = TRUE), .groups = "drop")
  return(curve_summary)
}

get_fits <- function(curve_summary){
  curve_summary %>%
    group_by(groupid, cluster) %>%
    group_split() %>%
    set_names(map_chr(., ~ paste(unique(.x$groupid), unique(.x$cluster), sep = "__"))) %>%
    map(fit_saturation_curve)
}

get_predicted_at_ref <- function(fits, curve_summary, ref_n){
  max_n_per_cluster <- curve_summary %>%
    group_by(groupid, cluster) %>%
    summarize(max_n = max(n_cells), .groups = "drop")
  
  predicted_at_ref <- imap_dfr(fits, function(fit, key) {
    parts <- strsplit(key, "__")[[1]]
    groupid_i <- parts[1]; cluster_i <- parts[2]
    
    max_observed_n <- max_n_per_cluster$max_n[max_n_per_cluster$groupid == groupid_i &
                                                max_n_per_cluster$cluster == cluster_i]
    if (is.null(fit)) {
      return(data.frame(groupid = groupid_i, cluster = cluster_i,
                        pred_degs = NA, max_observed_n = max_observed_n,
                        reliable = FALSE))
    }
    pred <- predict(fit, newdata = data.frame(n_cells = ref_n))
    data.frame(groupid = groupid_i, cluster = cluster_i,
               pred_degs = pred, max_observed_n = max_observed_n,
               reliable = ref_n <= max_observed_n)  # FALSE = extrapolation, treat cautiously
  })
  predicted_at_ref <- predicted_at_ref %>% arrange(desc(pred_degs))
  return(list("max_n"=max_n_per_cluster,
              "predicted"=predicted_at_ref))
}

get_curve_fitted <- function(fits){
  imap_dfr(fits, function(fit, key) {
    if (is.null(fit)) return(NULL)
    parts <- strsplit(key, "__")[[1]]
    new_n <- seq(1, ref_n * 1.5, length.out = 100)
    data.frame(groupid = parts[1], cluster = parts[2], n_cells = new_n,
               pred_degs = predict(fit, newdata = data.frame(n_cells = new_n)))
  })
}



plot_saturation_curves <- function(curve_summary, curve_fitted){
  ggplot() +
    geom_point(data = curve_summary, aes(n_cells, mean_degs, color = cluster)) +
    geom_line(data = curve_fitted, aes(n_cells, pred_degs, color = cluster)) +
    geom_vline(xintercept = ref_n, linetype = "dashed") +
    facet_wrap(~groupid, scales = "free_y") +
    labs(x = "Cells per condition (downsampled)", y = "Mean # DEGs (padj<0.05)",
         title = paste0("Saturation curves, reference N = ", ref_n)) +
    theme_minimal()
}
