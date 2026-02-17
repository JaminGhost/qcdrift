#' Generate PCA Scores Plot
#'
#' @param pcs List of type `prcomp`
#' @param title_prefix Character figure title label
#' @param pc_x Integer identifying the PC to be plotted on the x-axis
#' @param pc_y Integer identifying the PC to be plotted on the y-axis
#' @param qc_starts_with Character specifying the prefix used to identify QC samples
#'
#' @return A ggplot2 object representing the PCA scores plot.
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 aes ggplot geom_point labs scale_color_manual scale_shape_manual scale_size_manual theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @importFrom tidyr as_tibble
#' @export
generate_pca_plot <- function(pcs, title_prefix="PCA",  pc_x = 1, pc_y = 2, qc_starts_with = 'QC')
{

  toplot <- as_tibble(pcs$x) |>
    mutate(sample = rownames(pcs$x),

           # Determine Sample Type (QC vs Sample) - sample names should be rownames of pcs$x
           qc = ifelse(grepl(paste0('^', qc_starts_with), sample),
                       "QC",
                       "Sample"))


  # Variance Explained
  var_expl <- round(pcs$sdev^2 / sum(pcs$sdev^2) * 100, 2)

  # Plot
  p <- ggplot(toplot, aes(x=.data[[paste0("PC", pc_x)]], y=.data[[paste0("PC", pc_y)]],
                          color=qc, shape=qc)) +
    geom_point(aes(size=qc), alpha=0.8) +

    # New Layer: Labels only for QC points
    geom_text_repel(
        data = subset(toplot, qc == "QC"),
        aes(label = sample), # Replace 'sample_id' with your actual ID column name
        size = 3,
        show.legend = FALSE
    ) +

    # point characteristics
    scale_size_manual(values=c("QC"=4, "Sample"=2.5)) +
    scale_color_manual(values=c("QC"="blue", "Sample"="salmon")) +
    scale_shape_manual(values=c("QC"=18, "Sample"=16)) +

    labs(title = paste(title_prefix, "Scores Plot"),
         x = paste0("PC", pc_x, " (", var_expl[pc_x], "%)"),
         y = paste0("PC", pc_y, " (", var_expl[pc_y], "%)"),
         color = 'Sample type', shape = 'Sample type', size = 'Sample type') +
    theme_minimal()

  return(p)
}

#' Generate CV Plot
#'
#' @param data A data frame containing metabolite abundance data.
#' @param vbl The name of the column in `data` to use for calculating the coefficient of variation.
#' @param title_prefix A character string to use as the prefix for the plot title.
#'
#' @return A ggplot2 object representing the CV plot.
#'
#' @importFrom dplyr filter group_by summarise arrange desc
#' @importFrom ggplot2 ggplot aes geom_point geom_hline labs theme_minimal theme element_blank
#' @export
generate_cv_plot <- function(data, vbl, title_prefix="Metabolite CV (QC Samples)")
{
  # Use the boolean 'qc' column from the data directly
  # Group by metabolite and calculate CV on the variable 'vbl'
  cv_df <- data |>
    filter(qc == TRUE) |>
    group_by(metabolites) |>
    summarise(
      Mean = mean({{vbl}}, na.rm=TRUE),
      SD = sd({{vbl}}, na.rm=TRUE),
      CV = (SD/Mean) * 100,
      .groups = "drop"
    ) |>
    filter(is.finite(CV)) |>
    arrange(desc(CV))

  cv_df$metabolites <- factor(cv_df$metabolites, levels=cv_df$metabolites)

  p <- ggplot(cv_df, aes(x=metabolites, y=CV)) +
    geom_point(color="steelblue", alpha=0.6, size=1) +
    geom_hline(yintercept = 20, linetype="dashed", color="red") +
    labs(title=title_prefix, y="Coefficient of Variation (%)", x="Metabolites") +
    theme_minimal() +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank())

  return(p)
}

#' Generate Violin Plot
#'
#' @param data A data frame containing metabolite abundance data.
#' @param vbl The name of the column in `data` to use for plotting.
#' @param title_prefix A character string to use as the prefix for the plot title.
#'
#' @return A ggplot2 object representing the violin plot.
#'
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot labs theme_minimal theme element_text
#' @export
generate_qc_violin <- function(data, vbl, title_prefix="QC Abundance Distribution")
{
  # Filter for QCs and transform the variable for plotting
  plot_data <- data |>
    filter(qc == TRUE) |>
    mutate(log_abund = log10({{vbl}} + 1))

  p <- ggplot(plot_data, aes(x=sample, y=log_abund, fill=sample)) +
    geom_violin(trim=FALSE, alpha=0.6) +
    geom_boxplot(width=0.1, fill="white", outlier.shape=NA) +
    labs(title=title_prefix, y="Log10(Abundance)") +
    theme_minimal() +
    theme(legend.position="none", axis.text.x = element_text(angle=45, hjust=1))

  return(p)
}

#' Generate Heatmap
#'
#' @param data A data frame containing metabolite abundance data.
#' @param vbl The name of the column in `data` to use for the heatmap.
#' @param title_prefix A character string to use as the prefix for the plot title.
#'
#' @return A ggplot2 object representing the heatmap.
#'
#' @importFrom dplyr filter group_by mutate ungroup
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 labs theme_minimal theme element_blank element_text
#' @export
generate_qc_heatmap <- function(data, vbl, title_prefix="QC Heatmap")
{
  # Filter QCs -> Scale per metabolite -> Plot
  # This uses group_by to scale properly without pivoting!
  plot_data <- data |>
    filter(qc == TRUE) |>
    group_by(metabolites) |>
    mutate(ZScore = as.numeric(scale({{vbl}}))) |>
    ungroup()

  p <- ggplot(plot_data, aes(x = sample, y = metabolites, fill = ZScore)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = title_prefix, x = "QC Sample", y = "Metabolite") +
    theme_minimal() +
    theme(axis.text.y = element_blank(), axis.text.x = element_text(angle=90))

  return(p)
}


#' Generate a pairs plot from a PCA data object
#' 
#' @param pcs List of type `prcomp`
#' @param npcs Number of principal components to plot
#' @param qc_starts_with Character specifying the prefix used to identify QC samples
#' 
#' @return A ggplot2 object with the pairs plot
#' 
#' @importFrom dplyr mutate
#' @importFrom GGally ggpairs wrap
#' @importFrom ggplot2 aes scale_color_manual scale_fill_manual scale_shape_manual scale_size_manual theme_minimal
#' @importFrom tidyr as_tibble
#' 
#' @export
generate_pca_pairs <- function(pcs, npcs = 4, qc_starts_with = 'QC')
{
  var_expl <- round(pcs$sdev^2 / sum(pcs$sdev^2) * 100, 2)
  toplot <- as_tibble(pcs$x) |>
    mutate(sample = rownames(pcs$x),
           qc = ifelse(grepl(paste0('^', qc_starts_with), sample), "QC", "Sample"))

  p <- ggpairs(toplot, columns = 1:npcs,
               aes(color = qc, shape = qc, fill = qc, size = qc, alpha = 0.8),
               diag = list(continuous = wrap("densityDiag", linewidth = 0.5)),
               columnLabels = paste0(colnames(pcs$x)[1:npcs], " (", var_expl[1:npcs], "%)")) +
    scale_size_manual(values=c("QC"=3, "Sample"=2)) +
    scale_color_manual(values=c("QC"="blue", "Sample"="salmon")) +
    scale_fill_manual(values=c("QC"="blue", "Sample"="salmon")) +
    scale_shape_manual(values=c("QC"=18, "Sample"=16)) +
    theme_minimal()

  return(p)
}
