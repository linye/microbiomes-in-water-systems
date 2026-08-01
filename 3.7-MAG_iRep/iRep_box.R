# Load required packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggpubr)
library(rstatix)
install.packages("patchwork")  # Install if not already available

# Read data
data <- read.table("iRep_box.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Inspect data
head(data)
str(data)

# Data cleaning: handle missing values
data <- data %>%
  mutate(across(c(Order, Family, Genus), ~ifelse(. == "" | . == "o__" | . == "f__" | . == "g__", "Unclassified", .)))

# 1. Plot iRep distribution per group (histogram + boxplot) - corrected version
create_group_plots <- function(data) {
  # Ensure all groups are present; if missing, create empty data
  all_groups <- c("DW", "AS", "WW", "NW", "GW")
  
  plots <- list()
  for (group in all_groups) {
    if (group %in% data$Group) {
      group_data <- data %>% filter(Group == group)
    } else {
      group_data <- data[0, ] # empty data frame
    }
    
    # Histogram
    hist_plot <- ggplot(group_data, aes(x = iRep)) +
      geom_histogram(aes(y = after_stat(density)), bins = 20, fill = "lightblue", color = "black", alpha = 0.7) +
      geom_density(alpha = 0.2, fill = "red") +
      labs(title = paste("Group:", group),
           x = "iRep", y = "Density") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5),
            axis.title.x = element_blank())  # Remove x-axis title to let boxplot show it
    
    # Horizontal boxplot sharing the x-axis with the histogram
    box_plot <- ggplot(group_data, aes(x = iRep, y = "")) +
      geom_boxplot(fill = "lightgreen", alpha = 0.7, width = 0.2) +
      stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
      labs(x = "iRep", y = "") +
      theme_minimal() +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank())
    
    # Combine plots, ensuring x-axes align
    combined_plot <- hist_plot / box_plot + 
      plot_layout(heights = c(3, 1)) &
      scale_x_continuous(limits = range(data$iRep, na.rm = TRUE))  # Ensure same x-axis range
    
    plots[[group]] <- combined_plot
  }
  
  # Combine all group plots into one figure
  wrap_plots(plots, ncol = 2) +
    plot_annotation(title = "iRep Distribution by Group")
}

# Generate group plots
group_plot <- create_group_plots(data)
print(group_plot)

# Save group plots
ggsave("iRep_by_group.png", group_plot, width = 16, height = 12, dpi = 300)
ggsave("iRep_by_group.pdf", group_plot, width = 16, height = 12, dpi = 300)

# 2. Create boxplots by taxonomic level with significance testing
create_taxonomy_boxplots <- function(data, tax_level) {
  # Select taxonomic units with at least 3 samples
  tax_counts <- data %>%
    filter(!.data[[tax_level]] %in% c("Unclassified", "o__", "f__", "g__", "")) %>%
    group_by(.data[[tax_level]]) %>%
    filter(n() >= 3) %>%
    ungroup()
  
  if (nrow(tax_counts) == 0) {
    message(paste("No taxonomy units with sufficient data at level:", tax_level))
    return(NULL)
  }
  
  # Compute summary statistics per taxon
  tax_summary <- tax_counts %>%
    group_by(.data[[tax_level]]) %>%
    summarise(
      n = n(),
      mean_iRep = mean(iRep),
      .groups = 'drop'
    ) %>%
    arrange(desc(n))
  
  # Select top 15 most abundant taxa for visualization
  top_taxa <- head(tax_summary[[tax_level]], 15)
  plot_data <- tax_counts %>% 
    filter(.data[[tax_level]] %in% top_taxa)
  
  # Statistical tests: compare AS and WW groups against others
  significance_results <- data.frame()
  
  for (taxon in unique(plot_data[[tax_level]])) {
    taxon_data <- plot_data %>% filter(.data[[tax_level]] == taxon)
    
    if (length(unique(taxon_data$Group)) >= 2) {
      # Check if AS group exists and has data
      as_data <- taxon_data %>% filter(Group == "AS")
      other_data <- taxon_data %>% filter(Group != "AS")
      
      if (nrow(as_data) >= 2 && nrow(other_data) >= 2) {
        tryCatch({
          test_as <- wilcox.test(as_data$iRep, other_data$iRep, alternative = "greater")
          if (test_as$p.value < 0.05) {
            significance_results <- rbind(significance_results, 
                                          data.frame(taxon = taxon, 
                                                     group = "AS", 
                                                     p_value = test_as$p.value))
          }
        }, error = function(e) NULL)
      }
      
      # Check if WW group exists and has data
      ww_data <- taxon_data %>% filter(Group == "WW")
      other_data_ww <- taxon_data %>% filter(Group != "WW")
      
      if (nrow(ww_data) >= 2 && nrow(other_data_ww) >= 2) {
        tryCatch({
          test_ww <- wilcox.test(ww_data$iRep, other_data_ww$iRep, alternative = "greater")
          if (test_ww$p.value < 0.05) {
            significance_results <- rbind(significance_results, 
                                          data.frame(taxon = taxon, 
                                                     group = "WW", 
                                                     p_value = test_ww$p.value))
          }
        }, error = function(e) NULL)
      }
    }
  }
  
  # Create boxplot
  p <- ggplot(plot_data, aes(x = .data[[tax_level]], y = iRep, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), 
               alpha = 0.6, size = 1) +
    stat_compare_means(aes(group = Group), 
                       method = "wilcox.test", 
                       label = "p.signif",
                       hide.ns = TRUE) +
    labs(title = paste("iRep by", tax_level),
         x = tax_level,
         y = "iRep") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  # Add significance annotations
  if (nrow(significance_results) > 0) {
    message(paste("\nSignificant results for", tax_level, ":"))
    print(significance_results)
    
    # Add text annotations on the plot
    for (i in 1:nrow(significance_results)) {
      taxon <- significance_results$taxon[i]
      group <- significance_results$group[i]
      p_value <- significance_results$p_value[i]
      
      # Find x-axis position of this taxon
      taxon_levels <- levels(factor(plot_data[[tax_level]]))
      x_pos <- which(taxon_levels == taxon)
      
      if (length(x_pos) > 0) {
        p <- p + 
          annotate("text", 
                   x = x_pos, 
                   y = max(plot_data$iRep, na.rm = TRUE) * 1.05,
                   label = paste(group, "*"), 
                   color = "red", 
                   fontface = "bold",
                   size = 4)
      }
    }
  }
  
  return(p)
}

# Create plots for each taxonomic level
order_plot <- create_taxonomy_boxplots(data, "Order")
family_plot <- create_taxonomy_boxplots(data, "Family") 
genus_plot <- create_taxonomy_boxplots(data, "Genus")

# Display and save taxonomic plots
if (!is.null(order_plot)) {
  print(order_plot)
  ggsave("iRep_by_Order.png", order_plot, width = 14, height = 8, dpi = 300)
  ggsave("iRep_by_Order.pdf", order_plot, width = 14, height = 8, dpi = 300)
}

if (!is.null(family_plot)) {
  print(family_plot)
  ggsave("iRep_by_Family.png", family_plot, width = 14, height = 8, dpi = 300)
  ggsave("iRep_by_Family.pdf", family_plot, width = 14, height = 8, dpi = 300)
}

if (!is.null(genus_plot)) {
  print(genus_plot)
  ggsave("iRep_by_Genus.png", genus_plot, width = 14, height = 8, dpi = 300)
  ggsave("iRep_by_Genus.pdf", genus_plot, width = 14, height = 8, dpi = 300)
}

# 3. Generate summary statistics
generate_summary_stats <- function(data) {
  # Statistics by group
  group_stats <- data %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      mean_iRep = mean(iRep, na.rm = TRUE),
      sd_iRep = sd(iRep, na.rm = TRUE),
      median_iRep = median(iRep, na.rm = TRUE),
      min_iRep = min(iRep, na.rm = TRUE),
      max_iRep = max(iRep, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Overall statistics
  overall_stats <- data %>%
    summarise(
      Group = "Overall",
      n = n(),
      mean_iRep = mean(iRep, na.rm = TRUE),
      sd_iRep = sd(iRep, na.rm = TRUE),
      median_iRep = median(iRep, na.rm = TRUE),
      min_iRep = min(iRep, na.rm = TRUE),
      max_iRep = max(iRep, na.rm = TRUE)
    )
  
  stats_combined <- bind_rows(group_stats, overall_stats)
  
  # Save results
  write.csv(stats_combined, "iRep_summary_statistics.csv", row.names = FALSE)
  
  return(stats_combined)
}

# Generate and display summary statistics
summary_stats <- generate_summary_stats(data)
print(summary_stats)

# 4. Overall comparison between groups
group_comparison <- function(data) {
  # Check if there are at least two groups
  if (length(unique(data$Group)) < 2) {
    message("Not enough groups for comparison")
    return(NULL)
  }
  
  # Kruskal-Wallis test
  kruskal_test <- kruskal.test(iRep ~ Group, data = data)
  message("Kruskal-Wallis test for group differences:")
  print(kruskal_test)
  
  # If overall significant, perform pairwise comparisons
  if (kruskal_test$p.value < 0.05) {
    pairwise_test <- pairwise.wilcox.test(data$iRep, data$Group, 
                                          p.adjust.method = "BH")
    message("\nPairwise Wilcoxon tests:")
    print(pairwise_test)
  }
  
  # Plot group comparison
  p <- ggplot(data, aes(x = Group, y = iRep, fill = Group)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.6) +
    stat_compare_means(method = "kruskal.test", 
                       label.y = max(data$iRep, na.rm = TRUE) * 1.1) +
    labs(title = "iRep Comparison Across Groups",
         x = "Group", 
         y = "iRep") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

# Perform group comparison
comparison_plot <- group_comparison(data)
if (!is.null(comparison_plot)) {
  print(comparison_plot)
  ggsave("group_comparison.png", comparison_plot, width = 10, height = 8, dpi = 300)
  ggsave("group_comparison.pdf", comparison_plot, width = 10, height = 8, dpi = 300)
}

cat("Analysis completed! Check the generated PNG files and CSV statistics file.\n")