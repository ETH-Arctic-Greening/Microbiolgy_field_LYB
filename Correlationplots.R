#correlation plot

# Load required libraries
library(ggplot2)
library(corrplot)
library(reshape2)

# Example Dataset
# Replace 'mtcars' with your dataset
#outliers ITS qPCR removed!
data <- read.csv("metadata231124.csv")

filtered_data <- data[data$core == "yes", ]
selected_data <- filtered_data[, grepl("^m_", colnames(filtered_data))]
selected_data <- data[, grepl("^m_", colnames(data))]
# Step 1: Compute the Correlation Matrix
cor_matrix <- cor(selected_data, use = "pairwise.complete.obs") # Handles NA values
cor_matrix <- cor(filtered_data, use = "pairwise.complete.obs")
# Step 3: Compute R^2 Values
r2_matrix <- cor_matrix^2  # Square each correlation to get R^2
rmatrix <- cor_matrix
# Step 4: Visualize R^2 Matrix
# Using `corrplot` for a heatmap-style visualization
corrplot(r2_matrix, method = "color", type = "upper", tl.col = "black",
         tl.srt = 45, addCoef.col = "black", number.cex = 0.7, is.corr = FALSE)
corrplot(rmatrix, method = "color", type = "upper", tl.col = "black",
         tl.srt = 45, addCoef.col = "black", number.cex = 0.7, is.corr = FALSE)

# Step 5: Alternative: Heatmap with ggplot2
# Melt R^2 matrix for ggplot2
melted_r2 <- melt(r2_matrix)
melted_p <- melt()
# Plot heatmap
ggplot(data = melted_r2, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                       name = "R²") +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3) + # Add R² values
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_fixed() +
  labs(title = "R² Heatmap (Columns Starting with 'm_')",
       x = "",
       y = "")






# Step 3: Generate Correlation Plots
# Iterate over all pairs of columns
combinations <- combn(colnames(filtered_data), 2, simplify = FALSE)  # All pairs of columns

# Create and save plots for each pair
plots <- lapply(combinations, function(cols) {
  x_col <- cols[1]
  y_col <- cols[2]
  
  # Perform correlation test
  cor_test <- cor.test(selected_data[[x_col]], selected_data[[y_col]], use = "pairwise.complete.obs")
  
  # Extract r^2 and p-value
  r2 <- cor_test$estimate^2
  p_value <- cor_test$p.value
  
  # Create scatter plot with annotation for r^2 and p-value
  plot <- ggplot(selected_data, aes_string(x = x_col, y = y_col)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    labs(title = paste("Correlation between", x_col, "and", y_col),
         subtitle = paste("R² =", round(r2, 2), "| p =", signif(p_value, 3)),
         x = x_col,
         y = y_col) +
    theme_minimal()
  
  return(plot)
})

# Step 4: Display Plots
# Display the plots one by one
for (plot in plots) {
  print(plot)
}

# Optional: Save plots as individual files
# lapply(seq_along(plots), function(i) {
#   ggsave(filename = paste0("correlation_plot_", i, ".png"), plot = plots[[i]], width = 6, height = 4)
# })

selected_data <- data[, grepl("^m_", colnames(data))]

# Step 3: Define columns of interest
target_columns <- c("m_graminoid")
other_columns <- setdiff(colnames(selected_data), target_columns)

# Step 4: Generate Correlation Plots
plots <- lapply(target_columns, function(target_col) {
  lapply(other_columns, function(other_col) {
    # Perform correlation test
    cor_test <- cor.test(selected_data[[target_col]], selected_data[[other_col]], use = "pairwise.complete.obs")
    
    # Extract R^2 and p-value
    r2 <- cor_test$estimate^2
    p_value <- cor_test$p.value
    
    # Create scatter plot with annotation for R^2 and p-value
    plot <- ggplot(selected_data, aes_string(x = other_col, y = target_col)) +
      geom_point(color = "blue", alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
      labs(title = paste("Correlation:", target_col, "vs", other_col),
           subtitle = paste("R² =", round(r2, 2), "| p =", signif(p_value, 3)),
           x = other_col,
           y = target_col) +
      theme_minimal()
    
    return(plot)
  })
})

# Step 5: Flatten Plot List and Display
# Flatten list of plots
plots <- do.call(c, plots)

# Display the plots one by one
for (plot in plots) {
  print(plot)
}

# Optional: Save plots as individual files
lapply(seq_along(plots), function(i) {
   ggsave(filename = paste0("correlation_plot_", i, ".png"), plot = plots[[i]], width = 6, height = 4)
})


# Load required libraries
library(ggplot2)

# Example Dataset
# Replace 'mtcars' with your dataset
data <- mtcars

# Add a simulated `core` column and example columns for Shannon indices
data$core <- sample(c("yes", "no"), nrow(data), replace = TRUE)  # Simulating a "core" column
colnames(data)[1:2] <- c("m_plant_shannon", "m_microbe_shannon")  # Rename columns for this example

# Step 1: Filter data for rows where core == "yes"
filtered_data <- data[data$core == "yes", ]

# Step 2: Ensure the columns of interest exist
if (!all(c("m_plant_shannon", "m_microbe_shannon") %in% colnames(filtered_data))) {
  stop("Columns m_plant_shannon and m_microbe_shannon must be present in the data!")
}

# Step 3: Perform correlation test
cor_test <- cor.test(data$m_plant_shannon, data$m_microbe_shannon, use = "pairwise.complete.obs")

# Extract R^2 and p-value
r2 <- cor_test$estimate^2
p_value <- cor_test$p.value

# Step 4: Create scatter plot with R^2 and p-value
plot <- ggplot(data, aes(x = m_plant_shannon, y = m_microbe_shannon)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(title = "Correlation: m_plant_shannon vs m_microbe_shannon",
       subtitle = paste("R² =", round(r2, 2), "| p =", signif(p_value, 3)),
       x = "m_plant_shannon",
       y = "m_microbe_shannon") +
  theme_minimal()

# Step 5: Display the plot
print(plot)

# Optional: Save the plot to a file
# ggsave("correlation_m_plant_shannon_vs_m_microbe_shannon.png", plot, width = 6, height = 4)
md <- read.csv("metadata130225.csv")
#data <- read.csv("metadata070225.csv")

filtered_data <- data[md$core == "yes", ]
selected_data <- filtered_data[, grepl("^m_", colnames(filtered_data))]
#for the graminoids vs the co2 flux
cor_test <- cor.test(filtered_data$m_graminoid, filtered_data$CO2fluxdark, use = "pairwise.complete.obs")
# Extract R^2 and p-value
r2 <- cor_test$estimate^2
p_value <- cor_test$p.value
r2
p_value
# Step 3: Perform correlation test
cor_test <- cor.test(filtered_data1$m_microbe_shannon, filtered_data1$m_TOC, use = "pairwise.complete.obs")

# Extract R^2 and p-value
r2 <- cor_test$estimate^2
p_value <- cor_test$p.value


plotTOCs <- ggplot(filtered_data1, aes_string(x = filtered_data1$m_TOC, y = filtered_data1$m_microbe_shannon)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(title = paste("Correlation: TOC vs. shannon"),
       subtitle = paste("R² =", round(r2, 2), "| p =", signif(p_value, 3)),
       x = "TOC [%]",
       y = "Prokaryotic Shannon index") +
  theme_minimal()
plotTOCs

cor_test <- cor.test(filtered_data1$m_mean_qPCR_16s, filtered_data1$m_TOC, use = "pairwise.complete.obs")

# Extract R^2 and p-value
r2 <- cor_test$estimate^2
p_value <- cor_test$p.value


plotTOCa <- ggplot(filtered_data1, aes_string(x = filtered_data1$m_TOC, y = filtered_data1$m_mean_qPCR_16s)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(title = paste("Correlation: TOC vs. 16S qPCR"),
       subtitle = paste("R² =", round(r2, 2), "| p =", signif(p_value, 3)),
       x = "TOC [%]",
       y = "16S abundance [# copies/g of soil]") +
  theme_minimal()
plotTOCa

###all samples
md <- read.csv("metadata130225.csv")
#data <- read.csv("metadata070225.csv")

#filtered_data <- data[md$core == "yes", ]
selected_data <- md[, grepl("^m_", colnames(md))]
# Step 3: Perform correlation test
cor_test <- cor.test(md$m_microbe_shannon, md$m_TOC, use = "pairwise.complete.obs")

# Extract R^2 and p-value
r2 <- cor_test$estimate^2
p_value <- cor_test$p.value


plotTOCs <- ggplot(md, aes_string(x = md$m_TOC, y = md$m_microbe_shannon)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(title = paste("Correlation: TOC vs. shannon"),
       subtitle = paste("R² =", round(r2, 2), "| p =", signif(p_value, 3)),
       x = "TOC [%]",
       y = "Prokaryotic Shannon index") +
  theme_minimal()
plotTOCs

cor_test <- cor.test(md$m_mean_qPCR_16s, md$m_TOC, use = "pairwise.complete.obs")

# Extract R^2 and p-value
r2 <- cor_test$estimate^2
p_value <- cor_test$p.value


plotTOCa <- ggplot(md, aes_string(x = md$m_TOC, y = md$m_mean_qPCR_16s)) +
  geom_point(color = "blue", alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  labs(title = paste("Correlation: TOC vs. 16S qPCR"),
       subtitle = paste("R² =", round(r2, 2), "| p =", signif(p_value, 3)),
       x = "TOC [%]",
       y = "16S abundance [# copies/g of soil]") +
  theme_minimal()
plotTOCa
