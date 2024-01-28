#### Width figure ----
## Code by A.M. Dobson
# UPDATED January 26 2024

d4 <- read_csv("Data.csv")

# Load the necessary libraries
if (!require(lme4)) {
  install.packages("lme4")
}
if (!require(lmerTest)) {
  install.packages("lmerTest")
}
library(lme4)
library(lmerTest)
library(readr)

# Scale the continuous predictor variables
d4$deer <- as.factor(d4$deer)
# Rename the factor levels
levels(d4$deer) <- c("Fenced", "Unfenced")
d4$worms <- as.factor(d4$worms)
d4$Year.order <- as.numeric(d4$Year.order)
d4$Row <- as.factor(d4$Row)

# List of species
species_list <- c("Actaea", "Adiantum", "Allium","Caulophyllum", "Polygonatum", "Quercus", 
                  "Sanguinaria", "Tiarella", "Trillium")

species_titles <- c(
  "Actaea rubra",
  "Adiantum pedatum",
  "Allium canadense",
  "Caulophyllum thalictroides",
  "Polygonatum biflorum",
  "Quercus rubra",
  "Sanguinaria canadensis",
  "Tiarella cordifolia",
  "Trillium erectum"
)

names(species_titles) <- species_list

# Loop through each species
for (species in species_list) {
  width_col <- paste0(species, ".W")
  I_col <- paste0("I.", species)
  
  # Check if the species-specific columns exist in the data frame
  if (!(width_col %in% names(d4)) || !(I_col %in% names(d4))) {
    cat("Skipping", species, ": column not found.\n")
    next
  }
  
  # Scale species-specific columns
  d4[[I_col]] <- scale(d4[[I_col]])
  
  # Initialize flag for model convergence
  model_converged <- TRUE
  
  tryCatch({
    # Full model formula without three-way interaction
    full_model_formula <- paste0(width_col, " ~ Year.order * deer + Year.order * worms + deer * worms + I.", species, " + (1 | site/Row)")
    full_model <- lmer(as.formula(full_model_formula), data = d4)
    
    reduced_models_formulas <- c(
      paste0(width_col, " ~ Year.order + deer + worms + I.", species, " + (1 | site/Row)"),
      paste0(width_col, " ~ Year.order * deer + Year.order * worms + I.", species, " + (1 | site/Row)"),
      paste0(width_col, " ~ Year.order * worms + deer * worms + I.", species, " + (1 | site/Row)"),
      paste0(width_col, " ~ deer * worms + Year.order + I.", species, " + (1 | site/Row)")
    )
    
    best_model <- NULL
    best_model_aic <- Inf
    best_model_formula <- ""
    
    for (formula in reduced_models_formulas) {
      reduced_model <- lmer(as.formula(formula), data = d4)
      aic_value <- AIC(reduced_model)
      
      if (aic_value < best_model_aic) {
        best_model <- reduced_model
        best_model_aic <- aic_value
        best_model_formula <- formula
      }
    }
    
    if (!is.null(best_model)) {
      model_name <- paste0(species, "_best_width_model")
      assign(model_name, best_model)
      cat("\nBest width model for", species, ":\n")
      cat(best_model_formula, "\n")
      summary_model <- summary(get(model_name))
      print(summary_model)
      # Extract and print p-values
      cat("P-values for the fixed effects:\n")
      print(summary_model$coefficients[, "Pr(>|t|)"])
    } else {
      cat("No model was found to be more parsimonious than the full model within the significance threshold for", species, ".\n")
      model_name <- paste0(species, "_best_width_model")
      assign(model_name, full_model)
      cat("The full width model summary for", species, "is:\n")
      summary_model <- summary(get(model_name))
      print(summary_model)
      # Extract and print p-values
      cat("P-values for the fixed effects:\n")
      print(summary_model$coefficients[, "Pr(>|t|)"])
    }
    
  }, warning = function(w) {
    cat("Warning in fitting model for", species, ": ", w$message, "\n")
  }, error = function(e) {
    cat("Error in fitting model for", species, ": ", e$message, "\n")
    model_converged <- FALSE
  })
  
  if (!model_converged) {
    cat("Skipping to next species due to non-convergence.\n")
    next
  }
}

###################
# Ensure necessary libraries are loaded
library(emmeans)
library(ggplot2)
library(gridExtra)

# Define groups with different formatting
group_x_axis <- c("Sanguinaria", "Tiarella", "Trillium")  # x and y
group_y_axis <- c("Actaea", "Adiantum", "Allium", "Caulophyllum", "Polygonatum", "Quercus", "Sanguinaria", "Tiarella", "Trillium")  # y only
group_grey_strip <- c("")  # Species with special strip background

# Initialize list to store plots
plot_list <- list()

# Loop through each species to generate plots
for (species in species_list) {
  # Dynamically create model name and retrieve the model
  model_name <- paste0(species, "_best_width_model")
  
  # Check if the model exists
  if (!exists(model_name)) {
    next
  }
  
  best_model <- get(model_name) # Retrieve model by name
  
  # Use the species title as the plot title
  species_title <- species_titles[species]
  
  # Generate the plot with emmip
  plot <- emmip(best_model, as.factor(worms) ~ Year.order | as.factor(deer), 
                at = list(Year.order = c(0:5)), type = "response",
                CIs = TRUE, engine = "ggplot") + 
    scale_color_manual(values = c("#000000", "#E69F00"))
  
  # Apply a different theme when species is not in group_grey_strip
  if (species %in% group_grey_strip) {
    plot <- plot + theme(strip.background = element_rect(fill = "lightgrey"))
  } else {
    plot <- plot + theme(strip.background = element_blank(), strip.text.x = element_blank())
  }
  
  # Apply conditional formatting based on groups
  if (species %in% group_x_axis) {
    plot <- plot + scale_x_continuous(breaks = c(0:5), labels = c("0", "1", "2", "3", "4", "5"))
  } else if (species %in% group_y_axis) {
    plot <- plot + scale_x_continuous(breaks = c(0:5), labels = NULL)
  }
  
  # Common theme adjustments
  plot <- plot + 
    labs(x = NULL, y = NULL, title = bquote(italic(.(species_title)))) + # Use species_title here
    theme(plot.title = element_text(size = 10),
          legend.position = "none",
          plot.margin = unit(c(0.2, 0.3, 0.1, 0.1), "cm"),                          
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(size = .4),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  # Store the plot in the list
  plot_list[[species]] <- plot
}

# Combine all plots into a single multiplot
multiplot <- grid.arrange(grobs = plot_list, ncol = 3) # Adjust the number of columns as needed

# Save as a high-quality TIFF image
ggsave("/Users/annisedobson/Dropbox (YSE)/Data/Lumb x deer/DeerxEW/Final figures/Width.tiff", 
       plot = multiplot, 
       width = 10, 
       height = 8, 
       dpi = 300, 
       device = "tiff",
       compression = "none")
