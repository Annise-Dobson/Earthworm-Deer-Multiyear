#### Insect figure ----
## Code by A.M. Dobson
# UPDATED January 26 2024

d4 <- read_csv("Data.csv")

library(emmeans)
library(ggplot2)
library(gridExtra)

# Species list
##Need to run once without Deer x Year to Let Drop converge and vary 
# List of species
species_list <- c("Actaea", "Adiantum", "Allium", "Arisaema", "Brachyelytrum", "Caulophyllum","Dryopteris",
                  "Maianthemum", "Polygonatum","Polystichum", "Quercus","Sanguinaria","Thalictrum", "Tiarella","Trillium")

species_titles <- c(
  "Actaea rubra",
  "Adiantum pedatum",
  "Allium canadense",
  "Arisaema triphyllum",
  "Brachyelytrum erectum",
  "Caulophyllum thalictroides",
  "Dryopteris sp.",
  "Maianthemum racemosum",
  "Polygonatum biflorum",
  "Polystichum acrostichoides",
  "Quercus rubra",
  "Sanguinaria canadensis",
  "Thalictrum dioicum",
  "Tiarella cordifolia",
  "Trillium erectum"
)

names(species_titles) <- species_list

# Load the necessary library (if not already loaded)
if (!require(lme4)) {
  install.packages("lme4")
}
library(lme4)

# Scale the continuous predictor variables
d4$deer <- as.factor(d4$deer)
# Rename the factor levels
levels(d4$deer) <- c("Fenced", "Unfenced")
d4$worms <- as.factor(d4$worms)
d4$Year.order <- as.numeric(d4$Year.order)
d4$Row <- as.factor(d4$Row)
d4$pH <- as.factor(d4$pH)

# Loop through each species
for (species in species_list) {
  final_col <- paste0(species, ".i")
  I_col <- paste0("I.", species)
  
  # Check if the species-specific columns exist in the data frame
  if (!(final_col %in% names(d4)) || !(I_col %in% names(d4))) {
    cat("Skipping", species, ": column not found.\n")
    next
  }
  
  # Factorize and scale species-specific columns
  d4[[final_col]] <- factor(d4[[final_col]])
  d4[[I_col]] <- scale(d4[[I_col]])
  
  # Initialize flag for model convergence
  model_converged <- TRUE
  
  tryCatch({
    # Model formula with two-way interactions
    full_model_formula <- paste0(species, ".i ~ Year.order * deer + Year.order * worms + deer * worms + I.", species, " + (1 | site)")
    full_model <- glmer(as.formula(full_model_formula), data = d4, family = binomial, control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1000000)))
    
    # Reduced models without the three-way interaction
    reduced_models_formulas <- c(
      paste0(species, ".i ~ Year.order + deer + worms + I.", species, " + (1 | site)"),
      paste0(species, ".i ~ Year.order * deer + Year.order + worms + I.", species, " + (1 | site)"),
      paste0(species, ".i ~ Year.order * worms + deer + worms + I.", species, " + (1 | site)"),
      paste0(species, ".i ~ deer * worms + Year.order + deer + I.", species, " + (1 | site)")
    )
    
    best_model <- NULL
    best_model_p_value <- 0
    best_model_formula <- ""
    
    for (formula in reduced_models_formulas) {
      reduced_model <- glmer(as.formula(formula), data = d4, family = binomial)
      lr_test <- anova(reduced_model, full_model)
      p_value <- lr_test$`Pr(>Chisq)`[2] # get the p-value for the test
      
      if (!is.na(p_value) && p_value > best_model_p_value && p_value > 0.05) {
        best_model <- reduced_model
        best_model_p_value <- p_value
        best_model_formula <- formula
      }
    }
    
    if (!is.null(best_model)) {
      model_name <- paste0(species, "_best_model")
      assign(model_name, best_model)
      cat("\nBest model for", species, ":\n")
      cat(best_model_formula, "\n")
      print(summary(get(model_name)))
    } else {
      cat("No model was found to be more parsimonious than the full model within the significance threshold for", species, ".\n")
      model_name <- paste0(species, "_best_model")
      assign(model_name, full_model)
      cat("The full model summary for", species, "is:\n")
      print(summary(get(model_name)))
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

# Initialize an empty list to store plots
plot_list <- list()

# Define groups with different formatting
group1 <- c("Actaea", "Arisaema", "Dryopteris","Polystichum")  # y-axis numbers only
group2 <- c("Thalictrum")  # y-axis and x-axis (years) numbers
group3 <- c("Tiarella", "Trillium")  # x-axis (years) numbers only
special_strip_group <- c("") # Species with special strip background

# Loop through each species to generate plots
for(species in species_list) {
  # Retrieve the model by dynamically created model name
  model_name <- paste0(species, "_best_model")
  best_model <- get(model_name) # Retrieve model by name
  
  # Use the species title as the plot title
  species_title <- species_titles[species]
  
  # Initialize plot with common settings
  plot <- emmip(best_model, as.factor(worms) ~ Year.order | as.factor(deer), 
                at = list(Year.order = c(0:5)), type = "response",
                CIs = TRUE, engine = "ggplot") + 
    scale_color_manual(values = c("#000000", "#E69F00")) +
    guides(shape = guide_legend(override.aes = list(size = 5))) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0.5, 6)) +    
    labs(x = NULL, y = NULL, title = bquote(italic(.(species_title)))) + # Use species_title here
    theme(plot.title = element_text(size = 10),
          legend.position = "none",
          plot.margin = unit(c(0.2, 0.3, 0.1, 0.1), "cm"),   
          plot.background = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(), 
          panel.background = element_blank(),
          axis.line = element_line(size = .4))
  
  # Conditional theme adjustments for strip
  if (species %in% special_strip_group) {
    plot <- plot + 
      theme(strip.background = element_rect(fill = "lightgrey"), # Light grey for special group
            strip.text = element_text())
  } else {
    plot <- plot + 
      theme(strip.background = element_blank(), # Remove for other groups
            strip.text = element_blank())
  }
  
  # Apply conditional formatting based on groups
  if (species %in% group1) {
    plot <- plot + 
      scale_x_discrete(limits = c(0:5), labels = NULL) +
      scale_y_continuous(labels = scales::number_format())
  } else if (species %in% group2) {
    plot <- plot + 
      scale_x_continuous(breaks = c(0:5), labels = c("0", "1", "2", "3", "4", "5")) +
      scale_y_continuous(labels = scales::number_format())
  } else if (species %in% group3) {
    plot <- plot + 
      scale_x_continuous(breaks = c(0:5), labels = c("0", "1", "2", "3", "4", "5")) +
      scale_y_continuous(labels = NULL)
  } else {
    plot <- plot + 
      scale_x_discrete(limits = c(0:5), labels = NULL) +
      scale_y_continuous(labels = NULL)
  }
  
  # Store the plot in the list
  plot_list[[species]] = plot
}

# Combine all plots into a single multiplot
multiplot <- grid.arrange(grobs = plot_list, ncol = 3) # Adjust the number of columns as needed

# Save as a high-quality TIFF image
ggsave("/Users/annisedobson/Dropbox (YSE)/Data/Lumb x deer/DeerxEW/Final figures/Insect.tiff", 
       plot = multiplot, 
       width = 10, 
       height = 8, 
       dpi = 300, 
       device = "tiff",
       compression = "none")
