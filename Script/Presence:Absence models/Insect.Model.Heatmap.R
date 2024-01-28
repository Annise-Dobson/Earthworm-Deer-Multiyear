#### Insect models - Presence/Absence (supplemental) ----
## Code by A.M. Dobson
# UPDATED January 26 2024

library(lme4)
library(lmerTest)

library(readr)
d4 <- read_csv("Data.csv")

# Scale the continuous predictor variables
d4$deer <- scale(d4$deer)
d4$worms <- scale(d4$worms)
d4$Year <- scale(d4$Year)
d4$Row <- as.factor(d4$Row)

# List of species
species_list <- c("Actaea", "Adiantum", "Allium", "Arisaema", "Brachyelytrum", "Caulophyllum","Dryopteris",
                  "Maianthemum", "Polygonatum","Polystichum", "Quercus","Sanguinaria","Thalictrum", "Tiarella","Trillium")

# Loop through each species
for (species in species_list) {
  insect_col <- paste0(species, ".i")
  I_col <- paste0("I.", species)
  
  # Check if the species-specific columns exist in the data frame
  if (!(insect_col %in% names(d4)) || !(I_col %in% names(d4))) {
    cat("Skipping", species, ": column not found.\n")
    next
  }
  
  # Scale species-specific columns
  d4[[I_col]] <- scale(d4[[I_col]])
  
  # Initialize flag for model convergence
  model_converged <- TRUE
  
  tryCatch({
    # Full model without the three-way interaction
    full_model_formula <- paste0(insect_col, " ~ Year * deer + Year * worms + deer * worms + I.", species, " + (1 | site/Row)")
    full_model <- glmer(as.formula(full_model_formula),family=binomial, data = d4)
    
    reduced_models_formulas <- c(
      paste0(insect_col, " ~ Year + deer + worms + I.", species, " + (1 | site/Row)"),
      paste0(insect_col, " ~ Year * deer + worms + I.", species, " + (1 | site/Row)"),
      paste0(insect_col, " ~ Year * worms + deer + I.", species, " + (1 | site/Row)"),
      paste0(insect_col, " ~ deer * worms + Year + I.", species, " + (1 | site/Row)")
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
      model_name <- paste0(species, "_best_insect_model")
      assign(model_name, best_model)
      cat("\nBest insect model for", species, ":\n")
      cat(best_model_formula, "\n")
      summary_model <- summary(get(model_name))
      print(summary_model)
      # Extract and print p-values
      cat("P-values for the fixed effects:\n")
      print(summary_model$coefficients[, "Pr(>|t|)"])
    } else {
      cat("No model was found to be more parsimonious than the full model within the significance threshold for", species, ".\n")
      model_name <- paste0(species, "_best_insect_model")
      assign(model_name, full_model)
      cat("The full insect model summary for", species, "is:\n")
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

library(lme4)
library(ggplot2)
library(reshape2)
library(lmerTest)  # For summary with p-values

format_estimate <- function(estimate, p.value) {
  if (p.value < 0.001) {
    return(paste0(formatC(estimate, format = 'f', digits = 2), "***"))
  } else if (p.value < 0.01) {
    return(paste0(formatC(estimate, format = 'f', digits = 2), "**"))
  } else if (p.value < 0.05) {
    return(paste0(formatC(estimate, format = 'f', digits = 2), "*"))
  } else {
    return(paste0(formatC(estimate, format = 'f', digits = 2)))  # Return the estimate even if non-significant
  }
}

get_formatted_estimate <- function(variable_name, parameter_estimates) {
  if (variable_name %in% rownames(parameter_estimates)) {
    estimate <- parameter_estimates[variable_name, "Estimate"]
    p.value <- parameter_estimates[variable_name, "Pr(>|t|)"]
    
    if (is.na(estimate)) {
      return("")  # Return blank for NAs
    } else if (p.value < 0.001) {
      return(paste0(formatC(estimate, format = 'f', digits = 2), "***"))
    } else if (p.value < 0.01) {
      return(paste0(formatC(estimate, format = 'f', digits = 2), "**"))
    } else if (p.value < 0.05) {
      return(paste0(formatC(estimate, format = 'f', digits = 2), "*"))
    } else {
      return(formatC(estimate, format = 'f', digits = 2))  # Return the estimate for non-significant estimates
    }
  } else {
    return("")  # Return blank for non-existing variables
  }
}



# Function to extract numeric values from formatted estimates
extract_numeric_value <- function(value) {
  as.numeric(gsub("[^0-9.-]", "", value))
}

# Initialize an empty dataframe for all species
all_species_data <- data.frame()

# Loop through each species for model extraction and data preparation
for (species in species_list) {
  model_name <- paste0(species, "_best_insect_model")
  if (!exists(model_name)) {
    next
  }
  best_model <- get(model_name)
  parameter_estimates <- summary(best_model)$coefficients
  
  data <- data.frame(
    Species = species,
    Deer = get_formatted_estimate("deer", parameter_estimates),
    Worms = get_formatted_estimate("worms", parameter_estimates),
    `Year x Deer` = get_formatted_estimate("Year:deer", parameter_estimates),
    `Year x Worms` = get_formatted_estimate("Year:worms", parameter_estimates),
    `Deer x Worms` = get_formatted_estimate("deer:worms", parameter_estimates)
  )
  all_species_data <- rbind(all_species_data, data)
}

# Convert formatted estimates to numeric, treating NAs and non-numeric characters as zeros
all_species_data_numeric <- all_species_data
all_species_data_numeric[,-1] <- lapply(all_species_data_numeric[,-1], function(x) ifelse(is.na(x), 0, as.numeric(gsub("[^0-9.-]", "", x))))

# Calculate impact
all_species_data_numeric$impact <- rowSums(all_species_data_numeric[,-1], na.rm = TRUE)

# Add the impact column to the original all_species_data
all_species_data$impact <- all_species_data_numeric$impact

# Reorder species by 'impact' in descending order
all_species_data <- all_species_data[order(all_species_data$impact, decreasing = TRUE), ]

# Reshape data for plotting
data_long <- melt(all_species_data, id.vars = c("Species", "impact"))

# Adjust the factor levels of Species based on the reverse order in all_species_data
data_long$Species <- factor(data_long$Species, levels = rev(all_species_data$Species))

# Create and reorder the formatted variable names
data_long$formatted_variable <- gsub("Year:deer", "Year x Deer", 
                                     gsub("Year:worms", "Year x Worms", 
                                          gsub("deer:worms", "Deer x Worms", 
                                               gsub("\\.", " ", data_long$variable))))
ordered_variables <- c("Deer", "Worms", "Year x Deer", "Year x Worms", "Deer x Worms")
data_long$formatted_variable <- factor(data_long$formatted_variable, levels = ordered_variables)

# Create a new column for significance
data_long$is_significant <- grepl("\\*", data_long$value)

# Extract numeric values for the heatmap
data_long$value_numeric <- as.numeric(gsub("[^0-9.-]", "", data_long$value))

library(scales)  # Load the scales package

ggplot(data_long, aes(x = formatted_variable, y = Species, fill = ifelse(is_significant, value_numeric, NA))) +
  geom_tile() +
  geom_text(aes(label = value), color = "black", na.rm = TRUE) +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, na.value = "white", 
                       limits = c(min(data_long$value_numeric, na.rm = TRUE), max(data_long$value_numeric, na.rm = TRUE)),
                       labels = label_number(accuracy = 0.1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(face = "italic", size = 12),
        legend.position = "top",
        legend.text = element_text(size = 10)) +
  xlab("") +
  ylab("") +
  labs(fill = "Scaled Coefficient")
