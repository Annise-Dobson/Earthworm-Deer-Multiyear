#### Heatmap of survival models ----
## Code by A.M. Dobson
# UPDATED January 26 2024

library(lme4)
library(lmerTest)
library(MuMIn)
library(readr)
d4 <- read_csv("/Users/annisedobson/Dropbox (YSE)/Data/Lumb x deer/DeerxEW/AllMeasurements.2017.final.csv")

# Scale the continuous predictor variables
d4$deer <- scale(d4$deer)
d4$Drybiomass.mean <- scale(d4$Drybiomass.mean)
d4$Year <- scale(d4$Year)
d4$Row <- as.factor(d4$Row)

# List of species
species_list <- c("Actaea","Adiantum","Geum","Agrimonia","Allium","Arisaema", "Brachyelytrum","Carex","Caulophyllum",
                  "Dryopteris","Geranium","Maianthemum","Polygonatum","Polygonum","Polystichum","Quercus",
                  "Sanguinaria","Thalictrum","Tiarella","Trillium")

species_titles <- c(
  "Actaea rubra",
  "Adiantum pedatum",
  "Allium canadense",
  "Agrimonia gryposepala",
  "Arisaema triphyllum",
  "Brachyelytrum erectum",
  "Carex radiata",
  "Caulophyllum thalictroides",
  "Dryopteris sp.",
  "Geranium maculatum",
  "Geum canadense",
  "Maianthemum racemosum",
  "Polygonatum biflorum",
  "Polygonum virginianum",
  "Polystichum acrostichoides",
  "Quercus rubra",
  "Sanguinaria canadensis",
  "Thalictrum dioicum",
  "Tiarella cordifolia",
  "Trillium erectum"
)

# Loop through each species
for (species in species_list) {
  survival_col <- paste0(species, ".final")
  I_col <- paste0("I.", species)
  
  # Check if the species-specific columns exist in the data frame
  if (!(survival_col %in% names(d4)) || !(I_col %in% names(d4))) {
    cat("Skipping", species, ": column not found.\n")
    next
  }
  
  # Scale species-specific columns
  d4[[I_col]] <- scale(d4[[I_col]])
  
  # Initialize flag for model convergence
  model_converged <- TRUE
  
  tryCatch({
    # Full model without the three-way interaction
    full_model_formula <- paste0(survival_col, " ~ Year * deer + Year * Drybiomass.mean + deer * Drybiomass.mean + I.", species, " + (1 | site/Row)")
    full_model <- glmer(as.formula(full_model_formula), family=binomial, data = d4)
    
    reduced_models_formulas <- c(
      paste0(survival_col, " ~ Year + deer + Drybiomass.mean + I.", species, " + (1 | site/Row)"),
      paste0(survival_col, " ~ Year * deer + Drybiomass.mean + I.", species, " + (1 | site/Row)"),
      paste0(survival_col, " ~ Year * Drybiomass.mean + deer + I.", species, " + (1 | site/Row)"),
      paste0(survival_col, " ~ deer * Drybiomass.mean + Year + I.", species, " + (1 | site/Row)")
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
      model_name <- paste0(species, "_best_survival_model")
      assign(model_name, best_model)
      cat("\nBest survival model for", species, ":\n")
      cat(best_model_formula, "\n")
      summary_model <- summary(get(model_name))
      print(summary_model)
      # Extract and print p-values
      cat("P-values for the fixed effects:\n")
      print(summary_model$coefficients[, "Pr(>|t|)"])
    } else {
      cat("No model was found to be more parsimonious than the full model within the significance threshold for", species, ".\n")
      model_name <- paste0(species, "_best_survival_model")
      assign(model_name, full_model)
      cat("The full survival model summary for", species, "is:\n")
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
  model_name <- paste0(species, "_best_survival_model")
  if (!exists(model_name)) {
    next
  }
  best_model <- get(model_name)
  parameter_estimates <- summary(best_model)$coefficients
  
  data <- data.frame(
    Species = species,
    Deer = get_formatted_estimate("deer", parameter_estimates),
    Worms = get_formatted_estimate("Drybiomass.mean", parameter_estimates),
    `Year x Deer` = get_formatted_estimate("Year:deer", parameter_estimates),
    `Year x Worms` = get_formatted_estimate("Year:Drybiomass.mean", parameter_estimates),
    `Deer x Worms` = get_formatted_estimate("deer:Drybiomass.mean", parameter_estimates)
  )
  all_species_data <- rbind(all_species_data, data)
}

# Convert formatted estimates to numeric, treating NAs and non-numeric characters as zeros
all_species_data_numeric <- all_species_data
all_species_data_numeric[,-1] <- lapply(all_species_data_numeric[,-1], function(x) ifelse(is.na(x), 0, extract_numeric_value(x)))

# Calculate impact
all_species_data_numeric$impact <- rowSums(all_species_data_numeric[,-1], na.rm = TRUE)

# Add the impact column to the original all_species_data
all_species_data$impact <- all_species_data_numeric$impact

# Reorder species by 'impact' in descending order
all_species_data <- all_species_data[order(all_species_data$impact, decreasing = FALSE), ]

# Reshape data for plotting
library(reshape2)
data_long <- melt(all_species_data, id.vars = c("Species", "impact"))

# Map species titles to the Species factor for use in plots
species_title_mapping <- setNames(species_titles, species_list)
data_long$Species_Title <- species_title_mapping[data_long$Species]

# Adjust the factor levels of Species_Title based on the order in all_species_data
data_long$Species_Title <- factor(data_long$Species_Title, levels = species_title_mapping[all_species_data$Species])

# Create and reorder the formatted variable names
data_long$formatted_variable <- gsub("Year:deer", "Year x Deer", 
                                     gsub("Year:Drybiomass.mean", "Year x Worms", 
                                          gsub("deer:Drybiomass.mean", "Deer x Worms", 
                                               gsub("\\.", " ", data_long$variable))))
ordered_variables <- c("Deer", "Worms", "Year x Deer", "Year x Worms", "Deer x Worms")
data_long$formatted_variable <- factor(data_long$formatted_variable, levels = ordered_variables)

# Create a new column for significance
data_long$is_significant <- grepl("\\*", data_long$value)

# Extract numeric values for the heatmap
data_long$value_numeric <- as.numeric(gsub("[^0-9.-]", "", data_long$value))

# Load required libraries
library(ggplot2)
library(scales)

# Plotting
ggplot(data_long, aes(x = formatted_variable, y = Species_Title, fill = ifelse(is_significant, value_numeric, NA))) +
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

############r2
# Initialize a dataframe for R-squared values
r_squared_df <- data.frame(Species = character(), Marginal_R2 = numeric(), Conditional_R2 = numeric())

# Loop through each species
for (species in species_list) {
  survival_col <- paste0(species, ".final")
  I_col <- paste0("I.", species)
  
  # Check if the species-specific columns exist in the data frame
  if (!(survival_col %in% names(d4)) || !(I_col %in% names(d4))) {
    cat("Skipping", species, ": column not found.\n")
    next
  }
  
  # Scale species-specific columns
  d4[[I_col]] <- scale(d4[[I_col]])
  
  # Full model without the three-way interaction
  full_model_formula <- paste0(survival_col, " ~ Year * deer + Year * Drybiomass.mean + deer * Drybiomass.mean + I.", species, " + (1 | site/Row)")
  full_model <- lmer(as.formula(full_model_formula), data = d4)
  
  # Compute R-squared values
  tryCatch({
    r_squared <- r.squaredGLMM(full_model)
    if (is.numeric(r_squared) && length(r_squared) >= 2) {
      # Append the R2 values to the dataframe
      r_squared_df <- rbind(r_squared_df, data.frame(Species = species, Marginal_R2 = r_squared[1], Conditional_R2 = r_squared[2]))
    } else {
      warning("Unexpected format of R-squared values for species: ", species)
    }
  }, error = function(e) {
    cat("Error in computing R-squared for", species, ": ", e$message, "\n")
  })
}

# Round the R-squared values to three significant digits
r_squared_df$Marginal_R2 <- round(r_squared_df$Marginal_R2, 3)
r_squared_df$Conditional_R2 <- round(r_squared_df$Conditional_R2, 3)

# Print the R-squared values in tab-delimited format for copy-pasting
write.table(r_squared_df, row.names = FALSE, sep = "\t", quote = FALSE)

######
###########
library(lme4)
library(lmerTest)
species_list <- c("Actaea","Adiantum","Geum","Agrimonia","Allium","Arisaema", "Brachyelytrum","Carex","Caulophyllum",
                  "Dryopteris","Geranium","Maianthemum","Polygonatum","Polygonum","Polystichum","Quercus",
                  "Sanguinaria","Thalictrum","Tiarella","Trillium")
# Define the function to extract model results
extract_model_results <- function(model_name) {
  if (!exists(model_name)) {
    stop("Model does not exist: ", model_name)
  }
  
  # Retrieve the model
  model <- get(model_name)
  
  # Extract the summary of the model
  model_summary <- summary(model)
  
  # Check if the model is a glmer or lmer model
  if ("z value" %in% colnames(model_summary$coefficients)) {
    # glmer model
    stat_col = "z value"
    p_val_col = "Pr(>|z|)"
  } else {
    # lmer model
    stat_col = "t value"
    p_val_col = "Pr(>|t|)"
  }
  
  # Create a dataframe of the fixed effects results
  fixed_effects <- model_summary$coefficients
  data.frame(
    Parameter = rownames(fixed_effects),
    Estimate = fixed_effects[, "Estimate"],
    SE = fixed_effects[, "Std. Error"],
    Statistic = fixed_effects[, stat_col],
    P_value = fixed_effects[, p_val_col]
  )
}

# Initialize an empty data frame to store results for all models
all_model_results <- data.frame(Species = character(), Parameter = character(), 
                                Estimate = numeric(), SE = numeric(), 
                                Statistic = numeric(), P_value = numeric())

# Loop through each species in your species list
for (species in species_list) {
  model_name <- paste0(species, "_best_survival_model")
  model_results <- tryCatch({
    extract_model_results(model_name)
  }, error = function(e) {
    message("Error extracting results for ", species, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(model_results)) {
    model_results$Species <- species
    all_model_results <- rbind(all_model_results, model_results)
  }
}

# Print the results table
print(all_model_results, row.names = FALSE)

# Print the results in tab-delimited format for copying into Excel
write.table(all_model_results, row.names = FALSE, sep = "\t", quote = FALSE)

#########
# Initialize a data frame to store the results
random_effects_stats <- data.frame(Species = character(), 
                                   Random_Effect = character(),
                                   Variance = numeric(), 
                                   SD = numeric())

# Loop through each species
for (species in species_list) {
  model_name <- paste0(species, "_best_survival_model")
  
  if (exists(model_name)) {
    # Get the model
    model <- get(model_name)
    
    # Extract variance components
    variance_components <- VarCorr(model)
    
    # Loop through each random effect
    for (effect in names(variance_components)) {
      var <- variance_components[[effect]]@.Data[1] # Variance
      sd <- sqrt(var) # Standard deviation
      
      # Append to the data frame
      random_effects_stats <- rbind(random_effects_stats, 
                                    data.frame(Species = species, 
                                               Random_Effect = effect, 
                                               Variance = var, 
                                               SD = sd))
    }
  }
}

# Print or save the results
write.table(random_effects_stats, row.names = FALSE, sep = "\t", quote = FALSE)


# Assuming your data frame is named d4 and has the necessary variables
# Replace 'worms' with 'Drybiomass.mean' if that's your variable name for worm biomass
# Scale the continuous predictor variables

d4$deer <- scale(d4$deer)
d4$Drybiomass.mean <- scale(d4$Drybiomass.mean)
d4$Year <- scale(d4$Year)
d4$Row <- as.factor(d4$Row)
d4$I.Carex <- scale(d4$I.Carex)

# Fit the model (assuming you have already handled NA values and other preprocessing)
carex_model <- lmer(Carex.NF ~ I.Carex + deer + worms + Year +
                      (1 | site/Row), 
                    data = d4)

# Get the summary of the model
model_summary <- summary(carex_model)

# Extracting fixed effects
fixed_effects <- model_summary$coefficients
fixed_effects_df <- data.frame(Estimate = fixed_effects[, "Estimate"],
                               Std.Error = fixed_effects[, "Std. Error"],
                               t.value = fixed_effects[, "t value"],
                               Pr = fixed_effects[, "Pr(>|t|)"])

# Printing the fixed effects in a tab-delimited format
write.table(fixed_effects_df, row.names = TRUE, sep = "\t", quote = FALSE)

# Similarly, you can extract and print random effects and other components as needed
