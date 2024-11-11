
suppressMessages(library(tidyverse))
suppressMessages(library(usmap))
suppressMessages(library(scales))
suppressMessages(library(mice))
suppressMessages(library(glmnet))
suppressMessages(library(boot))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))



# Read data in and seprate State and County in to seperate columns, then edit any problematic entries and then create a new column of FIPS codes, that uniquelly identify a county
cancer_raw = read.csv("cancer_reg.csv")
cancer_edit = cancer_raw
cancer_edit <- cancer_edit %>% mutate(Target_div_Income = TARGET_deathRate/medIncome)
cancer_geo = cbind(cancer_edit, str_match(cancer_edit$Geography,"(.+), (.+)")[ ,-1])
colnames(cancer_geo)[37] ="State"
colnames(cancer_geo)[36] = "County"
cancer_geo[167,36] <- "Dona Ana County"
cancer_geo[821,36] <- "La Salle Parish"
codes <- rep(NULL, length(cancer_geo$County))

for (i in 1:length(cancer_geo$avgAnnCount)){
 codes[i] = fips(state = cancer_geo$State[i], county = cancer_geo$County[i])
}
cancer_final = cbind(cancer_geo, fips = codes)



moddat <- cancer_final

(colMeans(is.na(moddat)))*100



# Set reproducability seed and then impute data
set.seed(1)
trim = moddat[,-18]
imp <- mice(trim, m = 5, maxit = 50, meth = "pmm")
complete(imp)



imputed <- complete(imp)
imputed_new <- imputed



mod1 <- lm(data = imputed_new, TARGET_deathRate ~ povertyPercent + PctBlack + PctNoHS18_24 + PctHS18_24)



summary(mod1)



#First, code the Southeast variable for future use
new_england <- c("Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island", "Vermont")
mideast <- c("Delaware", "District of Columbia", "Maryland", "New Jersey", "New York", "Pennsylvania")
great_lakes <- c("Illinois", "Indiana", "Michigan", "Ohio", "Wisconsin")
plains <- c("Iowa", "Kansas", "Minnesota", "Missouri", "Nebraska", "North Dakota", "South Dakota")
southeast <- c("Alabama", "Arkansas", "Florida", "Georgia", "Kentucky", "Louisiana", "Mississippi", "North Carolina", "South Carolina", "Tennessee", "Virginia", "West Virginia")
southwest <- c("Arizona", "New Mexico", "Oklahoma", "Texas")
rocky_mountain <- c("Colorado", "Idaho", "Montana", "Utah", "Wyoming")
far_west <- c("Alaska", "California", "Hawaii", "Nevada", "Oregon", "Washington")

get_region <- function(state) {
  if (state %in% new_england) {
    return("New England")
  } else if (state %in% mideast) {
    return("Mideast")
  } else if (state %in% great_lakes) {
    return("Great Lakes")
  } else if (state %in% plains) {
    return("Plains")
  } else if (state %in% southeast) {
    return("Southeast")
  } else if (state %in% southwest) {
    return("Southwest")
  } else if (state %in% rocky_mountain) {
    return("Rocky Mountain")
  } else if (state %in% far_west) {
    return("Far West")
  } else {
    return(NA)
  }
}

imputed_new$Region <- sapply(imputed_new$State, get_region)

imputed_new$isSoutheast <- ifelse(imputed_new$Region == "Southeast", "Yes", "No")



#Create Lasso Lambda graph
set.seed(1)
y = imputed_new$TARGET_deathRate
x = data.matrix(imputed_new[, c('povertyPercent', 'PctBlack', 'PctHS18_24','PctNoHS18_24', 'isSoutheast','PctPublicCoverage','PctPublicCoverageAlone', "PctUnemployed16_Over")])
cv_model <- cv.glmnet(x, y, alpha = 1)
plot(cv_model)



#Assign 1se lambda and then run Lasso using it
min_lambda <- cv_model$lambda.min
se_lambda <- cv_model$lambda.1se
best_model <- glmnet(x, y, alpha = 1, lambda = se_lambda)
coef(best_model)



finmod <- lm(data = imputed_new, TARGET_deathRate ~ povertyPercent + PctHS18_24 + isSoutheast + PctPublicCoverage + PctPublicCoverageAlone + PctUnemployed16_Over)
summary(finmod)



plot(finmod, which =2)



shapiro.test(finmod$residuals)




nboot <- 10000
set.seed(1)

# Create a function to calculate the coefficients using the bootstrap
coef.boot <- function(data, indices) {
  model <- lm(TARGET_deathRate ~ povertyPercent + PctHS18_24 + isSoutheast + PctPublicCoverage + PctPublicCoverageAlone + PctUnemployed16_Over, data = data[indices, ])
  return(coef(model)[-1]) # exclude intercept column
}

# Perform the bootstrap using the defined function
boot.results <- boot(data = imputed_new, statistic = coef.boot, R = nboot)

# Convert bootstrap results to a data frame
boot.df <- as.data.frame(boot.results$t)
colnames(boot.df) <- c("povertyPercent", "PctHS18_24", "isSoutheastYes","PctPublicCoverage", "PctPublicCoverageAlone", "PctUnemployed16_Over")

# Get coefficient estimates from original model
finmod <- lm(TARGET_deathRate ~ povertyPercent + PctHS18_24 + isSoutheast + PctPublicCoverage + PctPublicCoverageAlone + PctUnemployed16_Over, data = imputed_new)
coef.estimates <- coef(finmod)[-1 , drop = TRUE]

# Create a function to plot histograms with quantile and coefficient lines and a title
plot.hist <- function(x, coef.est, varname) {
  p <- ggplot(data.frame(x), aes(x = x)) + 
    geom_histogram(binwidth = 0.05, color = "black", fill = "white") +
    geom_vline(xintercept = quantile(x, probs = c(0.004166667, 0.9958333)), linetype = "dashed") +
    geom_vline(xintercept = coef.est, color = "red", linetype = "dashed") +
    xlab("") + ylab("Frequency") +
    ggtitle(varname)+
    theme(plot.title = element_text(size = 9.5))
  return(p)
}

# Create a list of plots for each column in boot.df with titles
plot.list <- mapply(plot.hist, x = boot.df, coef.est = coef.estimates, varname = names(boot.df), SIMPLIFY = FALSE)

# Combine the plots into a single figure with a title
grid.arrange(grobs = plot.list, ncol = 3, top = textGrob('Histograms of Coefficient Estimates with 95% Confidence Intervals \n using a Bonferroni Correction for 6 Tests', gp=gpar(fontsize=13, fontface = 2))
)




# Create reduced data set and apply the Southeast column used previously
reduced <- na.omit(trim[,-21])
nrow(reduced) - nrow(imputed_new)

reduced$Region <- sapply(reduced$State, get_region)

reduced$isSoutheast <- ifelse(reduced$Region == "Southeast", "Yes", "No")



set.seed(1)
##DO NOT INCLUDE IN FINAL
nboot <- 10000

# Create a function to calculate the coefficients using the bootstrap
coef.boot <- function(data, indices) {
  model <- lm(TARGET_deathRate ~ povertyPercent + PctHS18_24 + isSoutheast + PctPublicCoverage + PctPublicCoverageAlone + PctUnemployed16_Over, data = data[indices, ])
  return(coef(model)[-1]) # exclude intercept column
}

# Perform the bootstrap using the defined function
boot.results.reduce <- boot(data = reduced, statistic = coef.boot, R = nboot)

# Convert bootstrap results to a data frame
boot.df.reduce <- as.data.frame(boot.results.reduce$t)
colnames(boot.df.reduce) <- c("povertyPercent", "PctHS18_24", "isSoutheastYes","PctPublicCoverage", "PctPublicCoverageAlone", "PctUnemployed16_Over")

# Get coefficient estimates from original model
finmod.reduce <- lm(TARGET_deathRate ~ povertyPercent + PctHS18_24 + isSoutheast + PctPublicCoverage + PctPublicCoverageAlone + PctUnemployed16_Over, data = reduced)
coef.estimates.reduce <- coef(finmod.reduce)[-1 , drop = TRUE]

# Create a function to plot histograms with quantile and coefficient lines and a title
plot.hist <- function(x, coef.est, varname) {
  p <- ggplot(data.frame(x), aes(x = x)) + 
    geom_histogram(binwidth = 0.05, color = "black", fill = "white") +
    geom_vline(xintercept = quantile(x, probs = c(0.004166667, 0.9958333)), linetype = "dashed") +
    geom_vline(xintercept = coef.est, color = "red", linetype = "dashed") +
    xlab("") + ylab("Frequency") +
    ggtitle(varname)+
    theme(plot.title = element_text(size = 9.5))
  return(p)
}

# Create a list of plots for each column in boot.df with titles
plot.list <- mapply(plot.hist, x = boot.df.reduce, coef.est = coef.estimates.reduce, varname = names(boot.df.reduce), SIMPLIFY = FALSE)

# Combine the plots into a single figure with a title
grid.arrange(grobs = plot.list, ncol = 3, top = textGrob('Reduced Data Histograms of Coefficient Estimates with 95% Confidence Intervals', gp=gpar(fontsize=13, fontface = 2))
)




# Combine the two bootstrapped data frames
boot.df_combined <- bind_rows(
  boot.df %>% mutate(group = "Original"),
  boot.df.reduce %>% mutate(group = "Reduced")
)

idx = 1:nrow(boot.df_combined)

# Create a function to plot histograms with quantile and coefficient lines and a title
plot.hist <- function(x, coef.est, coef.est.reduce, varname) {
  group <- ifelse(idx <= nrow(boot.df), "Original", "Reduced")
  p <- ggplot(data.frame(x, group = group), aes(x = x)) + 
    geom_histogram(binwidth = 0.05, alpha = 0.6, aes(fill = group), position = "identity") +
    geom_vline(xintercept = coef.est, color = "red", linetype = "dashed") +
    geom_vline(xintercept = coef.est.reduce, color = "blue", linetype = "dashed") +
    xlab("") + ylab("Frequency") +
    ggtitle(varname) +
    theme(plot.title = element_text(size = 9.5)) +
    scale_fill_manual(values = c("Original" = "red", "Reduced" = "blue")) + theme(legend.position = "none")
  return(p)
}

# Create a list of plots for each column in boot.df with titles
plot.list <- mapply(plot.hist, x = boot.df_combined[,-7], coef.est = coef.estimates, coef.est.reduce = coef.estimates.reduce, varname = names(boot.df_combined)[-7], SIMPLIFY = FALSE)

grid.arrange(grobs = plot.list, ncol = 3, top = textGrob('Overlayed Histograms of Coefficient Estimates, \nBlue is reduced data, Red is imputed data', gp=gpar(fontsize=13, fontface = 2)))



summary(finmod)



quantile(boot.df$povertyPercent, probs = c(0.004166667, 0.9958333))
quantile(boot.df$PctHS18_24, probs = c(0.004166667, 0.9958333))
quantile(boot.df$isSoutheastYes, probs = c(0.004166667, 0.9958333))
quantile(boot.df$PctPublicCoverage, probs = c(0.004166667, 0.9958333))
quantile(boot.df$PctPublicCoverageAlone, probs = c(0.004166667, 0.9958333))
quantile(boot.df$PctUnemployed16_Over, probs = c(0.004166667, 0.9958333))

