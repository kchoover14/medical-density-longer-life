######################## SETUP

library(httr)          # WHO API calls
library(jsonlite)      # parse API JSON response
library(countrycode)   # convert country names to ISO3 codes
library(dplyr)         # data wrangling
library(tidyr)         # drop_na
library(stringr)       # clean character values from API
library(ggplot2)       # static figures
library(plotly)        # interactive choropleth map
library(htmlwidgets)   # save plotly as standalone HTML
library(fitdistrplus)  # distribution fitting and visualization
library(car)           # QQ plot of residuals
library(psych)         # pairs panel for collinearity check
library(effectsize)    # eta-squared effect sizes
library(broom)         # tidy model output for export

######################## DATA ACQUISITION

# Note: The WHO API download and data merging were completed in 2023.
# The merged dataset is saved as finalData-Expectency.csv in the repo.
# Run DATA ACQUISITION and DATA CLEANING only to re-download fresh data.
# To reproduce analysis from the saved dataset, skip to TIME SERIES ANALYSIS.

# Get WHO indicator list (reference only -- not required for analysis)
indicators = GET('https://ghoapi.azureedge.net/api/Indicator')
rawToChar(indicators$content)
indicators = fromJSON(rawToChar(indicators$content))
indicators = indicators[-1]
list2env(indicators, envir = .GlobalEnv)
rm(indicators)

# Download dentist density (HWF_0010: density of dentists per 10,000)
data = GET('https://ghoapi.azureedge.net/api/HWF_0010')
rawToChar(data$content)
data = fromJSON(rawToChar(data$content))
dentists = data$value
dentists = dplyr::select(dentists, SpatialDim, TimeDim, NumericValue)
dentists = rename(dentists, dentists = NumericValue)
dentists$dentists = str_replace_all(dentists$dentists, pattern = " ", repl = "")
dentists$dentists = as.numeric(dentists$dentists)

# Download medic density (HWF_0001: density of medical doctors per 10,000)
data = GET('https://ghoapi.azureedge.net/api/HWF_0001')
rawToChar(data$content)
data = fromJSON(rawToChar(data$content))
medics = data$value
medics = dplyr::select(medics, SpatialDim, TimeDim, NumericValue)
medics = rename(medics, medics = NumericValue)
medics$medics = str_replace_all(medics$medics, pattern = " ", repl = "")
medics$medics = as.numeric(medics$medics)

# Download nurse and midwifery density (HWF_0006: density of nursing and midwifery per 10,000)
data = GET('https://ghoapi.azureedge.net/api/HWF_0006')
rawToChar(data$content)
data = fromJSON(rawToChar(data$content))
nurse = data$value
nurse = dplyr::select(nurse, SpatialDim, TimeDim, NumericValue)
nurse = rename(nurse, nursesMid = NumericValue)
nurse$nursesMid = str_replace_all(nurse$nursesMid, pattern = " ", repl = "")
nurse$nursesMid = as.numeric(nurse$nursesMid)

######################## DATA CLEANING AND MERGING

# Life expectancy from Kaggle (Gochiashvili): kaggle.com/datasets/lashagoch/life-expectancy-who-updated
lifeExp = read.csv('Life-Expectancy-Data-Updated.csv')
lifeExp = rename(lifeExp, lifeExpectancy = Life_expectancy)
lifeExp = rename(lifeExp, year = Year)
lifeExp = dplyr::select(lifeExp, Region, Country, year, lifeExpectancy)
lifeExp$iso3 = countrycode(lifeExp$Country, origin = "country.name", destination = "iso3c")
lifeExp = relocate(lifeExp, iso3, .after = Country)

# Merge WHO density indicators on country ISO3 code and year
expectancy = left_join(dentists, medics, by = c("SpatialDim", "TimeDim"))
expectancy = left_join(expectancy, nurse, by = c("SpatialDim", "TimeDim"))
expectancy = rename(expectancy, iso3 = SpatialDim)
expectancy = rename(expectancy, year = TimeDim)
expectancy = filter(expectancy, year %in% 2000:2015)
expectancy = left_join(expectancy, lifeExp, by = c("iso3", "year"))

# Tidy
rm(lifeExp, dentists, medics, nurse, data)

# Reorder columns
expectancy = relocate(expectancy, Country, .before = iso3)
expectancy = relocate(expectancy, Region, .before = Country)
expectancy = relocate(expectancy, lifeExpectancy, .before = dentists)

# Remove rows with missing life expectancy
expectancy = drop_na(expectancy, lifeExpectancy)

# Remove punctuation from country names
expectancy$Country = gsub('[[:punct:] ]+', ' ', expectancy$Country)

# Save merged dataset
write.csv(expectancy, 'who_density_life_expectancy_merged.csv', row.names = FALSE, quote = FALSE)

######################## TIME SERIES ANALYSIS

# Load pre-merged dataset (start here if skipping DATA ACQUISITION and DATA CLEANING)
expectancy = read.csv('who_density_life_expectancy_merged.csv')

# Summarize median life expectancy by region and year
timeSeries = expectancy |>
  group_by(Region, year) |>
  summarize(lifeExpectancy = median(lifeExpectancy, na.rm = TRUE)) |>
  arrange(desc(Region))

timeSeries = drop_na(timeSeries, lifeExpectancy)

# Line chart: annual life expectancy by region
timeSeries |>
  ggplot(aes(x = year, y = lifeExpectancy, color = Region)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Year",
    y = "Life Expectancy",
    caption = "Source: World Health Organization",
    color = "Country"
  ) +
  facet_wrap(vars(Region), scales = 'free', nrow = 3) +
  scale_color_viridis_d() +
  theme_classic()

ggsave('annual_life_expectancy_by_region.png', width = 10, height = 8, dpi = 300)

rm(timeSeries)

######################## SCATTERPLOTS

# Subset by life expectancy threshold
expectancyPlus = filter(expectancy, lifeExpectancy >= 60)  # model subset
expectancyMinus = filter(expectancy, lifeExpectancy < 60)  # primary aid targets (n=115)

# Countries with life expectancy below 60 -- primary aid targets
unique(expectancyMinus$Country)

# Full dataset scatterplots
ggplot(expectancy, aes(lifeExpectancy, dentists)) +
  geom_point(aes(color = Region)) +
  geom_smooth(method = lm, color = "black") +
  scale_color_viridis_d(alpha = .4) +
  labs(x = "Life Expectancy", y = "Density of Dentists") +
  theme_classic()

ggsave('scatterplot_dentists_all.png', width = 8, height = 6, dpi = 300)

ggplot(expectancy, aes(lifeExpectancy, medics)) +
  geom_point(aes(color = Region)) +
  geom_smooth(method = lm, color = "black") +
  scale_color_viridis_d(alpha = .4) +
  labs(x = "Life Expectancy", y = "Density of Medics") +
  theme_classic()

ggsave('scatterplot_medics_all.png', width = 8, height = 6, dpi = 300)

ggplot(expectancy, aes(lifeExpectancy, nursesMid)) +
  geom_point(aes(color = Region)) +
  geom_smooth(method = lm, color = "black") +
  scale_color_viridis_d(alpha = .4) +
  labs(x = "Life Expectancy", y = "Density of Nurses and Mid-Wives") +
  theme_classic()

ggsave('scatterplot_nurses_all.png', width = 8, height = 6, dpi = 300)

# Scatterplots for life expectancy >= 60 (model subset -- QMD figures)
ggplot(expectancyPlus, aes(lifeExpectancy, dentists)) +
  geom_point(aes(color = Region)) +
  geom_smooth(method = lm, color = "black") +
  scale_color_viridis_d(alpha = .4) +
  labs(x = "Life Expectancy", y = "Density of Dentists") +
  theme_classic()

ggsave('scatterplot_dentists_high_le.png', width = 8, height = 6, dpi = 300)

ggplot(expectancyPlus, aes(lifeExpectancy, medics)) +
  geom_point(aes(color = Region)) +
  geom_smooth(method = lm, color = "black") +
  scale_color_viridis_d(alpha = .4) +
  labs(x = "Life Expectancy", y = "Density of Medics") +
  theme_classic()

ggsave('scatterplot_medics_high_le.png', width = 8, height = 6, dpi = 300)

ggplot(expectancyPlus, aes(lifeExpectancy, nursesMid)) +
  geom_point(aes(color = Region)) +
  geom_smooth(method = lm, color = "black") +
  scale_color_viridis_d(alpha = .4) +
  labs(x = "Life Expectancy", y = "Density of Nurses and Mid-Wives") +
  theme_classic()

ggsave('scatterplot_nurses_high_le.png', width = 8, height = 6, dpi = 300)

# Scatterplots for life expectancy < 60 (aid target subset)
ggplot(expectancyMinus, aes(lifeExpectancy, dentists)) +
  geom_point(aes(color = Region)) +
  geom_smooth(method = lm, color = "black") +
  scale_color_viridis_d(alpha = .4) +
  labs(x = "Life Expectancy", y = "Density of Dentists") +
  theme_classic()

ggsave('scatterplot_dentists_low_le.png', width = 8, height = 6, dpi = 300)

ggplot(expectancyMinus, aes(lifeExpectancy, medics)) +
  geom_point(aes(color = Region)) +
  geom_smooth(method = lm, color = "black") +
  scale_color_viridis_d(alpha = .4) +
  labs(x = "Life Expectancy", y = "Density of Medics") +
  theme_classic()

ggsave('scatterplot_medics_low_le.png', width = 8, height = 6, dpi = 300)

ggplot(expectancyMinus, aes(lifeExpectancy, nursesMid)) +
  geom_point(aes(color = Region)) +
  geom_smooth(method = lm, color = "black") +
  scale_color_viridis_d(alpha = .4) +
  labs(x = "Life Expectancy", y = "Density of Nurses and Mid-Wives") +
  theme_classic()

ggsave('scatterplot_nurses_low_le.png', width = 8, height = 6, dpi = 300)

######################## DISTRIBUTION ANALYSIS

# Life expectancy is negatively skewed -- histogram confirms non-Gaussian distribution
ggplot(expectancy, aes(x = lifeExpectancy)) +
  geom_histogram(aes(y = ..density..), colour = "#005699", fill = "white") +
  geom_density(alpha = .2, fill = "#6799c0") +
  xlab("Life Expectancy") +
  theme_classic()

fit = drop_na(expectancy, lifeExpectancy)

# Describe distribution -- life expectancy follows Weibull distribution
descdist(fit$lifeExpectancy)
fit.weibull = fitdist(fit$lifeExpectancy, "weibull")
plot(fit.weibull, hist = TRUE, demp = TRUE)

# Log transformation increases negative skew -- not suitable
fit$expLog = log(fit$lifeExpectancy)

ggplot(fit, aes(expLog)) +
  geom_histogram(aes(y = ..density..), colour = "#005699", fill = "white") +
  geom_density(alpha = .2, fill = "#6799c0") +
  xlab("Life Expectancy (log)") +
  theme_classic()

# Square root transformation does not correct skew
fit$expSqrt = sqrt(fit$lifeExpectancy)

ggplot(fit, aes(expSqrt)) +
  geom_histogram(aes(y = ..density..), colour = "#005699", fill = "white") +
  geom_density(alpha = .2, fill = "#6799c0") +
  xlab("Life Expectancy (sqrt)") +
  theme_classic()

# Conclusion: untransformed variable used in linear model with acknowledgment
# of departure from normality. Model is for policy inference, not scientific precision.
rm(fit, fit.weibull)

######################## COLINEARITY

# Pairs panel for predictors and outcome -- check relationships and collinearity
pairs.panels(
  dplyr::select(expectancy, lifeExpectancy, dentists, medics, nursesMid),
  lm = TRUE,
  cor = TRUE,
  method = "spearman",
  pch = 21,
  stars = TRUE
)

######################## LINEAR MODEL

# Fit linear regression on countries with life expectancy >= 60
densityModel = lm(lifeExpectancy ~ dentists + medics + nursesMid, data = expectancyPlus)

# Model fit
AIC(densityModel)

# Residuals -- slight patterning expected due to non-normal outcome distribution
res1 = resid(densityModel)
plot(fitted(densityModel), res1)
abline(0, 0)
car::qqPlot(res1)

# Effect sizes
eta_squared(densityModel, partial = FALSE)

# Model summary
summary(densityModel)

# Export tidy model results
tidyDensityModel = tidy(densityModel)
write.csv(tidyDensityModel, 'healthcare_density_regression_results.csv', row.names = FALSE)

######################## CHOROPLETH MAP

dentDensity = expectancy |> dplyr::select(Country, iso3, dentists)
l = list(color = toRGB("grey"), width = 0.5)

dentistMap = plot_geo(dentDensity, height = 500) |>
  add_trace(
    z = ~dentists,
    color = ~dentists,
    colors = 'Blues',
    text = ~Country,
    locations = ~iso3,
    marker = list(line = l)
  ) |>
  colorbar(
    title  = 'Dentist Density',
    thickness = 15,
    len    = 0.6,
    x      = 1.0,
    y      = 0.5
  ) |>
  layout(
    geo    = list(
      showframe    = FALSE,
      showcoastlines = TRUE,
      projection   = list(type = "natural earth")
    ),
    margin = list(l = 0, r = 80, t = 50, b = 20)
  )

htmlwidgets::saveWidget(dentistMap, 'dentist_density_map.html', selfcontained = TRUE)