library(dplyr)
library(tidyr)
library(Gini)
library(ggplot2)

# Load data
#dfdeaths = read.csv("narkotikautloste_dodsfall_1996-2021.csv", sep = ",", header = TRUE) %>% 
dfdeaths = read.csv("narkotikautloste_dodsfall_1996-2022.csv", sep = ",", header = TRUE) %>% 
    filter(Bofylke != "Ukjent") %>% # No counts in Ukjent
    rename(
        year = Dødsår,
        county_name = Bofylke,
        number_of_deaths = Antall
)
# Fill forward the column "year" because the data only reports year for the first county in each year
dfdeaths = dfdeaths %>% fill(year, .direction = "down")

# Change "-" to 0 in number_of_deaths
dfdeaths$number_of_deaths[dfdeaths$number_of_deaths == "-"] = 0
# Change number_of_deaths to numeric
dfdeaths$number_of_deaths = as.numeric(dfdeaths$number_of_deaths)

# Get the population of each county by year
dfpop = fhidata::norway_population_by_age_b2020 %>% 
    filter(year %in% 1996:2022 & granularity_geo == "county") %>% 
    rename(county_code = location_code) %>%
    group_by(year, county_code) %>%
    summarise(pop = sum(pop)) %>%
    ungroup()

# The population data does not exist before 2005, so we assume that the population is the same as in 2005 for all years 1996-2004
for (i in 1996:2004) {
    dfpop_addyear = dfpop %>% filter(year == 2005) %>% mutate(year = i)
    dfpop = rbind(dfpop, dfpop_addyear)
}
dfpop = dfpop %>% arrange(year)

# # Write dfpop to csv and read it back in to not be reliant on fhidata package
write.csv(dfpop, "dfpop.csv", row.names = FALSE)
dfpop = read.csv("dfpop.csv", sep = ",", header = TRUE)

# Map county name to county code:
df_mapping = data.frame(
    county_code = paste0("county", c("03", "11", "15", "18", "30", "34", "38", "42", "46", "50", "54")),
    county_name = c("Oslo", "Rogaland", "Møre og Romsdal", "Nordland", "Viken", "Innlandet", "Vestfold og Telemark", "Agder", "Vestland", "Trøndelag", "Troms og Finnmark")
)
dfpop = dfpop %>% left_join(df_mapping, by = "county_code")

# Merge df and dfpop
df = dfdeaths %>% left_join(dfpop, by = c("year", "county_name"))

# Calculate the Gini coefficient for each year
df_gini = df %>% 
    group_by(year) %>% 
    summarise(gini = gini(pop, number_of_deaths))


# Make a 5-year smoothed average of the Gini coefficient
df_gini$gini_smoothed = zoo::rollmean(df_gini$gini, 5, fill = NA)


# Plot the Gini coefficient and the smoothed version as function of year
df_gini %>% 
    ggplot(aes(x = year, y = gini)) +
    geom_line() +
    geom_point() +
    geom_line(aes(y = gini_smoothed), color = "red") +
    ylim(0, 0.5) +
    labs(
        x = "Year",
        y = "Gini coefficient / 5-year rolling mean",
        title = "Gini coefficient for drug-related deaths in Norway 1996-2021"#,
        # subtitle = "Data source: https://www.fhi.no/hn/helseregistre-og-registre/dodsarsaksregisteret/"
    )
# Save plot
ggsave("gini_coefficient_drug_deaths.png", width = 10, height = 6, dpi = 300)


# Add an artificial GINI coefficient of 1 at 1972
df_gini = df_gini %>% 
    add_row(year = 1972, gini = 1, gini_smoothed = 1) %>% 
    arrange(year)

# Add NA for years 1973-1995
for (i in 1973:1995) {
    df_gini = df_gini %>% add_row(year = i, gini = NA, gini_smoothed = NA)
}


# Add years 2023-2030
for (i in 2023:2030) {
    df_gini = df_gini %>% add_row(year = i, gini = NA, gini_smoothed = NA)
}

df_gini = df_gini %>% arrange(year)


# # Fit a second-degree polynomial to the GINI coefficient
# fit = lm(gini ~ poly(year, 2, raw = TRUE), data = df_gini)
# summary(fit)
# Fit an exponential function to the GINI coefficient
fit = lm(log(gini) ~ year, data = df_gini)
summary(fit)
# Add the fit as a new column
df_gini$gini_fit_exp = exp(predict(fit, df_gini))

# Fit a linear function to the GINI coefficient
fit = lm(gini ~ year, data = df_gini)
summary(fit)
# Add the fit as a new column
df_gini$gini_fit_linear = predict(fit, df_gini)

# Plot the "gini" and "gini_fit" columns as function of year
p = ggplot(data = df_gini %>% filter(year > min(year)), aes(x = year, y = gini, color="Data")) +
    geom_line() +
    geom_point() +
    geom_line(data = df_gini, aes(y = gini_fit_exp, color="Fit")) +
    # geom_line(aes(y = gini_fit_linear, color="lin fit")) +
    ylim(0, 1.05) +
    labs(
        x = "Year",
        y = "Gini coefficient",
        # title = "Gini coefficient for drug-related deaths in Norway 1996-2022",
        # subtitle = "Data source: https://www.fhi.no/hn/helseregistre-og-registre/dodsarsaksregisteret/",
        color = ""
    ) +
    theme(legend.position = "bottom")
# Save plot
ggsave("gini_coefficient_drug_deaths_fit.png", width = 10, height = 6, dpi = 300)


# Write the df_gini dataframe to csv
write.csv(df_gini, "df_gini.csv", row.names = TRUE)
