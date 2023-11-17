library(dplyr)
library(ggplot2)

# Load mortality data
df_mort = read.csv("mortality.txt", header=F)
colnames(df_mort) = "mortality"
df_mort = df_mort %>% mutate(
    age = 0:(nrow(df_mort)-1)
)

# Plot mortality data
p = ggplot(df_mort, aes(x=age, y=mortality)) +
    geom_point() +
    # geom_smooth(method="lm", se=F) +
    labs(x="Age", y="Per-person yearly mortality rate") +
    # log scale on y axis
    scale_y_log10()
ggsave("mortality_by_age.png", p, width=6, height=4)


# Define mean age of PWID based on Meijerink
# Linear model based on two x, y pairs:
# (1972, 20) and (2020, 45)
meanage = function(x) {
    x1 = 1972
    y1 = 20
    x2 = 2020
    y2 = 45
    y = (y2-y1)/(x2-x1)*(x-x1) + y1
}

# Define average mortality output dataframe:
years = 1972:2040
df_avg_mort = data.frame(
    year = years,
    mean_age = meanage(years)
)

# Calculate average mortality rate for each year by assuming a gaussian distribution of age around the mean age and convolving with the mortality rate
df_avg_mort$mean_mortality = 0
df_avg_mort$sd = 5
for (i in 1:nrow(df_avg_mort)) {
    year = df_avg_mort$year[i]
    mean_age = df_avg_mort$mean_age[i]
    df_avg_mort$mean_mortality[i] = sum(df_mort$mortality * dnorm(df_mort$age, mean=mean_age, sd=df_avg_mort$sd))
}

# Save to file
write.csv(df_avg_mort, "mean_age_and_mortality_rate.csv", row.names=F)

# Plot average mortality rate (linear scale)
p = ggplot(df_avg_mort, aes(x=year, y=mean_mortality)) +
    geom_point() +
    # geom_smooth(method="lm", se=F) +
    labs(x="Year", y="Average yearly mortality rate")
ggsave("average_mortality_rate.png", p, width=6, height=4)

# Plot mean age
p = ggplot(df_avg_mort, aes(x=year, y=mean_age)) +
    geom_point() +
    # geom_smooth(method="lm", se=F) +
    labs(x="Year", y="Mean age of PWID")
ggsave("mean_age.png", p, width=6, height=4)


# Plot average mortality rate, log scale
p = ggplot(df_avg_mort, aes(x=year, y=mean_mortality)) +
    geom_point() +
    # geom_smooth(method="lm", se=F) +
    labs(x="Year", y="Average yearly mortality rate") +
    scale_y_log10()
ggsave("average_mortality_rate-log.png", p, width=6, height=4)