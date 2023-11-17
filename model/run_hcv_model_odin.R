rm(list=ls(all=TRUE))
library(tidyr)
library(ggplot2)
library(dplyr)
library(odin)
library(dust)
library(odin.dust)
library(mcstate)
library(reshape2)
library(zoo)

library(GGally)
library(coda)

# Set up scenario vector in case we want to run multiple scenarios
# -- they can be defined to be different in the code below
scenarios = c(
    "baseline"
)


for (is in 1:length(scenarios)) {
    scenario = scenarios[is]

    # Print scenario
    print(paste0("Running scenario: ", scenario))

    # Set up output directory
    output_dir = paste0("output_and_plots/scenario_", is, "-", scenario, "/")
    # Make directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Set model time step (in years)
    dt = 0.25 # 1 

    # Set MCMC and particle filter parameters
    n_particles = 1e4
    n_burnin = 2000
    n_steps = n_burnin + 5000
    n_chains = 2
    n_thin = 5
    n_threads = 40 # Number of CPU threads (OMP)
    

    ## ===== Load data =====
    df_pwid = read.csv("../data/number_of_PWID-2000_peak_modified.csv")
    df_pwid = df_pwid %>%
        mutate(
            N_active = as.integer(round(N_active)), # Integer for poisson likelihood
            N_ex_wnr = as.integer(round(N_ex_wnr)), # Integer for poisson likelihood
            N_ex_wr = as.integer(round(N_ex_wr)), # Integer for poisson likelihood
            N_active_lower = as.integer(round(N_active_lower)), # Integer for poisson likelihood
            N_active_upper = as.integer(round(N_active_upper)) # Integer for poisson likelihood
        )



    df_hcv = read.csv("../data/pwid_HCV_prevalence-all_surveys_combined.csv") %>%
        rename(N_sampled_hcv = N_sampled_rna) # Renaming for convenience -- hcv here means they're RNA positive for HCV, as opposed to antibody positive.

    file_df_treatment = "../data/treatments-93percent.csv"
    df_treatment = read.csv(file_df_treatment) %>% 
        mutate(
            year_int = year - 1972 + 1, # Note: year 1 is 1972.
            n_treated_total = as.integer(n_treated_total) # Needed for poisson likelihood
        )

    p = ggplot(data = df_treatment %>% filter(!is.na(n_treated_total))) +
        geom_point(aes(x=year, y=n_treated_total)) +
        xlab("Year") + ylab("Number of treatments")
    ggsave(plot=p, file=paste0("output_and_plots/scenario_", is, "-", scenario, "/data_treatments.png"))



    # NB! Assuming ("impute") a number where N_sampled_hcv is missing:
    df_hcv_imputed = df_hcv
    df_hcv_imputed$N_sampled_hcv[is.na(df_hcv_imputed$N_sampled_hcv)] = 200 
    df_hcv = df_hcv_imputed
                        
    # Calculate clopper-pearson uncertainty bands for rna and antibody positives
    df_hcv_errband_rna = binom::binom.confint(df_hcv$fraction_positive_hcv*df_hcv$N_sampled_hcv, df_hcv$N_sampled_hcv, conf.level = 0.95, method="exact")
    df_hcv = cbind(df_hcv, df_hcv_errband_rna)
    df_hcv_errband_ag = cbind(
        df_hcv[!is.na(df_hcv$fraction_antibody_positive),] %>% select(year_int, city), 
        binom::binom.confint(df_hcv[!is.na(df_hcv$fraction_antibody_positive), ]$fraction_antibody_positive*df_hcv[!is.na(df_hcv$fraction_antibody_positive), ]$N_sampled_antibody, df_hcv[!is.na(df_hcv$fraction_antibody_positive), ]$N_sampled_antibody, conf.level = 0.95, method="exact")
    )
    colnames(df_hcv_errband_ag)[1] = "year_int"
    df_hcv = left_join(
        df_hcv, 
        df_hcv_errband_ag %>% select(year_int, city, mean, lower, upper) %>% rename(mean_ag = mean, lower_ag = lower, upper_ag = upper),
        by = c("year_int", "city")
    )

    # Write df_hcv to file to use in plots_report.R
    write.csv(df_hcv, file="../data/pwid_HCV_prevalence-all_surveys_combined-with_uncertainty_columns_for_plotting.csv", row.names=F)

    # Add "standard deviation" column which is mean diff of clopper-pearson upper/lower
    df_hcv = df_hcv %>% mutate(sd_hcv_pos = 0.5*(abs(upper - mean) + abs(lower - mean)))

    # Plot HCV PWID prevalence data (including imputations and uncertainty)
    p = ggplot(data = df_hcv) +
        geom_point(aes(x=year_int+1971, y=mean, col=city)) +
        geom_errorbar(aes(x=year_int+1971, ymin=lower, ymax=upper, col=city)) + 
        ylab("Fraction HCV-RNA positive") +
        theme(legend.position="bottom") +
        xlab("Year")
    ggsave(plot = p, file=paste0("output_and_plots/scenario_", is, "-", scenario, "/data_HCV_prevalence.png"))

    # The df_hcv has multiple entries for the same year. This needs to be split to two columns for the particle filter data
    # Add a column that maps city to "Oslo" or "other"
    df_hcv = df_hcv %>% mutate(
        city_oslo = ifelse(city == "Oslo", "Oslo", "other")
    )

    # Select columns mean, lower and upper to pivot wider on city_oslo
    df_hcv_wide = df_hcv %>% 
        select(
            year_int, city_oslo, fraction_positive_hcv, N_sampled_hcv
        ) %>% 
        pivot_wider(
            names_from = city_oslo, 
            values_from = c(fraction_positive_hcv, N_sampled_hcv)
        )

    # Combine
    min_year_data = min(c(df_hcv$year_int, df_pwid$year_int))
    max_year_data = max(c(df_hcv$year_int, df_pwid$year_int))
    df_hcv_with_NA = data.frame(
            year_int = min_year_data:max_year_data
        ) %>% 
        left_join(
            # df_hcv, by="year_int"
            df_hcv_wide, by="year_int"
        )
    df_pwid_with_NA = data.frame(
            year_int = min_year_data:max_year_data
        ) %>%
        left_join(
            df_pwid, by="year_int"
        )
    df_treatment_with_NA = data.frame(
            year_int = min_year_data:max_year_data
        ) %>%
        left_join(
            df_treatment, by="year_int"
        )



    # Combine everything to one dataframe to go into likelihood function
    df_data = df_pwid_with_NA %>% 
        inner_join(df_hcv_with_NA, by="year_int") %>% 
        inner_join(df_treatment_with_NA, by="year_int")

    # Convert to mcstate data object
    data = mcstate::particle_filter_data(
        data = df_data,# %>% select(!year),
        time = "year_int",
        initial_time = 0,
        rate = 1/dt
    )


    ### Data to be input directly into model
    year_vector = seq(from=min_year_data, to=max_year_data, by=1)+1971

    # Treatment success probability
    # Make the vector of p_treatment_success and corresponding time_vector for input and interpolation
    p_treatment_success = df_data$percentage_successful_treatment_by_year/100
    p_treatment_success[is.na(p_treatment_success)] = min(p_treatment_success[!is.na(p_treatment_success)])
    time_step_vector = seq(from=min_year_data,to=max_year_data,by=dt)
    p_treatment_success_interpolated = approx(
        x=df_data$year_int, 
        y=p_treatment_success, 
        xout=time_step_vector,
        method="constant" 
    )

    # Treatment rates for the pre-data period
    rate_treatment_pre_data = ifelse(year_vector >= min(df_treatment$year), 0.02, 0) # 0.02 is the assumed treatment rate from 1990-2004 from whence we have data
    rate_treatment_input_interpolated = approx(
        x=year_vector-1971, 
        y=rate_treatment_pre_data, 
        xout=time_step_vector,
        method="constant" 
    )


    # Special treatment rates for active PWID
    df_treatment_active = data.frame(
        year = 1972:2022
    )
    hcv_treatment_rate_active_pwid_paper = c(0.34, 1.3, 1.4, 2.1, 2.2, 3.7, 4.6, 11.2, 24.1, 22.1) / 100
    # Custom make the rate_treatment_active 
    # - 0 from 1972 til 1989 inclusive
    # - exponential increase from 1990 to 2003 inclusive, with value 0.00034 in 1990 and 0.0034 in 2003
    # - constant from 2004 to 2009 inclusive, with value 0.0034
    # - data from vector hcv_treatment_rate_active_pwid_paper from 2010 to 2019
    # - constant from 2020 to 2022 inclusive, with value equal to 2019
    df_treatment_active$rate_treatment_active = ifelse(
        df_treatment_active$year < 1990, 
        0, 
        ifelse(
            df_treatment_active$year < 2004, 
            0.00034 * exp((df_treatment_active$year - 1990) * log(0.0034/0.00034)/(2004-1990)), 
            0.0034
        )
    )
    for (i in 1:length(hcv_treatment_rate_active_pwid_paper)) {
        df_treatment_active$rate_treatment_active[df_treatment_active$year == (2009 + i)] = hcv_treatment_rate_active_pwid_paper[i]
    }
    df_treatment_active$rate_treatment_active[df_treatment_active$year > 2019] = hcv_treatment_rate_active_pwid_paper[length(hcv_treatment_rate_active_pwid_paper)]


    # Interpolate it to the time step resolution
    rate_treatment_active_input_interpolated = approx(
        x=df_treatment_active$year-1971, 
        y=df_treatment_active$rate_treatment_active, 
        xout=time_step_vector,
        method="linear"
    )


    # NSP coverage: Define some coverage fractions, and interpolate 
    # Define the coverage fractions
    df_nsp_points = data.frame(
        year =     c(1972, 1988, 2016, 2017, 2018, 2019, 2020, 2021, 2022),
        coverage = c(   0,    0, 0.77, 0.80, 0.87, 0.92, 0.92, 0.92, 0.92)
    )
    # Interpolate the coverage fractions to get a continuous curve of coverage 1988-2022
    df_nsp = data.frame(
        year = seq(from=min(df_nsp_points$year), to=max(df_nsp_points$year), by=1)
    ) %>%
        mutate(
            year_int = year - 1971,
            # nsp_coverage = spline(
            #     x=df_nsp_points$year, 
            #     y=df_nsp_points$coverage, 
            #     xout=min(df_nsp_points$year):max(df_nsp_points$year),
            #     method="natural" # TODO does it make a difference if constant or linear?
            # )$y
            nsp_coverage = approx(
                x=df_nsp_points$year, 
                y=df_nsp_points$coverage, 
                xout=min(df_nsp_points$year):max(df_nsp_points$year),
                method="linear" # TODO does it make a difference if constant or linear?
            )$y
        )

    # Plot the NSP coverage
    # Y axis = "Percentage coverage of needle and syringe programmes", X axis = "Year"
    p = ggplot(df_nsp, aes(x=year, y=nsp_coverage)) + geom_line() +
        xlab("Year") + ylab("Percentage coverage of needle and syringe programmes")
    ggsave(filename=paste0("output_and_plots/scenario_", is, "-", scenario, "/nsp_coverage.png"), plot=p, width=10, height=5)


    # Calculate the reduction in risk on the overall population from NSP coverage
    max_rel_reduction_nsp = 0.56 # Updated 202308 to literature value 0.56 # Reduced to 0.3 20230329 because 0.5 and high NSP coverage was giving unrealistically low total incidence in 2022 0.5 # Updated to 0.5 in agreement with Rob  0.316
    df_nsp$risk_multiplier_nsp = 1 - df_nsp$nsp_coverage*max_rel_reduction_nsp


    ## LAR (legemiddelassistert rehabilitering = OST) coverage
    # Read LAR numbers
    df_lar = read.csv("../data/persons_in_LAR.csv")
    # Extend df_lar to year_vector
    df_lar_full = data.frame(
        year = year_vector
    )
    df_lar = df_lar_full %>% left_join(df_lar, by="year")
    # Plot the LAR coverage
    p = ggplot(df_lar, aes(x=year, y=persons_in_LAR)) + geom_line() +
        xlab("Year") + ylab("Number of patients")
    ggsave(filename=paste0("output_and_plots/scenario_", is, "-", scenario, "/persons_in_LAR.png"), plot=p, width=10, height=5)

    # Join active PWID on LAR numbers
    df_lar = df_lar %>% left_join(df_pwid %>% select(year, N_active), by="year")
    # Fill in NAs in N_active with previous value
    df_lar = df_lar %>% mutate(N_active = ifelse(is.na(N_active), lag(N_active), N_active))

    # Calculate fraction of active PWID on LAR by assuming 20 percent of LAR patients are active PWID
    active_PWID_LAR_correction_factor = 0.2
    df_lar = df_lar %>% mutate(
        lar_coverage = ifelse(persons_in_LAR < 1, 0, persons_in_LAR/N_active*active_PWID_LAR_correction_factor)
    )
    # Plot the LAR coverage
    p = ggplot(df_lar, aes(x=year, y=lar_coverage)) + geom_line()
    ggsave(filename=paste0("output_and_plots/scenario_", is, "-", scenario, "/lar_coverage.png"), plot=p, width=10, height=5)

    # Calculate the reduction in risk on the overall population from LAR coverage
    max_rel_reduction_lar = 0.5
    df_lar$risk_multiplier_lar = 1 - df_lar$lar_coverage*max_rel_reduction_lar

    # Combined reduction in risk from LAR and NSP
    df_lar_nsp = df_lar %>% left_join(df_nsp %>% select(year, risk_multiplier_nsp), by="year")
    df_lar_nsp = df_lar_nsp %>% mutate(
        risk_multiplier = risk_multiplier_nsp*risk_multiplier_lar # The reduction in risk from the combined LAR and NSP coverage is 1 - risk_multiplier
    )

    # Plot the risk reduction from each individual measure and combined
    df_risk_reductions = df_lar_nsp %>% select(year, risk_multiplier_nsp, risk_multiplier_lar, risk_multiplier) %>%
        mutate(
            risk_reduction_nsp = 1-risk_multiplier_nsp,
            risk_reduction_lar = 1-risk_multiplier_lar,
            risk_reduction_combined = 1-risk_multiplier
        )
    df_risk_reductions_plot = df_risk_reductions %>% select(c(year, risk_reduction_nsp, risk_reduction_lar, risk_reduction_combined)) %>% pivot_longer(cols = c(risk_reduction_nsp, risk_reduction_lar, risk_reduction_combined), names_to = "measure", values_to = "risk_reduction")
    # Rename measure column for plotting:
    # risk_reduction_nsp = Needle and syringe programmes
    # risk_reduction_lar = Opioid substitution therapy
    # risk_reduction_combined = Combined
    df_risk_reductions_plot$measure[df_risk_reductions_plot$measure == "risk_reduction_nsp"] = "Needle and syringe programmes"
    df_risk_reductions_plot$measure[df_risk_reductions_plot$measure == "risk_reduction_lar"] = "Opioid substitution therapy"
    df_risk_reductions_plot$measure[df_risk_reductions_plot$measure == "risk_reduction_combined"] = "Combined"
    # Order the measure column
    df_risk_reductions_plot$measure = factor(df_risk_reductions_plot$measure, levels=c("Needle and syringe programmes", "Opioid substitution therapy", "Combined"))

    # Plot the risk reduction from each individual measure and combined
    p = ggplot(df_risk_reductions_plot, aes(x=year, y=risk_reduction, color=measure)) + geom_line() +
        xlab("Year") + ylab("Reduction in the modelled force of infection") +
        # Put legend below and remove legend title
        theme(legend.position="bottom") +
        guides(color=guide_legend(title=NULL))
    ggsave(filename=paste0("output_and_plots/scenario_", is, "-", scenario, "/risk_reduction.png"), plot=p, width=10, height=5)

    # Interpolate the LAR and NSP coverage
    df_lar_nsp = df_lar_nsp %>% mutate(
        year_int = year - 1971
    )
    df_lar_nsp_interpolated = approx(
        x=df_lar_nsp$year_int, 
        y=df_lar_nsp$risk_multiplier, 
        xout=time_step_vector,
        method="linear"
    )


    # ## GINI coefficient drug deaths
    df_gini = read.csv("../data/gini/df_gini.csv") %>%
        select(year, gini_fit_exp) %>% 
        rename(gini_drug_deaths = gini_fit_exp) %>%
        mutate(
            year_int = year - 1971
        )

    # Plot the GINI coefficient
    p = ggplot(df_gini, aes(x=year, y=gini_drug_deaths)) + geom_line()
    ggsave(filename=paste0("output_and_plots/scenario_", is, "-", scenario, "/gini_drug_deaths.png"), plot=p, width=10, height=5)

    df_gini_interpolated = approx(
        x=df_gini$year_int,
        y=df_gini$gini_drug_deaths, 
        xout=time_step_vector,
        method="linear"
    )

    # Estimate immigration using detailed SSB data:
    df_immigrated_by_country = read.csv("../data/immigrants/ssb_immigrant_prevalence-for_model.csv") %>%
        rename(Country = Country_english) %>%
        mutate(
            Country = ifelse(Country == "Other", "Total", Country) # TODO this is a hack that should be fixed in the data by adding a region mapping
        )

    # Load HCV prevalence by country
    df_hcvprev_country = read.csv(paste0("../data/immigrants/hcv_prevalence_by_country-for_model-mock_values.csv")) %>%
        mutate(
            hcv_country_prevalence = Prevalence/100 # Convert to fraction
        )

    # Left join hcv prevalence on immigration data by country and multiply
    # Change "Other" to "Total" in Country_english
    df_immigrated_by_country = df_immigrated_by_country %>%
        left_join(df_hcvprev_country, by=c("Country"="Country")) %>%
            mutate(
            inc_immigrants_hcvpos = (hcv_country_prevalence * incidence), # Can be positive or negative, includes removal by deaths
            prev_immigrants_hcvpos = (hcv_country_prevalence * prevalence)
        )

    # Sum countries that are part of other countries (e.g. USSR + Russia)
    df_immigrated_by_country = df_immigrated_by_country %>%
        group_by(
            Country,
            Year
        ) %>%
        summarise(
            incidence = sum(incidence),
            prevalence = sum(prevalence),
            inc_immigrants_hcvpos = sum(inc_immigrants_hcvpos),
            prev_immigrants_hcvpos = sum(prev_immigrants_hcvpos)
        ) %>% 
        ungroup()

    # ==== SOME PLOTTING ====
    ## Plot yearly immigration by country
    p = ggplot(
            df_immigrated_by_country,
            aes(
                x = Year,
                y = inc_immigrants_hcvpos,
                fill = Country
            )
        ) +
        geom_bar(position="stack", stat="identity") +
        theme_bw(
            base_size = 10
        ) +
        theme(legend.position="bottom") +
        guides(fill=guide_legend(ncol=8))
    ggsave(plot=p, file=paste0("output_and_plots/scenario_", is, "-", scenario, "/immigrants_hcv_pos_by_country-all_countries.png"), width=20, height=20)

    ## Plot yearly immigration by country, aggregated to top N_countries plus rest
    N_countries = 12
    df_immigrated_totbycountry = df_immigrated_by_country %>%
        group_by(Country) %>%
        summarise(
            inc_immigrants_hcvpos = sum(inc_immigrants_hcvpos),
            incidence = sum(incidence)
            ) %>%
        ungroup() %>%
        arrange(desc(inc_immigrants_hcvpos))
    countries_notagg = df_immigrated_totbycountry$Country[1:N_countries]
    countries_notagg = countries_notagg[countries_notagg != "Total"]
    df_immigrated_by_country$Country_aggregate = ifelse(df_immigrated_by_country$Country %in% countries_notagg, df_immigrated_by_country$Country, "Other")

    # Bar plot stacked
    p = ggplot(
            df_immigrated_by_country %>% group_by(Year, Country_aggregate) %>% summarise(inc_immigrants_hcvpos=sum(inc_immigrants_hcvpos)) %>% ungroup() %>%
                mutate(Country_aggregate = factor(Country_aggregate, levels=c(sort(unique(Country_aggregate[Country_aggregate!="Other"])), "Other"))),
            aes(
                x = Year,
                y = inc_immigrants_hcvpos,
                fill = Country_aggregate
            )
        ) +
        geom_bar(position="stack", stat="identity") +
        theme_bw(
            base_size = 10
        ) +
        labs(
            x = "Year", 
            y = "Estimated net immigration by persons with chronic Hepatitis C"
        ) +
        scale_fill_manual(
            # name = "Land",
            name = "Country",
            # values = c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695','grey')
            values = c('#a50026','dark green','#f46d43','#fdae61','light green','#ffffbf','purple','#abd9e9','#74add1','#4575b4','black',"yellow", 'grey') 
        ) +
        theme(legend.position="bottom")
    ggsave(plot=p, file=paste0("output_and_plots/scenario_", is, "-", scenario, "/immigrants_hcv_pos_by_country-top_countries_and_rest-bar_plot_stacked.png"), width=8, height=6)
    # ==== END PLOTTING ====

    # Sum all countries
    df_immigrated = df_immigrated_by_country %>% 
        rename(year = Year) %>%
        group_by(year) %>%
        summarise(
            inc_immigrants = sum(incidence),
            inc_immigrants_hcvpos = sum(inc_immigrants_hcvpos)
        )

    # Distribute the first rows that apply to multiple years:
    df_immigrated_1970 = df_immigrated[1, ]
    df_immigrated_1980 = df_immigrated[2, ]
    df_immigrated_1986 = df_immigrated[3, ]
    df_immigrated_1987_on = df_immigrated[4:nrow(df_immigrated), ]

    # 1980 should be distributed on 1971-1980:
    df_immigrated_1971_1980 = data.frame(
        year = 1971:1980
    ) %>% mutate(
        inc_immigrants = df_immigrated_1980$inc_immigrants[1] / length(year),
        inc_immigrants_hcvpos = df_immigrated_1980$inc_immigrants_hcvpos[1] / length(year)
    )

    # 1986 should be distributed on 1981-1986:
    df_immigrated_1981_1986 = data.frame(
        year = 1981:1986
    ) %>% mutate(
        inc_immigrants = df_immigrated_1986$inc_immigrants[1] / length(year),
        inc_immigrants_hcvpos = df_immigrated_1986$inc_immigrants_hcvpos[1] / length(year)
    )

    # Rbind all to get a dataframe with all years, which we then will cut at 1972
    df_immigrated_allyr = rbind(
        df_immigrated_1970,
        df_immigrated_1971_1980,
        df_immigrated_1981_1986,
        df_immigrated_1987_on %>% select(year, inc_immigrants, inc_immigrants_hcvpos)
    )

    # Make the initial condition consisting of all counts up to 1972 inclusive
    n_immigrated_hcvpos_1972 = sum(df_immigrated_allyr$inc_immigrants_hcvpos[df_immigrated_allyr$year < 1972])
    df_immigrated_post1972 = df_immigrated_allyr %>% filter(year >= 1972)

    # Divide the post-1972 counts on time steps
    df_immigrated_input = df_immigrated_post1972 %>%
        slice(rep(1:n(), each = as.integer(1/dt))) %>% # Repeat each row 1/dt times
        mutate(
            inc_immigrants = inc_immigrants * dt, # Divide counts equally on timesteps
            inc_immigrants_hcvpos = inc_immigrants_hcvpos * dt # Divide counts equally on timesteps
        )

    # Verify preservation of total counts
    assertthat::assert_that(
        abs((n_immigrated_hcvpos_1972 + sum(df_immigrated_input$inc_immigrants_hcvpos)) - sum(df_immigrated$inc_immigrants_hcvpos)) < 1e-3,
        msg = "In preprocessing, the number of immigrants is not the same after distributing on all years"
    )


    ## Load mortality data
    # df_mortality = read.csv("../data/mortality/mean_age_and_mortality_rate-sd10.csv") %>%
    filename_df_mortality = "../data/mortality/mean_age_and_mortality_rate-sd5.csv"
    df_mortality = read.csv(filename_df_mortality) %>%
        mutate(
            year_int = year-1971
        )
    # Interpolate mortality rate
    mortality_rate_interpolated = approx(
        x = df_mortality$year_int,
        y = df_mortality$mean_mortality,
        xout = time_step_vector,
        method = "linear"
    )



    # ==== Define data comparison functions (likelihood) ====
    ll_pois = function(obs, model ) { # Follwing example at the end of https://mrc-ide.github.io/mcstate/articles/sir_models.html
        # exp_noise = 1e6 # add some small noise to avoid zero counts
        if (is.na(obs)) {
            ll_obs = numeric(length(model))
        } else {
            lambda = model + rexp(n = length(model), rate=1e6)
            # print(paste0("ll_pois: obs = ", obs, ", model = ", model, ", lambda = ", lambda))
            ll_obs = dpois(x = obs, lambda = lambda, log = TRUE)
        }
        ll_obs
    }

    ll_binom = function(x, size, p) {
        # Wrapper function for binomial distribution to take into account possibility of NA values
        # x is the observed number of hcv positives,
        # size is the sample size in hcv prevalence survey
        # p is the modelled fraction of positives

        exp_noise = 1e6 # add some small noise to avoid zero values for p parameter, which seems to trip the dbinom() call

        if (is.na(x)) {
            ll_obs = numeric(length(p))
        } else {
            # Had issue: sometimes there's an "if (accept) { : missing value where TRUE/FALSE is needed" error raised inside the dbinom() call. We seem to be sending in some invalid vector.
            # This seems to have been solved by adding the small noise component to p!
            ll_obs = dbinom(x=x, size=size, p=p+rexp(n = length(p), rate=exp_noise), log=TRUE)
        }
        ll_obs
    }

    # Define total likelihood function, calling it n_compare for some reason
    n_compare = function(state, observed, pars=NULL) {
        # Print model state input (debug)
        # print(paste0("dim(state) = ", dim(state)))

        # Variables to inform PWID population size:
        N_active_modelled = state[1, , drop=TRUE] + state[2, , drop=TRUE] + state[3, , drop=TRUE] + state[4, , drop=TRUE] # Adding together AN+AA+AC+AR. (Use model$info() to get index of state variables)
        N_ex_wr_modelled = state[5, , drop=TRUE] + state[6, , drop=TRUE] + state[7, , drop=TRUE] + state[8, , drop=TRUE] # Adding together EN+EA+EC+ER.
        N_ex_wnr_modelled  = state[9, , drop=TRUE] + state[10, , drop=TRUE] + state[11, , drop=TRUE] + state[12, , drop=TRUE] # Adding together QN+QA+QC+QR.
        N_active_observed = observed$N_active
        N_ex_wnr_observed = observed$N_ex_wnr
        N_ex_wr_observed = observed$N_ex_wr
        # Combine to poisson likelihood
        ll_active = ll_pois(N_active_observed, N_active_modelled) 
        ll_ex_wnr = ll_pois(N_ex_wnr_observed, N_ex_wnr_modelled)
        ll_ex_wr  = ll_pois(N_ex_wr_observed, N_ex_wr_modelled)
        
        # Variables to inform HCV positive population:
        N_active_hpos_modelled = state[2, , drop=TRUE] + state[3, , drop=TRUE] # Using only active PWID to calculate share, compare with N_active_modelled
        N_sampled_hcv_Oslo = observed$N_sampled_hcv_Oslo # Total number of people sampled in survey
        frac_pos_observed_Oslo = observed$fraction_positive_hcv_Oslo # Fraction positive in sample
        N_sampled_hcv_other = observed$N_sampled_hcv_other # Total number of people sampled in survey
        frac_pos_observed_other = observed$fraction_positive_hcv_other # Fraction positive in sample
        # Compare in binomial likelihood
        frac_pos_modelled = N_active_hpos_modelled / N_active_modelled                      
        N_active_hpos_observed_Oslo = as.integer(frac_pos_observed_Oslo * N_sampled_hcv_Oslo)
        N_active_hpos_observed_other = as.integer(frac_pos_observed_other * N_sampled_hcv_other)
        ll_hcv_Oslo = ll_binom(x=N_active_hpos_observed_Oslo, size=N_sampled_hcv_Oslo, p=frac_pos_modelled) # Assuming frac_pos_modelled is uncertainty-free due to large sample
        ll_hcv_other = ll_binom(x=N_active_hpos_observed_other, size=N_sampled_hcv_other, p=frac_pos_modelled) # Assuming frac_pos_modelled is uncertainty-free due to large sample
        
        # Treatments
        incidence_treated_observed = observed$n_treated_total
        incidence_treated_modelled = state[41, , drop=TRUE] + state[42, , drop=TRUE] + state[43, , drop=TRUE]
        ll_treatment = ll_pois(incidence_treated_observed, incidence_treated_modelled)

        # Add together and return total loglikelihood
        loglik = ll_active + ll_ex_wnr + ll_ex_wr + ll_hcv_Oslo + ll_hcv_other + ll_treatment

        return(loglik)
    }


    ### Compile model
    gen = odin.dust::odin_dust("hcv_model_odin.R", workdir=paste0("dust_model-", scenario))
    # gen = odin.dust::odin_dust("hcv_model_odin.R", workdir="dust", gpu=TRUE) # Tested runnig with GPU, but it's not faster


    # ==== Define model parameters ====
    # Set up some parameters that inform the model
    rate_IC_AC = 0.06
    spontaneous_clearance_probability = 0.26
    gini_drug_deaths_interpolated = df_gini_interpolated$y
    OD_mortality_rate = 0.03
    duration_acute = 0.5 

    # The total list of model parameters. Some of these are overridden by MCMC in the below transform() function
    model_params = list(
            rate_debut_ini = 300, # Larger rate for poisson
            rate_quitting_wr_ini = 0.15,   # Locked to values from the mortality multiplier model
            rate_quitting_wnr_ini = 0.025, 
            rate_relapsing_ini = 0.15,     
            sd_rw_debut = 0.1, 
            beta_ini = 0.45,
            AN_ini = 378,
            AA_ini = 150,
            AC_ini = 0,
            AR_ini = 0,
            EA_ini = 0,
            EN_ini = 35,
            EC_ini = 0,
            ER_ini = 0,
            QA_ini = 0,
            QN_ini = 0,
            QC_ini = 0,
            QR_ini = 0,
            IC_ini = n_immigrated_hcvpos_1972,
            IR_ini = 0,
            ID_ini = 0,
            GC_ini = 0,
            GR_ini = 0,
            GD_ini = 0,
            dt = dt,
            rate_treatment_rw_ini = 0.01, 
            rate_treatment_input_step = rate_treatment_input_interpolated$y,
            time_start_rw_treatment = 1990 - 1972,
            sd_rtreat = 0.35,
            rate_treatment_max = 2.0, # Absolute maximum rate of treatment initiation, used to set upper limit on random walk through a logit transform
            rate_treatment_active_input_step = rate_treatment_active_input_interpolated$y, 
            p_treatment_succcess_step = p_treatment_success_interpolated$y,
            reduction_factor_treatment_active = 1, 
            rate_death_active_step = mortality_rate_interpolated$y + OD_mortality_rate,
            rate_death_ex_step = mortality_rate_interpolated$y, 
            risk_multiplier_lar_nsp_step = df_lar_nsp_interpolated$y, # The combined effect of NSP and LAR on the risk of HCV transmission
            gini_drug_deaths_step = gini_drug_deaths_interpolated,
            rate_death_imm = 0, # Deaths are handled by the incidence directly
            rate_IC_AC = rate_IC_AC,
            immigrated_hcvpos_step = df_immigrated_input$inc_immigrants_hcvpos,
            rate_death_genpop = 0.001,
            genpop_hcvpos_step = rep(0, length(df_immigrated_input$inc_immigrants_hcvpos)),
            time_forecast = max(df_data$year_int), # No more random walk after this time
            rate_recovery_ex_ini = spontaneous_clearance_probability*(1/duration_acute), 
            rate_acute_chronic = (1-spontaneous_clearance_probability)*(1/duration_acute)
        )

    model = gen$new( # Arguments are (<list of model args>, <initial time step>, <number of particles>, <number of cpu threads>)
        model_params,
        0, # Initial step 
        n_particles,
        n_threads = n_threads, # Number of cpus for OMP parallelisation
        seed = 42L
        ) 

    # Print model state
    # model$state()
    # Print model info (compartment indices etc)
    # print(model$info())



    ### Test run the model to check compartment number conservation
    print("== Test running model ==")
    print(paste0("Running model with ", n_particles, " particles to test compartment number conservation"))
    t_start = Sys.time()
    # Run model for a set number of steps and return full time series
    t = 1972:2020
    t_year_int = t - t[1] + 1 # Time t_steps=1 is 1972
    t_steps = 1:(as.integer(length(t)/dt))
    # Inspect test output
    # dftest = filter$history()
    dftest = model$simulate(t_steps) # output is n state x n particles x n steps
    print(paste0("Model run took ", difftime(Sys.time(), t_start, units = "secs"), " seconds"))
    n_max_particles_plot = min(n_particles, 100, dim(dftest)[2]) # The output from a particle-filter only run is = n_particles
    index_particles_plot = sample(1:dim(dftest)[2], n_max_particles_plot) 
    dftest = dftest[, index_particles_plot, , drop=FALSE]
    # Reformat to long format
    dftest_list = list()
    for (i in 1:dim(dftest)[2]) {
        dftest_curr = as.data.frame(t(dftest[, i, ]))
        colnames(dftest_curr) = names(model$info()$index)
        dftest_curr$sim = i
        dftest_list[[i]] = dftest_curr
    }
    dftest = data.table::rbindlist(dftest_list)


    ## Sanity checks / "unit tests"
    # Check that total number of PWID is conserved
    dftest = dftest %>% mutate( # Add variables lhs for inflow and rhs for total current size + deaths
        lhs = cum_debut + 
            model_params$AN_ini + 
            # model_params$AA_ini + # These are subtracted from the total number of PWID in the model
            # model_params$AC_ini +
            model_params$AR_ini +
            model_params$EA_ini +
            model_params$EN_ini +
            model_params$EC_ini +
            model_params$ER_ini +
            model_params$QA_ini +
            model_params$QN_ini +
            model_params$QC_ini +
            model_params$QR_ini +
            # model_params$IC_ini + # Included in cum_immigrated_hcvpos
            model_params$IR_ini +
            model_params$ID_ini +
            model_params$GC_ini +
            model_params$GR_ini +
            model_params$GD_ini +
            cum_immigrated_hcvpos,
        rhs = AN + AA + AC + AR + EN + EA + EC + ER + QN + QA + QC + QR + D + IC + IR + ID + GC + GR + GD,
        rhs_lhs_diff = rhs - lhs
    )

    # # Also plot some verification/test plots
    # # Diff rhs - lhs
    # p = ggplot(data = dftest) +
    #     geom_line(
    #         aes(x = time_check, y = rhs - lhs, col=sim, group = sim),
    #         alpha = 0.3
    #     ) 
    # ggsave(plot = p, file = paste0("output_and_plots/scenario_", is, "-", scenario, "/verification/check_compartment_number_conservation-lhs_rhs_diff.png"), width = 10, height = 5)
    # # Separate plots rhs and lhs
    # p = ggplot(data = dftest) +
    #     geom_line(
    #         aes(x = time_check, y = rhs, col=sim, group = sim),
    #         alpha = 0.3
    #     ) +
    #     geom_line(
    #         aes(x = time_check, y = rhs, col=sim, group = sim),
    #         alpha = 0.3
    #     ) 
    # ggsave(plot = p, file = paste0("output_and_plots/scenario_", is, "-", scenario, "/verification/check_compartment_number_conservation-lhs_and_rhs.png"), width = 10, height = 5)
    # # rhs vs lhs
    # p = ggplot(data = dftest) +
    #     geom_line(
    #         aes(x = lhs, y = rhs, col=sim, group = sim),
    #         alpha = 0.3
    #     ) 
    # ggsave(plot = p, file = paste0("output_and_plots/scenario_", is, "-", scenario, "/verification/check_compartment_number_conservation-lhs_vs_rhs.png"), width = 10, height = 5)


    # TEST that total number of PWID is conserved
    assertthat::assert_that(
        all(abs(dftest$rhs_lhs_diff) < 1e-10),
        msg = "Total number of PWID is not conserved"
    )
    print("Model passed the checks for compartment number conservation")



    ## Instantiate particle filter
    filter = mcstate::particle_filter$new(
        data = data,
        model = gen,
        n_particles = n_particles,
        compare = n_compare,
        seed = 52L,#42L,
        n_threads = n_threads #40
    )



    ### Run pMCMC

    # Use transform function to input a mix of MCMC and constant user-defined parameters:
    transform = function(p){
        # p are the parameters under inference by the MCM
        params = list( 
            # All the user-defined input parameters.
            # This list needs to be all the required inputs minus the ones defined as MCMC
            # parameters below.
            rate_quitting_wnr_ini = model_params$rate_quitting_wnr_ini,
            rate_quitting_wr_ini = model_params$rate_quitting_wr_ini,
            rate_relapsing_ini = model_params$rate_relapsing_ini,
            sd_rw_debut = model_params$sd_rw_debut,
            AN_ini = model_params$AN_ini,
            AC_ini = model_params$AC_ini,
            EN_ini = model_params$EN_ini,
            EA_ini = model_params$EA_ini,
            EC_ini = model_params$EC_ini,
            QN_ini = model_params$QN_ini,
            QA_ini = model_params$QA_ini,
            QC_ini = model_params$QC_ini,
            IC_ini = model_params$IC_ini,
            IR_ini = model_params$IR_ini,
            ID_ini = model_params$ID_ini,
            GC_ini = model_params$GC_ini,
            GR_ini = model_params$GR_ini,
            GD_ini = model_params$GD_ini,
            dt = model_params$dt,
            rate_treatment_rw_ini = model_params$rate_treatment_rw_ini,
            rate_treatment_input_step = model_params$rate_treatment_input_step,
            time_start_rw_treatment = model_params$time_start_rw_treatment,
            sd_rtreat = model_params$sd_rtreat,
            rate_treatment_max = model_params$rate_treatment_max,
            rate_treatment_active_input_step = model_params$rate_treatment_active_input_step,
            p_treatment_succcess_step = model_params$p_treatment_succcess_step,
            reduction_factor_treatment_active = model_params$reduction_factor_treatment_active,
            rate_death_active_step = model_params$rate_death_active_step,
            rate_death_ex_step = model_params$rate_death_ex_step,
            nsp_coverage_step = model_params$nsp_coverage_step,
            risk_multiplier_lar_nsp_step = model_params$risk_multiplier_lar_nsp_step,
            gini_drug_deaths_step = model_params$gini_drug_deaths_step,
            rate_death_imm = model_params$rate_death_imm,
            rate_IC_AC = model_params$rate_IC_AC,
            immigrated_hcvpos_step = model_params$immigrated_hcvpos_step,
            rate_death_genpop = model_params$rate_death_genpop,
            genpop_hcvpos_step = model_params$genpop_hcvpos_step,
            time_forecast = model_params$time_forecast
        )
        new_params = c(params, as.list(p))

        return(new_params)
    } 

    # Define MCMC parameters with priors:
    # Param 1: beta
    beta_ini = mcstate::pmcmc_parameter("beta_ini", initial=1.36, min=0.001, max=10, prior=function(x) { dnorm(x, mean=1.36, sd=5.0, log=TRUE) })
    # Param 2: rate_debut_ini
    rate_debut_ini = mcstate::pmcmc_parameter("rate_debut_ini", initial=300, min=0.01, prior=function(x) { dnorm(x, mean=300, sd=600, log=TRUE) }) # for poisson draws
    # Param 3: AA_ini -- NB is subtracted from AN_ini inside the model
    AA_ini = mcstate::pmcmc_parameter("AA_ini", initial=200, min=0, max=model_params$AN_ini, integer=TRUE, prior=function(x) { dnorm(x, mean=200, sd=300, log=TRUE) })
    
    # Proposal matrix for MCMC steps:
    # (To get a new proposal matrix from a test run: dput(cov(pmcmc_run$pars)))
    proposal_matrix = structure(c(0.0223547541999802, 0.377581795569, -0.934171175762488, 
    0.377581795569, 6621.86306583631, -1682.43495275493, -0.934171175762488, 
    -1682.43495275493, 1729.66620958751), dim = c(3L, 3L), dimnames = list(
        c("beta_ini", "rate_debut_ini", "AA_ini"), c("beta_ini", 
        "rate_debut_ini", "AA_ini")))

    # Make mcstate pmcmc_parameters object
    mcmc_pars = mcstate::pmcmc_parameters$new(
        list(
            beta_ini=beta_ini,
            rate_debut_ini=rate_debut_ini, 
            AA_ini=AA_ini#[2:4, 2:4],
        ), 
        proposal = proposal_matrix,
        transform = transform
    )

    # Then run the model with MCMC + particle filter inference:
    print("== Running particle filtered MCMC on model ==")
    t_start = Sys.time()
    control = mcstate::pmcmc_control(
        n_steps = n_steps,
        # n_burnin = 1,
        save_state = TRUE,
        save_trajectories = TRUE,
        progress = TRUE,
        n_chains = n_chains,
        rerun_every = 10 # Rerun the particle filter for the current AND proposed every <rerun_every> MCMC steps
    )
    pmcmc_run = mcstate::pmcmc(mcmc_pars, filter, control=control)
    print(paste0("Model run took ", difftime(Sys.time(), t_start, units = "secs"), " seconds"))


    ## === MCMC diagnostics =====
    processed_chains <- mcstate::pmcmc_thin(pmcmc_run, burnin = n_burnin, thin = n_thin)
    parameter_mean_hpd <- apply(processed_chains$pars, 2, mean)
    print(parameter_mean_hpd)

    mcmc1 <- coda::as.mcmc(cbind(pmcmc_run$probabilities, pmcmc_run$pars, pmcmc_run$chain))

    # Print info on effective sample size and rejection rate
    print("coda::effectiveSize(mcmc1)")
    print(coda::effectiveSize(mcmc1))
    print("1 - coda::rejectionRate(mcmc1)")
    print(1 - coda::rejectionRate(mcmc1))


    if (n_chains > 1) { # Only applies if there are multiple chains, otherwise pmcmc_run$chain = NULL
        colnames(mcmc1)[length(colnames(mcmc1))] <- "chain"
    }
    mcmc1_thinned = coda::as.mcmc(cbind(processed_chains$probabilities, processed_chains$pars, processed_chains$chain))
    if (n_chains > 1) {
        colnames(mcmc1_thinned)[length(colnames(mcmc1_thinned))] <- "chain"
    }
    print(summary(mcmc1))
    png(filename=paste0("output_and_plots/scenario_", is, "-", scenario, "/mcmc_diagnostics.png"), width=800, height=1000)
    plot(mcmc1, ask=FALSE)#
    dev.off()

    if (n_chains > 1)
    {
        # Print gelman-rubin statistic R
        dftrace = as.data.frame(cbind(pmcmc_run$pars, pmcmc_run$chain))
        colnames(dftrace)[length(colnames(dftrace))] = "chain"
        # Calculate within-chain variance
        L = length(dftrace$chain[dftrace$chain==1])
        W = dftrace %>% group_by(chain) %>% summarise_all(var) %>% summarise_all(mean) %>% ungroup()
        # Calculate between-chain variance
        B = dftrace %>% group_by(chain) %>% summarise_all(mean) %>% ungroup() %>% 
            mutate(
                beta_ini = (beta_ini - mean(beta_ini))^2, 
                rate_debut_ini = (rate_debut_ini - mean(rate_debut_ini))^2,
                AA_ini = (AA_ini - mean(AA_ini))^2
            ) %>%
            summarise_all(sum) * L / (n_chains - 1)

        R = ((L-1)/L) + (B/W)*(1/L)
        print(R)

    }

    # Traceplot of the mcmc traces
    df_mcmc = as.data.frame(mcmc1) %>% rename("beta" = "beta_ini") %>% mutate(iteration = rep(seq(1, n_steps, by=1), n_chains), chain=factor(rep(1:n_chains, each=nrow(.)/n_chains)))
    # Save trace to rds
    saveRDS(df_mcmc, paste0("output_and_plots/scenario_", is, "-", scenario, "/mcmc_trace_full.rds"))
    df_mcmcp = df_mcmc %>% pivot_longer(cols = -c(iteration, chain), names_to = "parameter", values_to = "value")
    p = ggplot(df_mcmcp, aes(x = iteration, y = value, colour = chain, group = chain)) +
        geom_line(alpha=0.8) +
        facet_wrap(~parameter, scales = "free")
    ggsave(p, filename=paste0("output_and_plots/scenario_", is, "-", scenario, "/mcmc_traceplot.png"), width=8, height=8)


    # Correlation plot of the mcmc traces
    corrplot = TRUE
    if (corrplot) {
        p = ggpairs(data.frame(mcmc1_thinned) %>% select(!c(log_prior, log_likelihood, log_posterior)), title="correlogram") 
        ggsave(p, filename=paste0("output_and_plots/scenario_", is, "-", scenario, "/correlation_plot-thinned.png"), width=10, height=10)
    }


    ## Get output from the model
    dff = processed_chains$trajectories$state
    # Reformat df to long format for plotting
    n_particles_plot = min(dim(dff)[2], 5000) # Number of particles to plot
    index_particles_plot = sample(1:dim(dff)[2], n_particles_plot, replace=FALSE) 

    dfpl = list()
    for (i in index_particles_plot) {
        dfpl_curr = as.data.frame(t(dff[, i, ]))
        colnames(dfpl_curr) = names(model$info()$index)
        dfpl_curr$sim = i
        dfpl[[i]] = dfpl_curr
    }
    # Combine to make the main model output dataframe, in wide format (one column per output variable, sims stacked)
    df_modeloutput = data.table::rbindlist(dfpl)

    # Write to file
    saveRDS(df_modeloutput, file=paste0("output_and_plots/scenario_", is, "-", scenario, "/results.rds"))


    # Forecast using the thinned MCMC samples
    run_forecast = TRUE
    if (run_forecast) {
        ### Run prediction forward
        print("==== Running forecast on the thinned MCMC samples ====")
        time_vector_forward = seq(processed_chains$predict$time, processed_chains$predict$time+as.integer(8/dt), 1) # Forecast until 2030
        forecast = predict(
            processed_chains,
            times = time_vector_forward,
            prepend_trajectories = TRUE,
            seed = processed_chains$predict$seed
        )

        # Recalculate n_particles_plot here to be sure it is the same as for the model run
        n_particles_plot = min(dim(forecast$state)[2], 5000) # Number of particles to plot
        index_particles_plot = sample(1:dim(forecast$state)[2], n_particles_plot, replace=FALSE) 

        # Reformat df
        dfpl = list()
        for (i in index_particles_plot) {
            dfpl_curr = as.data.frame(t(forecast$state[, i, ]))
            colnames(dfpl_curr) = names(model$info()$index)
            dfpl_curr$sim = i
            dfpl[[i]] = dfpl_curr
        }
        # Combine to make the main model output dataframe, in wide format (one column per output variable, sims stacked)
        df_forecast = data.table::rbindlist(dfpl)
        
        # Thin to integer time points
        df_forecast = df_forecast %>% filter((df_forecast$t %% 1) < 0.0001)

        # Write to file
        saveRDS(df_forecast, file=paste0("output_and_plots/scenario_", is, "-", scenario, "/results_with_forecast.rds"))

    }

    # Source to run plotting scripts
    source("plots.R")

} # End loop over scenarios