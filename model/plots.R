library(tidyr)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(reshape2)
library(fhiplot) # This is used to get the FHI color palette and stuff
# fhiplot can be disabled manually in the code below or the package installed as follows.
# In an R session paste:
# options(repos=structure(c(
#  FHI="https://folkehelseinstituttet.github.io/drat/",
#  CRAN="https://cran.rstudio.com"
#)))
# install.packages("fhiplot")



# ======= PLOTTING =======
print("===== Running plots.R =====")
# Particle filter output is in the file results.rds / results_with_forecast.rds in the folder
# specified by root_folder. 
# We put all the plots for report in the folder plots under root_folder
root_folder = file.path(getwd(), "output_and_plots", paste0("scenario_", is, "-", scenario)) # Assuming called from run_hcv_model.R under git model directory
print(paste0("root_folder: ", root_folder))

label_list = c(
            "total_chronic_hcv_prevalence" = "Totalt",
            "total_chronic_hcv_prevalence_pwid"="Personer som har inntatt rusmidler med sprøyte",
            "chronic_hcv_prevalence_active_pwid"="Personer som aktivt inntar rusmidler med sprøyte",
            "chronic_hcv_prevalence_ex_pwid"="Personer som tidligere har innt rusmidler med sprøyte",
            "IC"="Innvandrere",
            "inc_infections" = "Personer som har inntatt rusmidler med sprøyte (akutt infeksjon)",
            "inc_immigrated_hcvpos" = "Innvandrere (kronisk infeksjon)",
            "AN_AA_AC_AR" = "Aktivt injiserende", 
            "EN_EA_EC_ER" = "Midlertidig holdt opp", 
            "QN_QA_QC_QR" = "Permanent holdt opp"
        )



for (with_forecast in c(FALSE, TRUE)) {
    forecast_string = ifelse(with_forecast, "-with_forecast", "") # For file saving

    filename = paste0(root_folder, "/results", ifelse(with_forecast, "_with_forecast", ""), ".rds")
    # Make folder plots under root_folder if it does not exist
    if (!dir.exists(paste0(root_folder, "/plots"))) {
        dir.create(paste0(root_folder, "/plots"))
    }

    df = readRDS(filename)

    if (!with_forecast) {
        # Add immigration population size by year, only for non-forecast
        df_immigrants = read.csv("/home/jorgenem/gitrepos/HCV/data/immigrants/ssb_immigrant_prevalence-for_model.csv") %>%
            group_by(Year) %>%
            summarise(
                pop_immigrants = sum(prevalence)
            )
        # Interpolate prevalence linearly to the years in model, 1972:2023
        immigrants_interpolated = approx(
            x=df_immigrants$Year,
            y=df_immigrants$pop_immigrants,
            xout=1972:2023,
            method="linear"
        )
        df_immigrants = data.frame(time_check=immigrants_interpolated$x-1971, pop_immigrants=immigrants_interpolated$y)
        # Join immigrants prevalence to df
        df = df %>% left_join(df_immigrants, by=c("time_check"))

        # Add general population size by year, only for non-forecast
        df_genpop = fhidata::norway_population_b2020 %>% 
            filter(granularity_geo == "nation") %>%
            group_by(year) %>%
            summarise(pop_national = sum(pop)) %>% 
            ungroup() %>%
            mutate(time_check = year-1971)
        # Join general population prevalence to df
        df = df %>% left_join(df_genpop, by=c("time_check"))
    }

    # Add some convenience columns that are summed
    df = df %>% mutate(
        AA_AC = AA + AC,
        EA_EC = EA + EC,
        QA_QC = QA + QC,
        AN_AA_AC_AR = AN + AA + AC + AR,
        EN_EA_EC_ER = EN + EA + EC + ER,
        QN_QA_QC_QR = QN + QA + QC + QR,
        frac_pos_active = (AA + AC) / (AN + AA + AC + AR),
        frac_recovered_active = AR / (AN + AA + AC + AR),
        frac_recovered_ex = (ER + QR) / (EN + EA + EC + ER + QN + QA + QC + QR),
        frac_antibodies_active = (AA + AC + AR) / (AN + AA + AC + AR),
        frac_antibodies_ex = (EA + EC + ER + QA + QC + QR) / (EN + EA + EC + ER + QN + QA + QC + QR),
        inc_treatments = inc_treatments_pwid + inc_treatments_immigrants + inc_treatments_genpop,
        total_chronic_hcv_prevalence = AC + EC + QC + IC + GC,
        # total_hcv_prevalence_pwid = AA + AC + EA + EC + QA + QC,
        total_acute_hcv_prevalence_pwid = AA + EA + QA,
        total_chronic_hcv_prevalence_pwid = AC + EC + QC,
        chronic_hcv_prevalence_active_pwid = AC,
        chronic_hcv_prevalence_ex_pwid = EC + QC,
    ) %>%
    rename(t = time_check)

    # Add cumulative treatments to df by summing inc_treatments for each sim
    df = df %>% group_by(sim) %>% mutate(cumulative_treatments = cumsum(inc_treatments)) %>% ungroup()

    if (!with_forecast) {
        df = df %>% mutate(
            frac_pos_immigrants = IC / pop_immigrants,
            frac_chronic_tot_pop = total_chronic_hcv_prevalence / pop_national
        )
    }




    year1_list = c(1972)
    for (year1 in year1_list) { # Cut time at some initial cutoff point
        # Pivot to long format (one row per output variable, sims stacked)
        dfp = df %>% pivot_longer(!c(sim, t), values_to="count", names_to="type")
        dfp = dfp %>% filter(t >= year1 - 1971) # Cut time at some initial cutoff point
        dfp = dfp %>% mutate(grouping = paste0(type, "_", sim))
        dfp = dfp %>% mutate(type = factor(type, levels = c( # Fix ordering for labels below
            "total_chronic_hcv_prevalence", 
            "total_chronic_hcv_prevalence_pwid",
            unique(dfp$type[!dfp$type %in% c("total_chronic_hcv_prevalence", "total_chronic_hcv_prevalence_pwid")])  
        )))
        dfp = dfp %>% mutate(year = t + 1971)

        # Calculate median + confidence intervals (before subsampling!)
        dfpc = dfp %>% 
            group_by(type, t, year) %>% 
            summarise_at(vars(count),
            list(
                # min=min, 
                median=median, 
                lower_95perc=~quantile(., probs = 0.025),
                upper_95perc=~quantile(., probs = 0.975),
                lower_IQR=~quantile(., probs = 0.25),
                upper_IQR=~quantile(., probs = 0.75)
                # max=max
                )) %>%
            ungroup()

        # Subsample dfp to get a more manageable number of lines for plotting (note that error bands are made before subsampling)
        sims_to_plot = sample(unique(dfp$sim), min(length(unique(dfp$sim)), 100))
        dfp = dfp %>% filter(sim %in% sims_to_plot)
        # Save to file, but only for full timeperiod (should be first elemnet in year1_list)
        if (year1 == year1_list[1]) {
            output_types = c(
                "total_chronic_hcv_prevalence", "total_chronic_hcv_prevalence_pwid", "chronic_hcv_prevalence_active_pwid", "chronic_hcv_prevalence_ex_pwid", "total_acute_hcv_prevalence_pwid", "IC", "inc_infections", 
                "inc_immigrated_hcvpos", "inc_treatments", "frac_pos_active", "frac_recovered_active", "frac_antibodies_active"
            )
            write.csv(dfpc %>% filter(type %in% output_types), paste0(root_folder, "/plots/model_results", forecast_string, ".csv"))
        }


        # == Plot chronic HCV prevalence, total and by group
        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("total_chronic_hcv_prevalence", "total_chronic_hcv_prevalence_pwid", "IC")), 
                aes(x = t+1971, y = count, col = type, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/total_HCV_prevalence_and_split_by_group-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("total_chronic_hcv_prevalence", "total_chronic_hcv_prevalence_pwid", "IC"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col=type), alpha=1) +
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/total_HCV_prevalence_and_split_by_group-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)


        # == Plot chronic HCV prevalence, total only
        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("total_chronic_hcv_prevalence")), 
                aes(x = t+1971, y = count, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "none"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/total_HCV_prevalence-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("total_chronic_hcv_prevalence"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col=type), alpha=1) +
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
                # colour = ""
            ) +
            # guides(fill_ramp = "none") + # no legend for level
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "none"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/total_HCV_prevalence-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        if (!with_forecast) {
            # == Plot chronic HCV prevalence, total only, as share of total population
            # Line plot
            p = ggplot() + 
                geom_line(
                    data = dfp %>% filter(type %in% c("frac_chronic_tot_pop")), 
                    aes(x = t+1971, y = count*100, group = grouping), 
                    alpha = 0.3
                ) + 
                labs(
                    x = "År", 
                    y = "Estimert prevalens av kronisk hepatitt C i befolkninga (%)"
                    # colour = ""
                ) +
                fhiplot::scale_color_fhi(
                    name = "",
                    palette = "primary",
                    labels = label_list
                ) +
                fhiplot::theme_fhi_lines() +
                theme_bw(
                    base_size = 12
                ) +
                theme(
                    legend.position = "none"
                ) +
                guides(colour = guide_legend(override.aes = list(alpha = 1)))
            ggsave(plot = p, file = paste0(root_folder, "/plots/total_HCV_prevalence_as_share_of_national_population-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

            # Error band plot
            p = ggplot(data=dfpc %>% filter(type %in% c("frac_chronic_tot_pop"))) + 
                geom_ribbon(aes(x=t+1971, ymin=lower_95perc*100, ymax=upper_95perc*100, fill=type), alpha=0.3) +
                geom_ribbon(aes(x=t+1971, ymin=lower_IQR*100, ymax=upper_IQR*100, fill=type), alpha=0.3) +
                geom_line(aes(x=t+1971, y=median*100, col=type), alpha=1) +
                labs(
                    x = "År", 
                    y = "Estimert prevalens av kronisk hepatitt C i befolkninga (%)"
                    # colour = ""
                ) +
                # guides(fill_ramp = "none") + # no legend for level
                fhiplot::scale_color_fhi(
                    name = "",
                    palette = "primary",
                    labels = label_list
                ) +
                fhiplot::scale_fill_fhi(
                    name = "",
                    palette = "primary",
                    labels = label_list
                ) +
                fhiplot::theme_fhi_lines() +
                theme_bw(
                    base_size = 12
                ) +
                theme(
                    legend.position = "none"
                ) +
                guides(colour = guide_legend(override.aes = list(alpha = 1)))
            ggsave(plot = p, file = paste0(root_folder, "/plots/total_HCV_prevalence_as_share_of_national_population-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)
        }

        # == Plot chronic HCV prevalence, individual groups only
        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("total_chronic_hcv_prevalence_pwid", "IC")), 
                aes(x = t+1971, y = count, col = type, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_prevalence_split_by_group-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("total_chronic_hcv_prevalence_pwid", "IC"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col=type), alpha=1) +
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_prevalence_split_by_group-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # == Plot chronic HCV prevalence, PWID only
        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("total_chronic_hcv_prevalence_pwid")), 
                aes(x = t+1971, y = count, col = type, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "none"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_prevalence_PWID_only-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("total_chronic_hcv_prevalence_pwid"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col=type), alpha=1) +
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "none"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_prevalence_PWID_only-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # == Plot chronic HCV prevalence, PWID split by active and ex
        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("chronic_hcv_prevalence_active_pwid", "chronic_hcv_prevalence_ex_pwid")), 
                aes(x = t+1971, y = count, col = type, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_prevalence_PWID_only-split_by_ex_active-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("chronic_hcv_prevalence_active_pwid", "chronic_hcv_prevalence_ex_pwid"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col=type), alpha=1) +
            labs(
                x = "År", 
                y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_prevalence_PWID_only-split_by_ex_active-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # == Plot acute HCV prevalence, PWID only
        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("total_acute_hcv_prevalence_pwid")), 
                aes(x = t+1971, y = count, col = type, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Estimert prevalens av akutt hepatitt C (antall personer)"
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "none"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/acute_HCV_prevalence_PWID_only-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("total_acute_hcv_prevalence_pwid"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col=type), alpha=1) +
            labs(
                x = "År", 
                y = "Estimert prevalens av akutt hepatitt C (antall personer)"
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "none"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/acute_HCV_prevalence_PWID_only-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        if (!with_forecast) {
            # == Plot HCV prevalence, immigrants only
            # Line plot
            p = ggplot() + 
                geom_line(
                    data = dfp %>% filter(type %in% c("IC")), 
                    aes(x = t+1971, y = count, col = type, group = grouping), 
                    alpha = 0.3
                ) + 
                labs(
                    x = "År", 
                    y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
                    # colour = ""
                ) +
                fhiplot::scale_color_fhi(
                    name = "",
                    palette = "primary",
                    labels = label_list
                ) +
                fhiplot::theme_fhi_lines() +
                theme_bw(
                    base_size = 12
                ) +
                theme(
                    legend.position = "none"
                ) +
                guides(colour = guide_legend(override.aes = list(alpha = 1)))
            ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_prevalence_immigrants_only-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

            # Error band plot
            p = ggplot(data=dfpc %>% filter(type %in% c("IC"))) + 
                geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill=type), alpha=0.3) +
                geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill=type), alpha=0.3) +
                geom_line(aes(x=t+1971, y=median, col=type), alpha=1) +
                labs(
                    x = "År", 
                    y = "Estimert prevalens av kronisk hepatitt C (antall personer)"
                    # colour = ""
                ) +
                # guides(fill_ramp = "none") + # no legend for level
                fhiplot::scale_color_fhi(
                    name = "",
                    palette = "primary",
                    labels = label_list
                ) +
                fhiplot::scale_fill_fhi(
                    name = "",
                    palette = "primary",
                    labels = label_list
                ) +
                fhiplot::theme_fhi_lines() +
                theme_bw(
                    base_size = 12
                ) +
                theme(
                    legend.position = "none"
                ) +
                guides(colour = guide_legend(override.aes = list(alpha = 1)))
            ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_prevalence_immigrants_only-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)


            # == Plot fractional HCV prevalence, immigrants only
            # Line plot
            p = ggplot() + 
                geom_line(
                    data = dfp %>% filter(type %in% c("frac_pos_immigrants")), 
                    aes(x = t+1971, y = count*100, col = type, group = grouping), 
                    alpha = 0.3
                ) + 
                labs(
                    x = "År", 
                    y = "Estimert andel av innvandrere som har kronisk hepatitt C (%)"
                    # colour = ""
                ) +
                fhiplot::scale_color_fhi(
                    name = "",
                    palette = "primary",
                    labels = label_list
                ) +
                fhiplot::theme_fhi_lines() +
                theme_bw(
                    base_size = 12
                ) +
                theme(
                    legend.position = "none"
                ) +
                guides(colour = guide_legend(override.aes = list(alpha = 1)))
            ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_percent_prevalence_immigrants-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

            # Error band plot
            p = ggplot(data=dfpc %>% filter(type %in% c("frac_pos_immigrants"))) + 
                geom_ribbon(aes(x=t+1971, ymin=lower_95perc*100, ymax=upper_95perc*100, fill=type), alpha=0.3) +
                geom_ribbon(aes(x=t+1971, ymin=lower_IQR*100, ymax=upper_IQR*100, fill=type), alpha=0.3) +
                geom_line(aes(x=t+1971, y=median*100, col=type), alpha=1) +
                labs(
                    x = "År", 
                    y = "Estimert andel av innvandrere som har kronisk hepatitt C (%)"
                ) +
                fhiplot::scale_color_fhi(
                    name = "",
                    palette = "primary",
                    labels = label_list
                ) +
                fhiplot::scale_fill_fhi(
                    name = "",
                    palette = "primary",
                    labels = label_list
                ) +
                fhiplot::theme_fhi_lines() +
                theme_bw(
                    base_size = 12
                ) +
                theme(
                    legend.position = "none"
                ) +
                guides(colour = guide_legend(override.aes = list(alpha = 1)))
            ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_percent_prevalence_immigrants-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)
        } # End if not forecast

        # END prevalence



        # === Fractional prevalence ===
        # == Plot HCV RNA prevalence as fraction, PWID only
        # Load df_hcv
        df_hcv = read.csv("/home/jorgenem/gitrepos/HCV/data/pwid_HCV_prevalence-all_surveys_combined-with_uncertainty_columns_for_plotting.csv")

        # Error bar plot
        p = ggplot(data = dfpc %>% filter(type %in% c("frac_pos_active"))
            ) +
            geom_ribbon(aes(x=t+1971, ymin=100*lower_95perc, ymax=100*upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=100*lower_IQR, ymax=100*upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=100*median, linetype=type), alpha=1) +
            geom_point(
                data = df_hcv,
                aes(x=year_int+1971, y=100*mean, col=city)
            ) + 
            geom_errorbar(
                data = df_hcv,
                aes(x=year_int+1971, ymin=100*lower, ymax=100*upper, col=city)
            ) +
            labs(
                x = "År", 
                y = "Estimert prevalens av HCV-RNA-positive (% av aktivt injiserende personer)"
                # colour = ""
            ) +
            scale_fill_manual(
                name = "",
                values = c("grey"),
                labels = list("frac_pos_active" = "Modellestimat")
            ) +
            scale_linetype_manual(
                name = "",
                values = c("solid"),
                labels = list("frac_pos_active" = "Modellestimat")
            ) +
            scale_colour_manual(
                name = "Data",
                values = c("violet", "orange", "red", "blue"),
                labels = list("frac_pos_active" = "Modellestimat", "Oslo" = "Oslo", "Bergen_and_Stavanger" = "Bergen og Stavanger", "Trondheim" = "Trondheim", "Bergen" = "Bergen")
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(linetype = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_RNA_percent_prevalence-PWID-errorbar-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("frac_pos_active")), 
                aes(x = t+1971, y = 100*count, linetype = type, group = grouping), 
                alpha = 0.1
            ) + 
            geom_point(
                data = df_hcv,
                aes(x=year_int+1971, y=100*mean, col=city)
            ) + 
            geom_errorbar(
                data = df_hcv,
                aes(x=year_int+1971, ymin=100*lower, ymax=100*upper, col=city)
            ) +
            labs(
                x = "År", 
                y = "Estimert prevalens av HCV-RNA-positive (% av aktivt injiserende)"
                # colour = ""
            ) +
            scale_linetype_manual(
                name = "",
                values = c("solid"),
                labels = list("frac_pos_active" = "Modellestimat")
            ) +
            scale_colour_manual(
                name = "Data",
                values = c("violet", "orange", "red", "blue"),
                labels = list("frac_pos_active" = "Modellestimat", "Oslo" = "Oslo", "Bergen_and_Stavanger" = "Bergen og Stavanger", "Trondheim" = "Trondheim", "Bergen" = "Bergen")
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(linetype = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_RNA_percent_prevalence-PWID-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)


        # == Plot HCV antibody prevalence as fraction, PWID only
        # Error bar plot
        p = ggplot(data = dfpc %>% filter(type %in% c("frac_antibodies_active", "frac_antibodies_ex"))
            ) +
            geom_ribbon(aes(x=t+1971, ymin=100*lower_95perc, ymax=100*upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=100*lower_IQR, ymax=100*upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=100*median, linetype=type), alpha=1) +
            geom_point(
                data = df_hcv %>% filter(city %in% c("Oslo")),
                aes(x=year_int+1971, y=100*mean_ag, col=city)
            ) + 
            geom_errorbar(
                data = df_hcv %>% filter(city %in% c("Oslo")),
                aes(x=year_int+1971, ymin=100*lower_ag, ymax=100*upper_ag, col=city)
            ) +
            labs(
                x = "År", 
                y = "Estimert prevalens av HCV-antistoff-positive (% av aktivt injiserende personer)"
            ) +
            fhiplot::scale_fill_fhi(
                name = "Modellestimat",
                palette = "primary",
                labels = list("frac_antibodies_active"="Aktive", "frac_antibodies_ex"="Eks-aktive")
            ) +
            scale_linetype_manual(
                name = "Modellestimat",
                values = c("solid", "solid"),
                labels = list("frac_antibodies_active"="Aktive", "frac_antibodies_ex"="Eks-aktive")
            ) +
            scale_color_brewer(
                name = "Data",
                palette = "BrBG",
                labels = list("Oslo" = "Oslo", "Bergen_and_Stavanger" = "Bergen og Stavanger", "Trondheim" = "Trondheim")
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(linetype = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_antigen_percent_prevalence-PWID-errorbar-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("frac_antibodies_active", "frac_antibodies_ex")), 
                aes(x = t+1971, y = 100*count, color=type, group = grouping), 
                alpha = 0.3
            ) + 
            geom_point(
                data = df_hcv %>% filter(city %in% c("Oslo")),
                aes(x=year_int+1971, y=100*mean_ag, col=city)
            ) + 
            geom_errorbar(
                data = df_hcv %>% filter(city %in% c("Oslo")),
                aes(x=year_int+1971, ymin=100*lower_ag, ymax=100*upper_ag, col=city)
            ) +
            labs(
                x = "År", 
                y = "Estimert prevalens av HCV-antistoff-positive (% av aktivt injiserende personer)"
            ) +
            scale_linetype_manual(
                name = "",
                values = c("solid", "solid"),
                labels = list("frac_antibodies_active"="Modell, aktive", "frac_antibodies_ex"="Modell, eks-aktive")
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = list("frac_antibodies_active"="Modell, aktive", "frac_antibodies_ex"="Modell, eks-aktive", "Oslo" = "Data, Oslo", "Bergen_and_Stavanger" = "Bergen og Stavanger", "Trondheim" = "Trondheim")
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(
                linetype="none"
            )
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_antigen_percent_prevalence-PWID-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)





        # == Plot HCV incidence for PWID only
        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("inc_infections")), 
                aes(x = t+1971, y = count, col = type, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Estimert årlig insidens av nye infeksjoner med hepatitt C"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "none"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_incidence_PWID-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("inc_infections"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col=type), alpha=1) +
            labs(
                x = "År", 
                y = "Estimert årlig insidens av nye infeksjoner med hepatitt C"
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "none"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1)))
        ggsave(plot = p, file = paste0(root_folder, "/plots/HCV_incidence_PWID-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)


        
        # === Population size, PWID ===
        # Load and reformat PWID dataset for plotting
        df_pwid = read.csv("/home/jorgenem/gitrepos/HCV/data/number_of_PWID-2000_peak_modified.csv")
        df_pwid_plot = rbind(
            data.frame(
                year=df_pwid$year,
                type="active",
                count=df_pwid$N_active,
                upper=df_pwid$N_active_upper,
                lower=df_pwid$N_active_lower
            ),
            data.frame(
                year=df_pwid$year,
                type="ex_will_relapse",
                count=df_pwid$N_ex_wr,
                upper=NA,
                lower=NA
            ),
            data.frame(
                year=df_pwid$year,
                type="ex_will_not_relapse",
                count=df_pwid$N_ex_wnr,
                upper=NA,
                lower=NA
            )
        ) %>% mutate(
            type = factor(type, levels = c("active", "ex_will_relapse", "ex_will_not_relapse"))
        )

        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("AN_AA_AC_AR", "EN_EA_EC_ER", "QN_QA_QC_QR")), 
                aes(x = t+1971, y = count, col = type, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Antall personer"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "Modell",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(
                colour = guide_legend(
                    override.aes = list(alpha = 1),
                    nrow = 3
                )
            ) +
            # New colour scale for data points using ggnewscale
            ggnewscale::new_scale_colour() + 
            geom_point(
                data = df_pwid_plot,
                aes(x=year, y=count, col=type)
            ) +
            geom_errorbar(
                data = df_pwid_plot,
                aes(x=year, y=count, ymin=lower, ymax=upper, col=type)
            ) +
            fhiplot::scale_color_fhi(
                name = "Data",
                palette = "warning",
                labels = c(
                    "active" = "Aktivt injiserende", 
                    "ex_will_relapse" = "Midlertidig holdt opp", 
                    "ex_will_not_relapse" = "Permanent holdt opp"
                )
            ) +
            guides(
                colour = guide_legend(
                    nrow = 3
                )
            )
        ggsave(plot = p, file = paste0(root_folder, "/plots/PWID_population-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("AN_AA_AC_AR", "EN_EA_EC_ER", "QN_QA_QC_QR"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill=type), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill=type), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col=type), alpha=1) +
            labs(
                x = "År", 
                y = "Antall personer"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "none"
            ) +
            guides(colour = guide_legend(override.aes = list(alpha = 1))) +
            # New colour scale for data points using ggnewscale
            ggnewscale::new_scale_colour() + 
            geom_point(
                data = df_pwid_plot,
                aes(x=year, y=count, col=type)
            ) +
            geom_errorbar(
                data = df_pwid_plot,
                aes(x=year, y=count, ymin=lower, ymax=upper, col=type)
            ) +
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = c(
                    "active" = "Aktivt injiserende", 
                    "ex_will_relapse" = "Midlertidig holdt opp", 
                    "ex_will_not_relapse" = "Permanent holdt opp"
                )
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(
                colour = guide_legend(
                    nrow = 1
                )
            )
        ggsave(plot = p, file = paste0(root_folder, "/plots/PWID_population-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)


        # === Plot incidence treatments and data ===
        # Load treatment data
        df_treatment = read.csv("/home/jorgenem/gitrepos/HCV/data/treatments-93percent.csv") 
        # Add column for cumulative treatments
        df_treatment = df_treatment %>% mutate(
            cumulative_n_treated_total = cumsum(n_treated_total)
        )

        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("inc_treatments")), 
                aes(x = t+1971, y = count, colour=type, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Antall behandlinger per år"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "Modell",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "none"
            ) +
            geom_point(
                data = df_treatment,
                aes(x=year, y=n_treated_total, colour="data")
            ) 
        ggsave(plot = p, file = paste0(root_folder, "/plots/treatments-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("inc_treatments"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill="Modell"), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill="Modell"), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col="Modell"), alpha=1) +
            labs(
                x = "År", 
                y = "Antall behandlinger per år"
                # colour = ""
            ) +
            # guides(fill_ramp = "none") + # no legend for level
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(fill="none") +
            geom_point(
                data = df_treatment,
                aes(x=year, y=n_treated_total, col="Data")
            ) 
        ggsave(plot = p, file = paste0(root_folder, "/plots/treatments-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Also plot cumulative treatments
        # Line plot
        p = ggplot() + 
            geom_line(
                data = dfp %>% filter(type %in% c("cumulative_treatments")), 
                aes(x = t+1971, y = count, colour=type, group = grouping), 
                alpha = 0.3
            ) + 
            labs(
                x = "År", 
                y = "Kumulativt antall behandlinger"
                # colour = ""
            ) +
            fhiplot::scale_color_fhi(
                name = "Modell",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 10
            ) +
            theme(
                legend.position = "none"
            ) +
            geom_point(
                data = df_treatment,
                aes(x=year, y=cumulative_n_treated_total, colour="data")
            )
        ggsave(plot = p, file = paste0(root_folder, "/plots/cumulative_treatments-trajectories-cutoff_", year1, forecast_string, ".png"), width=8, height=6)

        # Error band plot
        p = ggplot(data=dfpc %>% filter(type %in% c("cumulative_treatments"))) + 
            geom_ribbon(aes(x=t+1971, ymin=lower_95perc, ymax=upper_95perc, fill="Modell"), alpha=0.3) +
            geom_ribbon(aes(x=t+1971, ymin=lower_IQR, ymax=upper_IQR, fill="Modell"), alpha=0.3) +
            geom_line(aes(x=t+1971, y=median, col="Modell"), alpha=1) +
            labs(
                x = "År", 
                y = "Kumulativt antall behandlinger"
                # colour = ""
            ) +
            # guides(fill_ramp = "none") + # no legend for level
            fhiplot::scale_color_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::scale_fill_fhi(
                name = "",
                palette = "primary",
                labels = label_list
            ) +
            fhiplot::theme_fhi_lines() +
            theme_bw(
                base_size = 12
            ) +
            theme(
                legend.position = "bottom"
            ) +
            guides(fill="none") +
            geom_point(
                data = df_treatment,
                aes(x=year, y=cumulative_n_treated_total, col="Data")
            )
        ggsave(plot = p, file = paste0(root_folder, "/plots/cumulative_treatments-errorband-cutoff_", year1, forecast_string, ".png"), width=8, height=6)



    } # End loop over cutoff year
} # End loop over with/without forecast