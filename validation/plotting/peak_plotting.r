{
    library(tidyverse)
    library(ggplot2)
    library(reshape2)
    library(stringr)
    setwd("validation/plotting")
}

# themes
{
    def_theme <- theme_classic() +
        theme(
            axis.text.x = element_text(size = 13, angle = 0, colour = "black"),
            axis.text.y = element_text(size = 13, face = "bold"),
            axis.title = element_text(size = 18, face = "bold"),
            plot.title = element_blank()
        )
}

{
    dset <- "b_sub_168"
    # dset <- "e_coli"
    df <- read.csv(str_glue("../alignment/outputs/{dset}/plots/peak_log.csv"))

    df$pval <- -log(df$coverage_delta_pval)
    df$pval[df$pval > 10] <- 10
}

{
    peakplot <- ggplot(df, aes(x = from_riboswith_end_relative, y = coverage_delta)) +
        # geom_smooth(⁠method = "gam", formula = y ~ s(x, bs = "cs"), se = TRUE) +
        geom_line(aes(color = transcriptome)) +
        geom_point(aes(alpha = pval)) +
        # stat_summary(aes(group = rowid), fun = mean, geom = "line", colour = "#000000", size = 2) +
        # stat_summary(aes(group = rowid), geom = "ribbon", fun.data = mean_cl_normal, width = 0.1, conf.int = 0.95, fill = "lightblue") +
        # scale_alpha(range = c(0, 1)) +
        # facet_wrap(~rowid, scales = "free") +
        facet_wrap(~rowid, scales = "fixed") +
        labs(
            title = str_glue("{dset} riboswitches across all available conditions")
        ) +
        def_theme +
        theme(
            plot.title = element_text(size = 14)
        )

    peakplot
}
