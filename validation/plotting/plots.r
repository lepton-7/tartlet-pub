{
    library(tidyverse)
    library(ggplot2)
    setwd("validation/plotting")
}

{
    aliased_results <- read.csv("data/results_alias.csv")
    rfam_regions <- read.csv("data/all_regions_14.10.csv")
}

# Colour palattes
{
    decision_pal <- c("pass" = "#6956e7", "fail" = "#eaef5b", "incon" = "#000000")
}

# -----------------------------------------------------------------------------
# Stacked decision matrix faceted by condition and dataset/strain riboswitch

{
    dset <- "s_sanguinis"
    df <- aliased_results[aliased_results[, "dataset"] == dset, ]

    df |>
        count(rowid_alias, condition_alias, decision) |>
        mutate(prop = n / sum(n), .by = c(rowid_alias, condition_alias)) |>
        ggplot(aes(x = 0, y = prop, fill = decision)) +
        geom_col(position = "stack") +
        scale_fill_manual(name = "decision", values = decision_pal) +
        scale_x_continuous(breaks = NULL, labels = NULL) +
        facet_grid(reorder(rowid_alias, -n) ~ reorder(condition_alias, -n)) +
        labs(x = NULL, y = "Proportion", fill = "Outcome") +
        theme_minimal() +
        theme(
            strip.text.y = element_text(size = 12, angle = 0),
            strip.text.x = element_text(size = 12, angle = 90)
        )
}

{
    save_path <- sprintf("plots/decision_matrix/%s.png", dset)
    ggsave(save_path, dpi = 320, units = "px")
}

# -----------------------------------------------------------------------------
# Box plot for riboswitch sizes across rfam families for all rfam riboswitches
{
    rfam_region_sizes <- ggplot(rfam_regions, aes(x = rfam_id, y = region_size)) +
        geom_point(
            position = position_jitterdodge(jitter.width = 0.5),
            aes(alpha = 0.01), color = "#eea8a8"
        ) +
        geom_boxplot(
            outlier.shape = NA
        ) +
        guides(alpha = "none") +
        theme(
            # plot.margin = unit(c(1, 1, 1, 1), "lines"),
            axis.title.x = element_blank(),
            axis.text.x = element_text(vjust = 1, hjust = 1, size = 13, angle = 60, colour = "black"),
            # axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 13, face = "bold"),
            axis.title = element_text(size = 18, face = "bold"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.title = element_blank()
        ) +
        labs(
            y = "Riboswitch size (nt)"
        )

    rfam_region_sizes
}

{
    save_path <- "plots/rfam_region_sizes.svg"
    ggsave(save_path, dpi = 320, units = "px")
}
