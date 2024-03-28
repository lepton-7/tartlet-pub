{
    library(tidyverse)
    library(ggplot2)
    library(rjson)
    library(reshape2)
    library(ggtree)
    library(stringr)
    library(patchwork)
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
    save_path <- "plots/rfam_region_sizes.pdf"
    ggsave(save_path, dpi = 320, units = "px")
}

# -----------------------------------------------------------------------------
# Density plot of sizes for cumulative rfam riboswitches
{
    rfam_size_density <- ggplot(rfam_regions, aes(x = region_size)) +
        # geom_histogram(aes(y = ..density..), colour = "black", fill = "white") +
        geom_density(alpha = .2, fill = "#FF6666") +
        theme_classic() +
        theme(
            axis.text.x = element_text(size = 13, angle = 0, colour = "black"),
            axis.text.y = element_text(size = 13, face = "bold"),
            axis.title = element_text(size = 18, face = "bold"),
            plot.title = element_blank()
        ) +
        labs(
            y = "Density",
            x = "Riboswitch size"
        )

    rfam_size_density
}

{
    save_path <- "plots/rfam_region_size_density.svg"
    ggsave(save_path, dpi = 320, units = "px", width = 2000, height = 1500)
}


# -----------------------------------------------------------------------------
# Density plots of rfam riboswitch sizes faceted by riboswitch IDs

{
    rfam_size_density_facet_id <- rfam_size_density +
        facet_wrap(~rfam_id, scales = "free") +
        theme(
            axis.text.x = element_text(size = 9, angle = 0, colour = "black"),
            axis.text.y = element_text(size = 9, face = "bold"),
        ) +
        theme_classic()

    rfam_size_density_facet_id
}

{
    save_path <- "plots/rfam_region_size_density_by_id.svg"
    ggsave(save_path, dpi = 320, units = "px", width = 4500, height = 3000)
}

# -----------------------------------------------------------------------------
# Density plot of read sizes
{
    combo <- dplyr::bind_rows(list(dat1 = rfam_regions, dat2 = read_sizes),
        .id = "dataset"
    )

    bw <- 10

    fragment_size_comparison <- ggplot(combo) +
        geom_histogram(aes(x = region_size, y = ..density..), colour = "black", fill = "#FF6666", binwidth = bw) +
        geom_histogram(aes(x = sizes, y = ..density.., alpha = 0.5), colour = "#626262", fill = "white", binwidth = bw) +
        # geom_density(alpha = .2, fill = "#FF6666") +
        def_theme +
        guides(alpha = "none") +
        labs(
            y = "Density",
            x = "Size (nt)"
        ) +
        xlim(0, 1000)

    fragment_size_comparison
}

{
    save_path <- "plots/fragment_switch_size_comparison.svg"
    ggsave(save_path, dpi = 320, units = "px", width = 3000, height = 1500)
}

# -----------------------------------------------------------------------
# Compare fragment size distributions across SRAs

{
    dset <- "b_sub_168"
    read_sizes <- read.csv("data/b_sub_168_fragment_sizes.csv")

    melt_sizes <- melt(read_sizes, id = c("size"), value.name = "freq", variable.name = "sra", na.rm = TRUE)
    melt_sizes$freq <- as.numeric(melt_sizes$freq)
}

{
    frag_size_distr <- ggplot(melt_sizes, aes(x = size, weight = freq)) +
        geom_histogram(
            aes(
                y = after_stat(density),
                alpha = 0.5
            ),
            colour = "#464646", fill = "#f96d6d", binwidth = 20
        ) +
        facet_wrap(~sra) +
        guides(alpha = "none") +
        def_theme +
        labs(
            x = "Fragment size (nt)",
            y = "Density"
        ) +
        scale_x_continuous(breaks = seq(250, 1000, 500), limits = c(0, 1000)) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 2))

    frag_size_distr
}

{
    alph <- 2
    save_path <- str_glue("plots/{dset}_fragment_sizes.png")
    ggsave(save_path, dpi = 320 * alph, units = "px", width = 5000 * alph, height = 3500 * alph)
}

# -----------------------------------------------------------------------
# The big results fig with a tax tree and inferred riboswitch mechs

{
    tree <- read.tree("data/big_fig/taxa_list.tree")
    phylo_dict <- fromJSON(file = "data/big_fig/fancy_tips.json")

    tax_tree <- ggtree(tree, size = 0.8, layout = "rectangular", ladderize = FALSE, branch.length = -10) +
        ggtitle("Validation set tax tree") +
        geom_tiplab(size = 8, linesize = .5, align = TRUE) +
        xlim_tree(5)

    for (i in seq_along(tree$tip.label)) {
        tax_tree$data$label[i] <- phylo_dict[[tax_tree$data$label[i]]][2]
    }

    tax_tree
}
