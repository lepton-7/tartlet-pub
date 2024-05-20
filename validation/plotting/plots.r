{
    library(tidyverse)
    library(ggplot2)
    library(rjson)
    library(reshape2)
    library(ggtree)
    library(stringr)
    library(patchwork)
    library(dplyr)
    library(ggforce)
    library(treeio)
    library(magick)
    setwd("validation/plotting")
}


{
    aliased_results <- read.csv("data/results_alias.csv")
    rfam_regions <- read.csv("data/all_regions_14.10.csv")
}

# Colour palattes
{
    decision_pal <- c("pass" = "#6956e7", "fail" = "#eaef5b", "incon" = "#000000")
    active_pal <- c("0" = "black", "1" = "#bb9601")
    active_pie_pal <- c("0" = "#8b8b8b", "1" = "#d9b937")
}

# themes
{
    def_theme <- theme(
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


# Funcs
{
    tree_y <- function(ggtree, data) {
        if (!inherits(ggtree, "ggtree")) {
            stop("not a ggtree object")
        }
        left_join(select(data, label), select(ggtree$data, label, y)) %>%
            pull(y)
    }

    # overwrite the default expand for continuous scales
    scale_y_tree <- function(expand = expand_scale(0, 0.6), ...) {
        scale_y_continuous(expand = expand, ...)
    }

    # get the range of the ggtree y-axis data
    tree_ylim <- function(ggtree) {
        if (!inherits(ggtree, "ggtree")) {
            stop("not a ggtree object")
        }
        range(ggtree$data$y)
    }


    # plot data next to a ggtree aligned by shared labels
    ggtreeplot <- function(ggtree, data = NULL, mapping = aes(), flip = FALSE,
                           expand_limits = expand_scale(0, .6), ...) {
        if (!inherits(ggtree, "ggtree")) {
            stop("not a ggtree object")
        }

        # match the tree limits
        limits <- tree_ylim(ggtree)
        limits[1] <- limits[1] + (limits[1] * expand_limits[1]) - expand_limits[2]
        limits[2] <- limits[2] + (limits[2] * expand_limits[3]) + expand_limits[4]

        if (flip) {
            mapping <- modifyList(aes_(x = ~x), mapping)
            data <- mutate(data, x = tree_y(ggtree, data))
            gg <- ggplot(data = data, mapping = mapping, ...) +
                scale_x_continuous(limits = limits, expand = c(0, 0))
        } else {
            mapping <- modifyList(aes_(y = ~y), mapping)
            data <- mutate(data, y = tree_y(ggtree, data))
            gg <- ggplot(data = data, mapping = mapping, ...) +
                scale_y_continuous(limits = limits, expand = c(0, 0))
        }
        gg
    }

    # get rid of superfluous axis - this works after coord_flip, so it also works
    # for the rotated histogram
    no_y_axis <- function() {
        theme(
            axis.line.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
    }

    pie_df <- function(d, r, incr) {
        ssize <- 10000
        x <- c()
        y <- c()
        idvec <- c()
        fval <- c()

        cntr <- 1

        for (i in seq_along(d[, 1])) {
            dd <- d[i, ]
            segs <- dd$total
            act <- dd$active
            x0 <- dd$x_
            y0 <- dd$y_

            if (segs < 1) next
            theta <- 2 * pi / segs

            # adjust incr to ensure are no unfilled slices
            incr2 <- incr + ((theta %% incr) / (theta %/% incr))

            for (idx in seq_along(1:segs)) {
                if (segs > 1) {
                    x <- c(x, x0)
                    y <- c(y, y0)
                    fval <- c(fval, as.numeric(idx <= act))
                    idvec <- c(idvec, toString(interaction(dd$microbe, dd$target_name, as.factor(idx))))

                    cntr <- cntr + 1
                }

                angvec <- seq.int(from = theta * (idx - 1), to = (theta * idx), by = incr2)
                x <- c(x, x0 + r * sin(angvec))
                y <- c(y, y0 + r * cos(angvec))
                fval <- c(fval, rep(as.numeric(idx <= act), length(angvec)))
                idvec <- c(idvec, rep(toString(interaction(dd$microbe, dd$target_name, as.factor(idx))), length(angvec)))
            }
        }
        fval <- as.factor(fval)
        polydf <- na.omit(data.frame(x, y, idvec, fval, stringsAsFactors = TRUE))

        return(polydf)
    }
}

{
    {
        tree <- read.tree("data/big_fig/taxa_list.tree")
        phylo_dict <- fromJSON(file = "data/big_fig/fancy_tips.json")

        to_drop <- c(
            "'s__Pseudomonas_E fluorescens'",
            "'s__Piscirickettsia salmonis'",
            "'s__Acinetobacter baumannii'"
        )

        tree <- drop.tip(tree, to_drop)
        tax_tree <- ggtree(tree, size = 0.8, layout = "rectangular", ladderize = FALSE) +
            # ggtitle("Validation set tax tree") +
            geom_tiplab(size = 7, linesize = .5, align = TRUE, as_ylab = FALSE) +
            xlim_tree(2.5)

        for (i in seq_along(tree$tip.label)) {
            tax_tree$data$label[i] <- phylo_dict[[tax_tree$data$label[i]]][3]
        }

        tax_tree
    }


    {
        lst <- list()
        for (i in phylo_dict) {
            lst[[i[2]]] <- i[3]
        }

        inf_df <- read.csv("data/big_fig/locus_inferences_sum.csv")
        inf_df$inactive <- as.numeric(inf_df$total) - as.numeric(inf_df$active)
        nl <- lapply(inf_df$microbe, function(x) lst[[x]])
        inf_df$label <- unlist(nl)
        # inf_df$label <- factor(inf_df$label)

        inf_df$microbe <- factor(inf_df$microbe)
        inf_df$target_name <- factor(inf_df$target_name)

        # lst[inf_df$microbe]
        # inf_df$label <- inf_df$microbe
        inf_df$x_ <- match(inf_df$target_name, levels(inf_df$target_name))
        inf_df$y_ <- tree_y(tax_tree, inf_df)
        inf_df$microbe <- fct_reorder(inf_df$microbe, inf_df$y_, min)

        polydf <- pie_df(inf_df, r = 0.4, incr = 2 * pi / 120)
        head(inf_df)
    }


    {
        test_pie <- ggtreeplot(tax_tree, inf_df, aes(x = x_)) +
            def_theme +
            theme(
                axis.text.x = element_text(size = 13, angle = 70, hjust = 1, colour = "black"),
                axis.ticks.y = element_blank(),
                # axis.text.y = element_text(size = 13, hjust = 1, colour = "black"),
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
                axis.line = element_blank(),
                # plot.margin = unit(c(0, 0, 0, 0), "inches"),
                panel.background = element_blank(),
                legend.position = "top",
            ) +
            geom_polygon(
                data = polydf,
                mapping = aes(
                    x = x, y = y, fill = fval, group = idvec, color = "a"
                    # linetype = "solid",# linewidth = 0.5
                ),
                linewidth = 0.1
            ) +
            geom_text(
                aes(
                    # x = floor(max(polydf$x)) + 1,
                    x = -1,
                    label = sras,
                    hjust = 1
                ),
                size = 7
            ) +
            scale_colour_manual(name = "d", values = c("black")) +
            scale_fill_manual(name = "Transcriptional activity:", labels = c("no evidence", "condition-dependent termination efficiency"), values = active_pie_pal) +
            coord_fixed(clip = "off") +
            scale_x_continuous(breaks = c(seq_along(levels(inf_df$target_name))), labels = levels(inf_df$target_name)) +
            # scale_y_continuous(breaks = c(seq_along(levels(inf_df$microbe))), labels = levels(inf_df$microbe)) +
            guides(colour = "none")

        test_pie
    }

    # head(ggplot_build(test_pie)$plot)
    # head(ggplot_build(test_pie)$data[[1]])
    # head(ggplot_build(test_pie)$data[[2]])

    tag_y <- 0.98
    patched <- tax_tree +
        theme(plot.tag.position = c(0, tag_y)) +
        test_pie +
        theme(plot.tag.position = c(0, tag_y)) +
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(size = 30, face = "bold"))
    patched

    # match(df$microbe, levels(df$microbe))
    # length(levels(inf_df$target_name))
}


{
    alph <- 0.6
    save_path <- str_glue("plots/big_fig_3.png")
    ggsave(save_path, plot = patched, dpi = 320 * alph, units = "px", width = 7000 * alph, height = 4000 * alph)
}

# -----------------------------------------------------------------------
# Conceptual fragmentation and convolution application figure
{
    frag_img <- image_read_svg("./plots/concept_figs/__conceptual_fragmentation.svg", width = 3000)
    conv_img <- image_read_svg("./plots/concept_figs/convolution_application.svg", width = 3000)

    frag_plot <- image_ggplot(frag_img)
    conv_plot <- image_ggplot(conv_img)

    patched <- frag_plot / conv_plot +
        plot_annotation(tag_levels = "a") &
        theme(plot.tag = element_text(size = 30, face = "bold"))

    patched
}


{
    alph <- 0.6
    save_path <- str_glue("plots/concept_fragmentation.png")
    ggsave(save_path, plot = patched, dpi = 320 * alph, units = "px", width = 7000 * alph, height = 4000 * alph)
}

