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

# Helper for geom: return a polygonGlob for each pie slice
# For now this creates triangles instead of radial segments
{
    slicer <- function(coor, dat, idx) {
        if (dat$segments < 2) {
            return(grid::nullGrob())
        }
        theta <- 2 * pi / dat$segments
        x0 <- coor$x
        y0 <- coor$y
        r <- coor$r

        # x <- c(x0)
        # print(r)
        # y <- c(y0)

        px <- c(x0, x0 + r * sin(theta * (idx - 1)), x0 + r * sin(theta * idx))
        py <- c(y0, x0 + r * cos(theta * (idx - 1)), x0 + r * cos(theta * idx))

        # append(x, px)
        # append(y, py)
        # print(py)
        grid::polygonGrob(
            x = px, y = py,
            # id.lengths = c(3),
            default.units = "native",
            gp = grid::gpar(
                col = coor$colour,
                fill = scales::alpha(coor$fill, coor$alpha),
                lwd = coor$linewidth * .pt,
                lty = coor$linetype
            )
        )
    }

    # xslicer <- function(coor, dat)


    # Custom geom lmao
    GeomPie <- ggproto("GeomPie", Geom,
        required_aes = c("x", "y", "r", "segments"),
        default_aes = aes(
            colour = NA, fill = NA, linewidth = 0.5,
            linetype = 1, alpha = 1
        ),
        draw_key = draw_key_polygon,
        draw_panel = function(data, panel_params, coord) {
            coords <- coord$transform(data, panel_params)
            # if (data$segments < 2) {
            #     return(grid::nullGrob())
            # }
            for (j in seq_along(data[, 1])) {
                for (i in seq_along(1:data[j, ]$segments)) {
                    # print(data[j, ])
                    # print(i)
                    slicer(coords[j, ], data[j, ], i)
                }
            }
        }
    )
    geom_pie <- function(mapping = NULL, data = NULL, stat = "identity",
                         position = "identity", na.rm = FALSE, show.legend = NA,
                         inherit.aes = TRUE, ...) {
        layer(
            geom = GeomPie, mapping = mapping, data = data, stat = stat,
            position = position, show.legend = show.legend, inherit.aes = inherit.aes,
            params = list(na.rm = na.rm, ...)
        )
    }
}

# Funcs
{
    tree_y <- function(ggtree, data) {
        if (!inherits(ggtree, "ggtree")) {
            stop("not a ggtree object")
        }
        left_join(select(data, label), select(ggtree$data, label, y)) %>%
            pull(y)
    }

    pie_df <- function(d, r, incr) {
        ssize <- 10000
        # x <- rep(NA, ssize)
        # y <- rep(NA, ssize)
        # idvec <- rep("PP", ssize)
        # fval <- rep(0, ssize)
        x <- c()
        y <- c()
        idvec <- c()
        fval <- c()

        cntr <- 1

        for (i in seq_along(d[, 1])) {
            dd <- d[i, ]
            # print(d[i, ])
            segs <- dd$total
            act <- dd$active
            x0 <- dd$x_
            y0 <- dd$y_

            if (segs < 1) next
            theta <- 2 * pi / segs

            # adjust incr to ensure are no unfilled slices
            incr2 <- incr + ((theta %% incr) / (theta %/% incr))

            for (idx in seq_along(1:segs)) {
                # print(idx)
                if (segs > 1) {
                    x <- c(x, x0)
                    y <- c(y, y0)
                    fval <- c(fval, as.numeric(idx <= act))
                    idvec <- c(idvec, toString(interaction(dd$microbe, dd$target_name, as.factor(idx))))
                    # print(idvec[cntr])
                    # print(interaction(dd$microbe, dd$target_name, as.factor(idx)))

                    cntr <- cntr + 1
                }
                # if (segs > 1) {
                #     x[cntr] <- x0
                #     y[cntr] <- y0
                #     fval[cntr] <- as.numeric(idx <= act)
                #     idvec[cntr] <- toString(interaction(dd$microbe, dd$target_name, as.factor(idx)))
                #     # print(idvec[cntr])
                #     # print(interaction(dd$microbe, dd$target_name, as.factor(idx)))

                #     cntr <- cntr + 1
                # }

                angvec <- seq.int(from = theta * (idx - 1), to = (theta * idx), by = incr2)
                x <- c(x, x0 + r * sin(angvec))
                y <- c(y, y0 + r * cos(angvec))
                fval <- c(fval, rep(as.numeric(idx <= act), length(angvec)))
                idvec <- c(idvec, rep(toString(interaction(dd$microbe, dd$target_name, as.factor(idx))), length(angvec)))

                # for (ang in seq.int(from = theta * (idx - 1), to = (theta * idx), by = incr2)) {
                #     x[cntr] <- x0 + r * sin(ang)
                #     y[cntr] <- y0 + r * cos(ang)
                #     fval[cntr] <- as.numeric(idx <= act)
                #     idvec[cntr] <- toString(interaction(dd$microbe, dd$target_name, as.factor(idx)))

                #     cntr <- cntr + 1
                # }
            }
        }
        fval <- as.factor(fval)
        polydf <- na.omit(data.frame(x, y, idvec, fval, stringsAsFactors = TRUE))

        return(polydf)
    }
    polydf <- pie_df(inf_df, r = 0.4, incr = 2 * pi / 60)

    # incr <- 0.10472 + ((0.571198 %% 0.10472) / (0.571198 %/% 0.10472))
    # incr <- 2 * pi / 60
    # theta <- 2 * pi / 11
    # incr2 <- incr + ((theta %% incr) / (theta %/% incr))
}
# seq.int(from = 0.1 * (5 - 1), to = 0.1 * 5, by = 0.013)

{
    {
        tree <- read.tree("data/big_fig/taxa_list.tree")
        phylo_dict <- fromJSON(file = "data/big_fig/fancy_tips.json")

        tree <- drop.tip(tree, c("'s__Piscirickettsia salmonis'", "'s__Acinetobacter baumannii'"))
        tax_tree <- ggtree(tree, size = 0.8, layout = "rectangular", ladderize = FALSE) +
            # ggtitle("Validation set tax tree") +
            geom_tiplab(size = 8, linesize = .5, align = TRUE, as_ylab = FALSE) +
            xlim_tree(2.5)

        for (i in seq_along(tree$tip.label)) {
            tax_tree$data$label[i] <- phylo_dict[[tax_tree$data$label[i]]][2]
        }

        tax_tree
    }

    {
        df <- read.csv("data/big_fig/locus_inferences.csv")
        df$is_active <- as.factor(df$is_active)

        df$microbe <- factor(df$microbe)
        df$target_name <- factor(df$target_name)
        df$label <- df$microbe

        pie_mat <- ggplot(df, aes(x = target_name, y = tree_y(tax_tree, df), colour = is_active)) +
            def_theme +
            theme(
                axis.text.x = element_text(size = 13, angle = 70, hjust = 1, colour = "black")
            ) +
            scale_colour_manual(name = "Transcriptionally active", values = active_pal) +
            geom_point(
                size = 2,
                position = position_jitter(width = 0.2, height = 0.2),
            ) +
            guides(colour = "none")

        pie_mat
        # ggplot2::coord_fixed() +
        # geom_circle(aes(x0 = 6, y0 = 9, r = 5, fill = 1))

        # levels(factor(df$target_name))
        # head(pie_mat$data)
    }

    {
        inf_df <- read.csv("data/big_fig/locus_inferences_sum.csv")
        inf_df$inactive <- as.numeric(inf_df$total) - as.numeric(inf_df$active)

        inf_df$microbe <- factor(inf_df$microbe)
        inf_df$target_name <- factor(inf_df$target_name)

        inf_df$label <- inf_df$microbe
        inf_df$x_ <- match(inf_df$target_name, levels(inf_df$target_name))
        inf_df$y_ <- tree_y(tax_tree, inf_df)
        inf_df$microbe <- fct_reorder(inf_df$microbe, inf_df$y_, min)

        polydf <- pie_df(inf_df, r = 0.4, incr = 2 * pi / 120)
    }
    {
        test_pie <- ggplot(inf_df, aes(x = x_, y = y_)) +
            def_theme +
            theme(
                axis.text.x = element_text(size = 13, angle = 70, hjust = 1, colour = "black"),
                # axis.ticks.y = element_blank(),
                axis.text.y = element_text(size = 13, hjust = 1, colour = "black"),
                # axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                axis.title.x = element_blank(),
            ) +
            geom_polygon(
                data = polydf,
                mapping = aes(
                    x = x, y = y, fill = fval, group = idvec, color = "a"
                    # linetype = "solid",# linewidth = 0.5
                ),
                linewidth = 0.1
            ) +
            scale_colour_manual(name = "d", values = c("black")) +
            scale_fill_manual(name = "Active?", values = active_pie_pal) +
            # geom_circle(mapping = aes(
            #     x0 = x_,
            #     y0 = y_,
            #     # fill = active
            #     # linewidth = 0.5
            #     r = 0.4
            # )) +
            coord_fixed() +
            # geom_pie(mapping = aes(r = 0.04, segments = total, colour = "black")) +
            scale_x_continuous(breaks = c(seq_along(levels(inf_df$target_name))), labels = levels(inf_df$target_name)) +
            scale_y_continuous(breaks = c(seq_along(levels(inf_df$microbe))), labels = levels(inf_df$microbe)) +
            guides(colour = "none")

        test_pie
    }


    # levels(as.factor(inf_df$y_))

    # head(ggplot_build(test_pie)$plot)
    # head(ggplot_build(test_pie)$data[[1]])
    # head(ggplot_build(test_pie)$data[[2]])

    # patched <- tax_tree + pie_mat + test_pie
    patched <- tax_tree + test_pie
    patched

    # match(df$microbe, levels(df$microbe))
    # length(levels(inf_df$target_name))
}

{
    alph <- 0.6
    save_path <- str_glue("plots/big_fig.png")
    ggsave(save_path, dpi = 320 * alph, units = "px", width = 7000 * alph, height = 3500 * alph)
}
