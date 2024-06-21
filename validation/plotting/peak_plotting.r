{
    library(tidyverse)
    library(ggplot2)
    library(reshape2)
    library(stringr)
    library(ggforce)
    library(patchwork)
    setwd("validation/plotting")
}

# themes
{
    def_theme <- theme_classic() +
        theme(
            axis.text.x = element_text(size = 13, angle = 0, colour = "black"),
            axis.text.y = element_text(size = 13, colour = "black"),
            axis.title = element_text(size = 18),
            # plot.title = element_blank()
        )
}

{
    rowidlabeller <- function(variable) {
        splist <- strsplit(variable, split = "#")[[1]]
        rswtch <- splist[1]

        magname <- strsplit(splist[2], split = "_")[[1]]
        magname2 <- paste(magname[1:length(magname) - 1], collapse = "_")
        if (magname2 != "") {
            tt <- sprintf("%s|%s|%s%s", rswtch, magname2, splist[3], splist[length(splist)])
        } else {
            tt <- sprintf("%s|%s%s", rswtch, splist[3], splist[length(splist)])
        }
        return(tt)
    }
}

peakplotmaker <- function(dset) {
    # dset <- "b_sub_168"
    # dset <- "Acidimicrobiia"
    df <- read.csv(str_glue("../alignment/outputs/{dset}/p2/peak_log.csv"))
    statdf <- read.csv(str_glue("../alignment/outputs/{dset}/p2/cluster_stats.csv"))

    # ppath <- str_glue("D:/School stuff/Bagby Lab/Projects/Riboswitches stuff/Code/Genomes/EMERGE_Sharing/tart/alignment/outputs/{dset}/plots")
    # df <- read.csv(str_glue("{ppath}/peak_log.csv"))
    # statdf <- read.csv(str_glue("{ppath}/cluster_stats.csv"))

    df <- df[df$decision != "", ]

    df$rowid <- sapply(df$rowid, rowidlabeller)
    statdf$rowid <- sapply(statdf$rowid, rowidlabeller)

    statdf$ann <- statdf$delta_mean < 0 &
        statdf$delta_mean_pval < 0.01 &
        statdf$delta_variance_pval < 0.05 &
        statdf$delta_variance > statdf$noiseset_delta_variance
    statdf <- statdf[statdf$ann, ]
    # statdf$ann_text <- "*✝"
    if (nrow(statdf) > 0) statdf$ann_text <- "*"
    # df$pval <- -log(df$coverage_delta_pval)
    # df$pval[df$pval > 10] <- 10



    peakplot <- ggplot(df, aes(x = from_riboswith_end_relative, y = coverage_delta_stable_relative)) +
        geom_vline(xintercept = 0, colour = "#959595", linetype = "dashed") +
        # peakplot <- ggplot(df, aes(x = from_riboswith_end_relative, y = 0)) +
        # peakplot <- ggplot(df, aes(x = from_riboswith_end_relative, y = coverage_delta_relative)) +
        geom_line(aes(color = transcriptome)) +
        geom_point() +
        # facet_wrap(~rowid, scales = "free_y") +
        facet_wrap(~rowid, scales = "fixed") +
        # ggforce::geom_mark_hull(aes(group = cluster)) +
        geom_mark_ellipse(aes(fill = as.factor(cluster)), expand = unit(0.5, "mm")) +
        labs(
            title = str_glue("{dset} riboswitches across all available conditions"),
            y = "Fractional coverage change",
            x = "Position relative to riboswitch 3′ end as a fraction of riboswitch length"
        ) +
        def_theme +
        theme(
            plot.title = element_text(size = 14),
            legend.position = "none",
            strip.text = element_text(size = 12),
            strip.background = element_rect(fill = "#cdcdcd", linewidth = 0),
        ) +
        coord_cartesian(ylim = c(-1.2, 0.5)) +
        guides(alpha = "none", color = "none")

    if (nrow(statdf) > 0) {
        peakplot <- peakplot + geom_text(data = statdf, mapping = aes(x = pos_mean, y = 0.3, label = ann_text), size = 12, fontface = "bold")
    }

    # peakplot

    print(str_glue("Saving plot: {dset}"))
    alph <- 0.7
    save_path <- str_glue("plots/peak_plots/{dset}_peak_plot_3.png")
    ggsave(save_path, plot = peakplot, dpi = 320 * alph, units = "px", width = 7000 * alph, height = 4000 * alph)
}

{
    dsets <- c(
        "a_fischeri_ES114",
        "a_kunk",
        "b_anth",
        "b_frag",
        "b_pseudo",
        "b_sub_168",
        "b_theta",
        "b_xyla",
        "c_diff",
        "c_vibrioides",
        "d_vulg",
        "e_coli",
        "e_fae",
        "e_limo",
        "k_pneum",
        "m_smeg",
        "m_tuber",
        "n_gonorr",
        "p_cholor_aureo3084",
        # "p_fluor",
        "s_coelicolor",
        "s_elon",
        "s_enter_typh",
        "s_epi",
        "s_meli",
        "s_sanguinis",
        "s_spcc6803",
        "x_albi",
        "x_ory"
    )
    # These were removed from ^
    #   "c_basil"
    #   "p_aeru
    #   "p_salmo" <- This one had no peak log?
}


# plot suffix:
# 1. original plots (both < 0.01)
# 2. new var stat (both < 0.01)
# 3. new var stat (mean < 0.01, var < 0.05)

sapply(dsets, peakplotmaker)
peakplotmaker("b_sub_168")
