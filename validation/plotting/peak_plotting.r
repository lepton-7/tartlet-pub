{
    library(tidyverse)
    library(ggplot2)
    library(reshape2)
    library(stringr)
    library(ggforce)
    setwd("validation/plotting")
}

# themes
{
    def_theme <- theme_classic() +
        theme(
            axis.text.x = element_text(size = 13, angle = 0, colour = "black"),
            axis.text.y = element_text(size = 13, face = "bold"),
            axis.title = element_text(size = 18, face = "bold"),
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

{
    dset <- "b_sub_168"
    # dset <- "Acidimicrobiia"
    df <- read.csv(str_glue("../alignment/outputs/{dset}/plots/peak_log.csv"))
    statdf <- read.csv(str_glue("../alignment/outputs/{dset}/plots/cluster_stats.csv"))

    # ppath <- str_glue("D:/School stuff/Bagby Lab/Projects/Riboswitches stuff/Code/Genomes/EMERGE_Sharing/tart/alignment/outputs/{dset}/plots")
    # df <- read.csv(str_glue("{ppath}/peak_log.csv"))
    # statdf <- read.csv(str_glue("{ppath}/cluster_stats.csv"))

    df <- df[df$decision != "", ]

    df$rowid <- sapply(df$rowid, rowidlabeller)
    statdf$rowid <- sapply(statdf$rowid, rowidlabeller)

    statdf$ann <- statdf$delta_mean < 0 &
        statdf$delta_mean_pval < 0.01 &
        statdf$delta_variance_pval < 0.01 &
        statdf$delta_variance > statdf$noiseset_delta_variance
    statdf <- statdf[statdf$ann, ]
    # statdf$ann_text <- "*✝"
    if (nrow(statdf) > 0) statdf$ann_text <- "*"
    # df$pval <- -log(df$coverage_delta_pval)
    # df$pval[df$pval > 10] <- 10



    peakplot <- ggplot(df, aes(x = from_riboswith_end_relative, y = coverage_delta_stable_relative)) +
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
            y = "Stable fractional coverage change",
            x = "Position relative to riboswitch 3′ end as a fraction of riboswitch length"
        ) +
        def_theme +
        theme(
            plot.title = element_text(size = 14),
            legend.position = "none"
        ) +
        coord_cartesian(ylim = c(-1.2, 0.5)) +
        guides(alpha = "none", color = "none")

    if (nrow(statdf) > 0) {
        peakplot <- peakplot + geom_text(data = statdf, mapping = aes(x = pos_mean, y = 0.3, label = ann_text), size = 12, fontface = "bold")
    }

    peakplot
}

{
    alph <- 2.5
    save_path <- str_glue("plots/peak_plots/{dset}_peak_plot.png")
    ggsave(save_path, dpi = 320 * alph, units = "px", width = 7000 * alph, height = 4000 * alph)
}
