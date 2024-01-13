install.packages("nycflights13")
install.packages("tidyverse")

library(tidyverse)
library(nycflights13)

setwd("validation/plotting")

cps <- quantile(flights$dep_delay, probs = c(0.25, 0.5, 0.9), na.rm = TRUE)

df <- flights |>
  mutate(dep_cut = cut(dep_delay, breaks = cps)) |>
  mutate(
    dep_cut = fct_na_value_to_level(dep_cut, level = "Missing"),
    dep_cut = fct_recode(dep_cut, Low = "(-5,-2]", Mid = "(-2,49]"),
    carrier = factor(carrier),
    month = factor(month)
  ) |>
  select(month, carrier, dep_cut)

pal <- c("pass" = "#6956e7", "fail" = "#eaef5b", "incon" = "#000000")

{
  dset <- "s_sanguinis"
  df <- read.csv("data/results_alias.csv")
  df <- df[df[, "dataset"] == dset, ]

  df |>
    count(rowid_alias, condition_alias, decision) |>
    mutate(prop = n / sum(n), .by = c(rowid_alias, condition_alias)) |>
    ggplot(aes(x = 0, y = prop, fill = decision)) +
    geom_col(position = "stack") +
    scale_fill_manual(name = "decision", values = pal) +
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
  save_path <- sprintf("plots/%s.png", dset)
  ggsave(save_path, dpi = 320, units = "px")
}
