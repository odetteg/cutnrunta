library(ggplot2)
library(tidyverse)
library(cowplot)
library(RColorBrewer)

frag_data <- read.csv(snakemake@input[[1]], header = TRUE)

frag_data$condition <- str_replace(frag_data$sample, "-[0-9].*$", "")
frag_data$replicate <- str_extract(frag_data$sample, "(?<=-)[0-9]+")
frag_data$genotype <- str_extract(frag_data$sample, "^(T)?BL3")
frag_data$mark <- str_replace(frag_data$sample, "^(T)?BL3-", "") %>%
                str_replace("_S[0-9].*$", "")


# plot_data <- frag_data %>% filter(size <= 500)

frag_ladder_plot <- ggplot(frag_data, aes(x = size, y = count, color = condition, linetype = replicate)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    scale_color_brewer(palette = "Set1") +
    theme_cowplot() +
    labs(
        title = "Fragment Length Distribution",
        subtitle = "CUT&RUN Nucleosomal Ladder",
        x = "Fragment Size (bp)",
        y = "Read Count",
        color = "Condition",
        linetype = "Replicate"
    ) +
    geom_vline(xintercept = c(120, 150), linetype = "dotted", color = "gray40")

frag_data <- frag_data %>%
  mutate(type = ifelse(grepl("IgG", sample, ignore.case = TRUE), 
                       "Control (IgG)", 
                       "Target Marks"))
facet_plot_data <- frag_data %>%
  filter(!(type == "Target Marks" & grepl("IgG", mark)))

faceted_plot <- ggplot(facet_plot_data, 
                       aes(x = size, y = count, color = mark, linetype = replicate)) +
    geom_line(linewidth = 0.8, alpha = 0.8) +
    facet_grid(type ~ genotype, scales = "free") + 
    scale_color_brewer(palette = "Set1") + 
    theme_cowplot() +
    labs(
        title = "CUT&RUN Fragment Size Periodicity",
        subtitle = "Comparing Histone Marks against IgG Background",
        x = "Fragment Size (bp)",
        y = "Read Count"
    ) +
    geom_vline(xintercept = c(120, 150), linetype = "dotted", color = "black", alpha=0.3) +
    theme(panel.spacing = unit(1, "lines"),
          strip.background = element_rect(fill = "gray90"))

ggsave(snakemake@output[["frag_len_dst_pdf"]], 
       plot = frag_ladder_plot, width = 10, height = 5)

ggsave(snakemake@output[["frag_len_dst_facet_pdf"]], 
       plot = faceted_plot, width = 12, height = 8)