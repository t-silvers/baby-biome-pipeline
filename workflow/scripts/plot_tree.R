.libPaths(c(file.path(Sys.getenv("GROUP_HOME"), "envs/renv/ggtree/library"), .libPaths()))

library(dplyr)
library(ggnewscale)
library(ggtree)
library(optparse)
library(readr)

option_list <- list(
    optparse::make_option(c("-t", "--tree")),
    optparse::make_option(c("-m", "--metadata")),
    optparse::make_option(c("-o", "--output"))
    )

# option_list <- list(
#     optparse::make_option(c("-t", "--tree"), default = "/ptmp/thosi/wv-multilib/results/trees/B002-Bifidobacterium_longum.raxml.bestTree"),
#     optparse::make_option(c("-m", "--metadata"), default = "/ptmp/thosi/wv-multilib/results/trees/B002-Bifidobacterium_longum/metadata.txt"),
#     optparse::make_option(c("-o", "--output"), default = "/ptmp/thosi/wv-multilib/results/trees/B002-Bifidobacterium_longum/phylo.png")
#     )

args <- optparse::parse_args(
    optparse::OptionParser(option_list=option_list)
    )

tree.hat <- ggtree::read.tree(args$tree)


# Parse metadata

metadata <- readr::read_tsv(args$metadata)

metadata.plot <- 
    metadata %>%
    dplyr::filter(id %in% tree.hat$tip.label) %>%
    dplyr::arrange(match(id, tree.hat$tip.label)) %>%
    tibble::column_to_rownames(var = "id") %>%
    dplyr::mutate(relationship = dplyr::case_when(
        relationship == 'M' ~ 'Mother',
        relationship == 'F' ~ 'Father',
        relationship == 'S' ~ 'Sibling',
        relationship == 'B' ~ 'Baby',
        TRUE ~ NA)
        ) %>%
    dplyr::select(relationship, timepoint_days)


# Plot

dpi <- 72
plot.dims <- c(800, 100 + nrow(metadata.plot) * 10)
scaling <- 2
hmap.x <- .3

options(
    bitmapType = "cairo",
    repr.plot.width = plot.dims[1] / dpi,
    repr.plot.height = plot.dims[2] / dpi
    )

png(
    args$output,
    type="cairo",
    width = scaling * plot.dims[1],
    height = scaling * plot.dims[2],
    res = scaling * dpi
    )

phylo.plot <-
    ggtree::ggtree(tree.hat) + 
    ggtree::geom_tiplab(align = T)

# msaplot(phylo.plot, "/viper/ptmp/thosi/wv-multilib/results/trees/B002-Bifidobacterium_longum/msa-filtered.fas", offset=3, width=2)

# phylo.plot <-
#     ggtree::ggtree(
#         tree.hat,
#         branch.length = "none",
#         layout = "circular"
#         ) +
#     ggtree::geom_tiplab()

phylo.plot <- 
    ggtree::gheatmap(
        phylo.plot,
        dplyr::select(metadata.plot, timepoint_days),
        offset = hmap.x,
        width = .025,
        colnames = FALSE,
        color = "black",
        ) + 
    ggplot2::scale_fill_gradient(limits = c(-1, 300), low = "white", high = "black", name = "Time Point\n(days post-delivery)")

phylo.plot <- 
    phylo.plot + 
    ggnewscale::new_scale_fill()

phylo.plot <-
    ggtree::gheatmap(
        phylo.plot,
        dplyr::select(metadata.plot, relationship),
        offset = hmap.x + .05,
        width = .025,
        colnames = FALSE,
        color = "black",
        ) + 
    ggplot2::scale_fill_manual(
        values = c("Mother" = "gold1", "Father" = "royalblue1", "Baby" = "seagreen1", "Sibling" = "seagreen4"),
        name = "Relationship",
        )

phylo.plot

dev.off()