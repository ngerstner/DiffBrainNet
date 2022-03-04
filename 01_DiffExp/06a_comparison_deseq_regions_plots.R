##################################################
## Project: DexStim Mouse Brain
## Date: 21.01.2021
## Author: Nathalie
##################################################
# Compare deseq SV DE genes between regions
# Plots for SAB meeting

setwd("~/Documents/ownCloud/DexStim_RNAseq_Mouse")

# library(rlist)
library(RColorBrewer)
# library(org.Mm.eg.db)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(gridExtra)
library(ggraph)
library(igraph)


regions <-
  c("AMY", "PFC", "PVN", "CER", "vDG", "dDG", "vCA1", "dCA1")

folder_plots <- paste0("figures")
folder_tables <- paste0("tables")


# 1. Read DE tables from all regions ----------

list_reg_sig <- list()
list_genes_sig <- list()

for (reg in regions) {
  res <-
    fread(
      file = paste0(
        folder_tables,
        "/02_",
        reg,
        "_deseq2_Dex_1_vs_0_lfcShrink.txt"
      ),
      sep = "\t"
    )
  na_indices <- which(is.na(res$padj))
  res$padj[na_indices] <- 1
  res_sig <- res[res$padj <= 0.1, ]
  # res_sig <- res[res$log2FoldChange >= 1]
  list_reg_sig[[reg]] <- res_sig
  list_genes_sig[[reg]] <- rownames(res_sig)
}



# 2. Concatenate DE tables -----------------

data <- bind_rows(list_reg_sig, .id = "region") %>%
  group_by(Ensembl_ID) %>%
  summarise(region = list(region))

data_unique <- data %>%
  mutate(nr_regions = lengths(region)) %>%
  mutate(unique = (nr_regions == 1)) %>%
  unnest(cols = c(region)) %>%
  mutate("combined_id" = paste0(region, "-", Ensembl_ID))

data_barplot <- data_unique %>%
  group_by(region, unique) %>%
  count() %>%
  group_by(region) %>%
  mutate(sum = sum(n))

matrix_heatmap <- matrix(
  data = 0,
  nrow = length(regions),
  ncol = length(regions),
  dimnames = list(regions, regions)
)
for (reg1 in regions) {
  for (reg2 in regions) {
    if (reg1 != reg2) {
      nr_genes <-
        length(which(
          sapply(data$region, function(x)
            reg1 %in% x) &
            sapply(data$region, function(x)
              reg2 %in% x)
        ))
      matrix_heatmap[reg1, reg2] <- nr_genes
    } else{
      nr_genes <-
        length(which(sapply(data$region, function(x) {
          (reg1 %in% x) & (length(x) == 1)
        })))
      matrix_heatmap[reg1, reg2] <- nr_genes
    }
  }
}
matrix_heatmap <- data.frame(matrix_heatmap)
matrix_heatmap$x.values <- rownames(matrix_heatmap)
matrix_heatmap.melted <-
  melt(matrix_heatmap, id.vars = c("x.values"))
matrix_heatmap.melted <-
  pivot_longer(
    matrix_heatmap,
    cols = AMY:dCA1,
    names_to = c("variable"),
    names_transform = list(variable = as.factor)
  )
matrix_heatmap.melted$x.values <-
  factor(x = matrix_heatmap.melted$x.values,
         levels = levels(matrix_heatmap.melted$variable))



# 3. Stacked barplot -------------------------

bp <- ggplot(data_barplot, aes(fill = unique, y = n, x = region)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(
    name = "",
    labels = c("DE in multiple regions", "DE unique"),
    values = c("tan1", "royalblue")
  ) +
  xlab("brain region") +
  ylab("# diff. exp. genes") +
  theme_light() +
  theme(
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    # legend.title = element_blank(),
    legend.position = "top"
  ) +
  geom_text(aes(label = paste0(round((n / sum) * 100, digits = 1
  ), "%")),
  position = position_stack(vjust = 0.5),
  size = 4)
ggsave(
  bp,
  filename = paste0(
    folder_plots,
    "/06_comparison_deseq_barplot_vertical_percentage.png"
  ),
  width = 8,
  height = 6
)

bp_horizontal <- bp +
  # scale_x_reverse() +
  coord_flip()

ggsave(
  bp_horizontal,
  filename = paste0(
    folder_plots,
    "/06_comparison_deseq_barplot_horizontal_percentage.png"
  ),
  width = 6,
  height = 8
)



# 4. Heatmap ---------------------------------

# Prepare heatmap
# Would it make sense to have total nr of DE genes in diagonal?
gm <-
  ggplot(data = matrix_heatmap.melted, aes(
    x = factor(x.values),
    y = variable,
    fill = value
  )) +
  geom_tile() +
  scale_fill_distiller(
    name = "Nr. of DE genes",
    palette = "Reds",
    direction = 1,
    na.value = "transparent"
  ) +
  theme_light() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Prepare x axis barplot
bp.x <-
  ggplot(data = data_barplot, aes(
    x = factor(region, levels = levels(matrix_heatmap.melted$x.values)),
    y = n,
    fill = unique
  )) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(
    name = "",
    labels = c("DE in multiple regions", "DE unique"),
    values = c("red3", "navy")
  ) +
  ylab("") +
  theme_light() +
  theme(
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.position = "top",
    legend.text = element_text(size = 10)
  ) +
  scale_x_discrete(position = "top")

# Put plots together
hm_comb <- grid.arrange(bp.x, gm, nrow = 2, ncol = 1)

ggsave(
  hm_comb,
  filename = paste0(folder_plots, "/06_comparison_deseq_heatmap.png"),
  height = 10,
  width = 7
)

# # 5. Circular plot -------------------------
#
# # create data frame giving hierarchical structure of regions and DE genes
# d1 <- data.frame(from="origin", to=regions)
# d2 <- data.frame(from=data_unique$region, to=data_unique$combined_id)
# edges <- rbind(d1, d2)
#
# # create a datafram with connection between genes
# all_leaves <- data_unique$combined_id
# comb_genes <- inner_join(x = data_unique, y = data_unique, by = "ensembl_id")
# connect <- data.frame(from = comb_genes$combined_id.x, to = comb_genes$combined_id.y)
# connect$value <- 1
#
# # create a vertices data-frame
# vertices <- data.frame(
#   name = unique(c(as.character(edges$from), as.character(edges$to))) ,
#   value = 1
# )
# vertices$group  <-  edges$from[ match( vertices$name, edges$to ) ]
#
#
# # Create a graph object
# mygraph <- igraph::graph_from_data_frame( edges, vertices=vertices )
#
# # The connection object must refer to the ids of the leaves:
# from  <-  match( connect$from, vertices$name)
# to  <-  match( connect$to, vertices$name)
#
# # Basic usual argument
# ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
#   geom_conn_bundle(data = get_con(from = from, to = to), alpha=0.2, colour="skyblue", tension = .5) +
#   geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group)) +
#   scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
#   theme_void()
#
# # ggsave(filename = paste0(folder_plots, "/06_comparison_deseq_circular.png"))
