library(ggforce)
library(reshape2)

### so : seurat object
### group.by : name of metadata column indicate experiment condition
### control : name of control group
### genes : genes of interest to get proportion and average expression
get_exp_proportion_vec <- function(so, group.by, control, genes) {
    control.cells = Cells(so)[so@meta.data[, group.by] == control]
    exp.cells = Cells(so)[so@meta.data[, group.by] != control]
    
    exp = c()
    proportion = c()
    for (gene in genes) {
        control.exp = mean(so@assays$RNA$data[gene, control.cells] %>% expm1)
        exp.exp = mean(so@assays$RNA$data[gene, exp.cells] %>% expm1)
        exp = c(exp, c(control.exp, exp.exp))

        control.proportion = sum(so@assays$RNA$data[gene, control.cells] > 0) / length(control.cells)
        exp.proportion = sum(so@assays$RNA$data[gene, exp.cells] > 0) / length(exp.cells)
        proportion = c(proportion, c(control.proportion, exp.proportion))
    }
    return(list(average_expression = exp, proportion = proportion))
}

####### load p_value_df
### row : gene
### column : celltype
### element : p-value of certain test
# write.csv(p_value_df, file = paste0(save_path, 'p_value_df.csv'), row.names = T)
p_value_df <- read.csv(paste0(save_path, 'p_value_df.csv'), row.names = 1)
pval_melt <- melt(p_value_df, id.vars = "gene", variable.name = "celltype", value.name = "p_value")

# assign *
pval_melt$signif <- cut(pval_melt$p_value,
                        breaks = c(0, 0.001, 0.01, 0.05, 1),
                        labels = c("***", "**", "*", ""),
                        include.lowest = TRUE)
pval_melt$x0 <- 0 
pval_melt$y0 <- match(pval_melt$celltype, Tcell_cols_order) + 0.3

### prepare proper data format for mushroom plot
so.subset <- list()
df.list <- list()
for (celltype in Tcell_cols_order) {
    ### so.cd4 : seurat object (CD4 T cell only)
    so.subset[[celltype]] = subset(so.cd4, cells = Cells(so.cd4)[so.cd4$dpi == '14dpi' & so.cd4$annotation %in% c(celltype)])
    
    exp_proportion = get_exp_proportion_vec(so.subset[[celltype]], 'condition', 'floxed', gene_of_interest)
    exp = exp_proportion$average_expression; proportion = exp_proportion$proportion
    ## dAT2 : name of experiment group
    df <- data.frame(condition = rep(c('floxed', 'dAT2'), length(gene_of_interest)),
                     gene = rep(gene_of_interest, each = 2), 
                     proportion = proportion,
                     avg.expression = exp,
                     celltype = rep(celltype, length(gene_of_interest) * 2))
    df.list[[celltype]] <- df
}
df <- do.call(rbind, df.list)

df$proportion <- ifelse(df$proportion != 0, df$proportion + min.size, df$proportion)
# write.csv(df, paste0(save_path, 'Areg.proportion.csv'), row.names = T)
df <- df %>%
  mutate(
    y0 = as.numeric(factor(celltype, levels = Tcell_cols_order)),
    start = ifelse(condition == "floxed", pi, 0),
    end   = ifelse(condition == "floxed", 2* pi, pi),
    r = proportion/scale_factor # scale factor for proportion visualization
  )

p <- ggplot(df, aes(x0 = 0, y0 = y0)) +
  geom_arc_bar(aes(start = start, end = end, r0 = 0, r = r, fill = avg.expression),
               color = "black", size = 0.1) +
  geom_point(aes(size = proportion), data = df[0,]) +
  scale_size_continuous(name = "Proportion") +
  scale_fill_gradient2(name = "Avg Expression", low = "white", mid = '#82cbff', high = "#060048", midpoint = (df$avg.expression %>% max + df$avg.expression %>% min)/2) +
  facet_wrap(~ gene, nrow = 1, strip.position = "bottom") +
  labs(x = "Gene", y = "Celltype") +
  scale_y_continuous(breaks = seq(1, as.numeric(factor(celltype, levels = Tcell_cols_order))),
                     labels = Tcell_cols_order) +
  coord_fixed() + 
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(color = "black"),
    legend.position = "right",
    axis.title = element_text(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )
# add significance *
p <- p + 
  geom_text(data = pval_melt,
            aes(x = x0, y = y0, label = signif),
            inherit.aes = FALSE,
            color = "black",
            size = 3)
            
ggsave(p, file = paste0(save_path, 'test_mushroom.png'), width = 5, height = 5)
ggsave(p, file = paste0(save_path, 'test_mushroom.pdf'), width = 5, height = 5)
ggplot2pptx(p, 5, 5, paste0(save_path, 'test_mushroom.pptx')
