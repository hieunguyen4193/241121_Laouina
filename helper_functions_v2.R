#####----------------------------------------------------------------------#####
##### generate DESeq2 dataset object
#####----------------------------------------------------------------------#####
generate_DESeq2_dataset <- function(path.to.salmon.quants, 
                                    path.to.tx2gene,
                                    metadata, 
                                    file.ext = "quant.sf",
                                    ERCC.spike.in = FALSE){
  # fetch all quant.sf files
  files <- file.path(path.to.salmon.quants, metadata$sample, file.ext)
  names(files) <- metadata$sample
  
  # read in the tx2gene list
  tx2gene <- read_tsv(path.to.tx2gene, 
                      col_names = c("transcript_id", "gene_id", "gene_name"), 
                      show_col_types = FALSE)
  txi <- tximport(files, type="salmon", tx2gene=tx2gene)
  output <- DESeqDataSetFromTximport( txi,
                                      colData = metadata,
                                      design = ~ full.condition) # make sure that group exists in metadata
  return(output)
}

#####----------------------------------------------------------------------#####
##### Generate interactive table in html file
#####----------------------------------------------------------------------#####
create_dt <- function(x){
  DT::datatable(x,
                extensions = 'Buttons',
                filter = "top",
                options = list(dom = 'Blfrtip',
                               buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                               lengthMenu = list(c(10,25,50,-1),
                                                 c(10,25,50,"All")),
                               columnDefs = list(list(
                                 targets = "_all",
                                 render = JS(
                                   "function(data, type, row, meta) {",
                                   "return type === 'display' && data != null && data.length > 100 ?",
                                   "'<span title=\"' + data + '\">' + data.substr(0, 100) + '...</span>' : data;",
                                   "}")
                               ))
                ))
}

#####----------------------------------------------------------------------#####
##### RUN DESEQ2
#####----------------------------------------------------------------------#####
run_DESeq2_and_preprocess <- function(  deseq.dataset, 
                                        tx2gene, 
                                        thresh.pval = 0.05,
                                        geneid.to.remove = NULL){
  
  # run DESeq2
  if (is.null(geneid.to.remove) == FALSE){
    deseq.dataset <- deseq.dataset[rownames(deseq.dataset) %in% geneid.to.remove == FALSE]
  }
  
  deseq.object <- DESeq(deseq.dataset)
  
  # normalized read counts
  norm.counts <- counts(deseq.object, normalized=TRUE)
  
  # preprocessing
  ## extract the Differential expression analysis results from DESeq2 object
  ## this object contains: "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  deseq.result <- results(deseq.object)
  
  # generate a dataframe containing gene.id and gene.name
  genedf <- data.frame(gene_id = tx2gene$gene_id, 
                       gene_name = tx2gene$gene_name) 
  
  ## remove duplicated gene.name / gene.id
  genedf <- genedf[!duplicated(genedf), ]
  
  ## merge the normalized count matrix and the gene information altogether
  norm.counts <- merge(genedf, norm.counts, by.x= "gene_id", by.y= "row.names")
  
  ## merge the normalized count matrix (with gene information) with DESeq2 results
  resultdf <- merge(genedf, as.data.frame(deseq.result), by.x= "gene_id", by.y="row.names", all.x=F, all.y=T)
  
  resultdf <- merge(resultdf, norm.counts, by=c("gene_id", "gene_name"), all.x=T, all.y=F)
  
  ## remove rows with NA values
  resultdf <- resultdf[complete.cases(resultdf), ]
  
  ## marking significant and non-significant genes by the DESeq2 p-value
  resultdf$sig <- "Non-sig."
  
  resultdf$sig[resultdf$padj < thresh.pval] <- "Sig. genes"
  
  ## Marking ERCC genes
  sel_ERCC <- str_detect(resultdf$gene_id, "^ERCC-*gene")
  
  resultdf$sig[sel_ERCC] <- "Spike in"
  
  ## remove ERCC spike-in genes
  resultdf <- resultdf[!sel_ERCC,]
  
  ## extract the dataframe with significant genes only
  resultdf.sig <- resultdf[resultdf$padj < thresh.pval,]
  
  ## extract the dataframe with significant genes only
  resultdf.nonsig <- resultdf[resultdf$padj >= thresh.pval,]
  
  # collect all outputs
  output <- list(norm.count = norm.counts,
                 all.resultdf = resultdf,
                 resultdf.sig = resultdf.sig,
                 resultdf.nonsig = resultdf.nonsig,
                 deseq.object = deseq.object)
  return(output)
}


##### MAIN FUNCTION TO RUN DESEQ2 AND PLOTS

run_pipeline_DESEQ2 <- function(salmon.output, 
                                path.to.save.output,
                                filtered.metadata){
  ##### DESEQ2
  path.to.tx2gene <- file.path(salmon.output, "tx2gene.tsv")
  deseq.dataset <- generate_DESeq2_dataset(path.to.salmon.quants = salmon.output,
                                           path.to.tx2gene = path.to.tx2gene,
                                           metadata = filtered.metadata)
  tx2gene <- read_tsv(path.to.tx2gene, 
                      col_names = c("transcript_id", "gene_id", "gene_name"), 
                      show_col_types = FALSE)
  deseq.output <- run_DESeq2_and_preprocess(deseq.dataset, tx2gene, thresh.pval = 0.05)
  
  sigdf <- deseq.output$resultdf.sig
  nonsigdf <- deseq.output$resultdf.nonsig
  
  writexl::write_xlsx(sigdf, file.path(path.to.save.output, "sig_DGE.xlsx"))
  writexl::write_xlsx(nonsigdf, file.path(path.to.save.output, "nonsig_DGE.xlsx"))
  
  ##### PCA plot
  input.df <- deseq.output$norm.count[filtered.metadata$sample]
  pca.object <- prcomp(t(input.df), rank. = 2, scale. = FALSE)
  
  pcadf <- data.frame(pca.object$x)
  
  row.names(pcadf) <- filtered.metadata$sample
  
  pcadf <- merge(pcadf, filtered.metadata, by.x = "row.names", by.y = "sample")
  
  pca.plot <- ggplot(pcadf, aes(x=PC1, y=PC2, color= condition, text = Row.names, label = Row.names)) +
    geom_point(size = 4) +
    ggtitle(sprintf("PCA plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    geom_label_repel() + 
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0)
  
  ggsave(plot = pca.plot, filename = "PCA.svg", path = path.to.save.output, dpi = 300, width = 10, height = 10)
  
  ##### Volcano plot
  cutoff.adjp <- 0.05
  
  input.df <- deseq.output$all.resultdf
  input.df <- input.df %>% mutate(abs.log2FoldChange = abs(log2FoldChange))
  input.df <- input.df %>% rowwise() %>%
    mutate(show.gene.name = ifelse(padj < cutoff.adjp, gene_name, NA))
  
  volcano.plot <- ggplot(data=input.df, 
                         aes(x=log2FoldChange, y=-log10(padj), col=sig, label=show.gene.name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    geom_vline(xintercept=c(-1, 1), col="#9a9fa6", linetype='dotted') +
    geom_hline(yintercept=-log10(cutoff.adjp), col="#9a9fa6", linetype='dotted') +
    geom_text_repel() +
    ggtitle(sprintf("Volcano plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
    xlim(c(-max(input.df$abs.log2FoldChange), max(input.df$abs.log2FoldChange))) +
    ggtitle(sprintf("Positive LogFC: %s > %s (right-handed side)", 
                        levels(filtered.metadata$full.condition)[[2]],
                        levels(filtered.metadata$full.condition)[[1]]
                        ))
  ggsave(plot = volcano.plot, filename = "volcano_plot.svg", path = path.to.save.output, dpi = 300, width = 14, height = 10)
  
  ##### MA plot
  ma.plot <- ggplot(data=input.df, 
                    aes(x=log2(baseMean), y=log2FoldChange, col=sig, label=show.gene.name)) + 
    geom_point() + 
    scale_color_manual(values=c("#c0d2f0", "#f28095")) +
    theme_minimal() +
    geom_text_repel() +
    ggtitle(sprintf("MA plot, Sample: %s vs. %s", condition1, condition2)) +
    theme_bw() + 
    theme(plot.title = element_text(hjust=0.5, face="bold", size = 12))
  ggsave(plot = ma.plot, filename = "MA_plot.svg", path = path.to.save.output, dpi = 300, width = 14, height = 10)
  
  ##### heatmap
  sig.genes.with.highlogFC.up <- subset(input.df, (input.df$sig == "Sig. genes") & (input.df$log2FoldChange > 0)) %>% arrange(desc(log2FoldChange)) %>% head(20)
  sig.genes.with.highlogFC.down <- subset(input.df, (input.df$sig == "Sig. genes") & (input.df$log2FoldChange < 0)) %>% arrange(desc(log2FoldChange)) %>% tail(20)
  sig.genes.with.highlogFC <- rbind(sig.genes.with.highlogFC.up, sig.genes.with.highlogFC.down)
  
  input.to.heatmap <- subset(sig.genes.with.highlogFC, select = c("gene_id", "gene_name", filtered.metadata$sample))
  if (nrow(input.to.heatmap) > 0){
    heatmap.values <- log10(input.to.heatmap[,3:(dim(input.to.heatmap)[2])] + 1)
    selected.genes.heatmap.plot <- heatmaply(heatmap.values, 
                                             main= sprintf("Heatmap, Sample: %s vs. %s", condition1, condition2),
                                             method = "plotly",labRow=input.to.heatmap$gene_name,
                                             xlab = "Samples", ylab = "Genes", width = 1000, height = 600,
                                             showticklabels = c(TRUE, FALSE), show_dendrogram = c(FALSE, TRUE),
                                             key.title = "log10 scale colormap",
                                             label_names = c("Gene", "Sample", "Expression"),
                                             k_col = 2, file = file.path(path.to.save.output, "heatmap.html"))
  }
  saveRDS(deseq.output, file.path(path.to.save.output, "deseq_output.rds"))
  return(deseq.output)
}