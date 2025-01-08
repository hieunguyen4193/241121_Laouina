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
##### RUN DESEQ2
#####----------------------------------------------------------------------#####
run_DESeq2_and_preprocess <- function(  deseq.dataset, 
                                        tx2gene, 
                                        thresh.pval = 0.05,
                                        geneid.to.remove = NULL,
                                        explicit.size.factor = TRUE,
                                        hk.gene = "Gapdh"){
  # run DESeq2
  if (is.null(geneid.to.remove) == FALSE){
    deseq.dataset <- deseq.dataset[rownames(deseq.dataset) %in% geneid.to.remove == FALSE]
  }
  print(length(unique(colData(deseq.dataset)$organ)))
  
  if (length(unique(colData(deseq.dataset)$organ)) < 2){
    print(sprintf("Only %s samples exist in this comparison, do not perform size factor estimation", unique(colData(deseq.dataset)$organ)))
    explicit.size.factor <- FALSE
    print(sprintf("Setting explicit size factor to %s", explicit.size.factor))
  } 
  
  if (explicit.size.factor == FALSE){
    deseq.object <- DESeq(deseq.dataset)    
  } else {
    print("Running DESEQ2 with input SIZE FACTOR...")
    meta.data <- colData(deseq.dataset)
    sample.list <- list()
    for (o in unique(meta.data$organ)){
      sample.list[[o]] <- subset(meta.data, meta.data$organ == o)$sample
    }
    all.hk.genes <- c("Rpl13a", "Tbp", "Hprt1", "Actb", "Gapdh")
    
    hkdf <- counts(deseq.dataset) %>%  as.data.frame() %>% rownames_to_column("gene_id")
    genedf <- subset(tx2gene, select = c(gene_id, gene_name)) 
    genedf <- genedf[!duplicated(genedf$gene_id), ]
    
    hkdf <- merge(hkdf, genedf, by.x = "gene_id", by.y = "gene_id")
    hkdf <- subset(hkdf, select = -c(gene_id)) %>% subset(gene_name %in% all.hk.genes)
    row.names(hkdf) <- NULL
    hkdf <- hkdf %>%
      column_to_rownames("gene_name")
    
    # Fold change = distal / proximal
    mean.hkdf.1 <- hkdf[hk.gene, sample.list$distal] %>% as.numeric() %>% mean()
    mean.hkdf.2 <- hkdf[hk.gene, sample.list$proximal] %>% as.numeric() %>% mean()
    fc <- mean.hkdf.1/mean.hkdf.2
    
    all.samples <- colData(deseq.dataset)$organ
    sizeFactors(deseq.dataset) <- to_vec(
      for (item in all.samples){
        if (item == "proximal"){
          fc
        } else {
          1
        }
      }
    )
    deseq.dataset <- estimateDispersions(deseq.dataset)
    deseq.object <- nbinomWaldTest(deseq.dataset)
  }
  
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

#####----------------------------------------------------------------------#####
##### PRELIMINARY PLOTS
#####----------------------------------------------------------------------#####

preliminary_plot_from_deseq_output <- function(countdf, outputdir){
  dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)
  #####----------------------------------------------------------------#####
  ##### HEATMAP OF SELECTED GENES
  #####----------------------------------------------------------------#####
  all.gene.list <- list(list1 = c("Lct", 
                                  "Lyz1", 
                                  "Tff3", 
                                  "Zg16", 
                                  "Chga", 
                                  "Olfm4", 
                                  "Dclk1"),
                        list2 = c("Pigr", 
                                  "Ccl28", 
                                  "Ccl25", 
                                  "Aprill", 
                                  "Baff", 
                                  "Tgf-beta", 
                                  "Il-6", 
                                  "Il-7", 
                                  "Cxcl9",
                                  "Cxcl10", 
                                  "Il-22")
  )
  
  for (i in seq(length(all.gene.list))){
    gene.list <- all.gene.list[[i]]
    countdf.filtered <- subset(countdf, gene_name %in% gene.list) %>%
      subset(select = -c(gene_id)) 
    row.names(countdf.filtered) <- NULL
    
    ##### heatmap norm by row
    heatmapdf <- countdf.filtered %>% column_to_rownames("gene_name") %>% t() %>% as.data.frame()
    
    heatmapdf.scaled <- (heatmapdf - rowMeans(heatmapdf))/rowSds(as.matrix(heatmapdf))
    heatmap.plot <- heatmapdf.scaled %>% rownames_to_column("SampleID") %>% 
      pivot_longer(!SampleID, names_to = "signature", values_to = "z_score") %>%
      ggplot(aes(y = SampleID, x = signature, fill = z_score)) + geom_tile() + 
      scale_fill_distiller(palette = "RdBu") + 
      theme(axis.text = element_text(size = 22))
    ggsave(plot = heatmap.plot, filename = sprintf("heatmap_scale_within_gene.list%s.svg", i), path = outputdir, dpi = 300, width = 10, height = 10)
    
    ##### heatmap norm by column
    heatmapdf <- countdf.filtered %>% column_to_rownames("gene_name") %>% as.data.frame()
    
    heatmapdf.scaled <- (heatmapdf - rowMeans(heatmapdf))/rowSds(as.matrix(heatmapdf))
    heatmap.plot <- heatmapdf.scaled %>% rownames_to_column("SampleID") %>% 
      pivot_longer(!SampleID, names_to = "signature", values_to = "z_score") %>%
      ggplot(aes(y = SampleID, x = signature, fill = z_score)) + geom_tile() + 
      scale_fill_distiller(palette = "RdBu") + 
      theme(axis.text.x = element_text(size = 22, angle = 90),
            axis.text.y = element_text(size = 22))
    ggsave(plot = heatmap.plot, filename = sprintf("heatmap_scale_within_sample.list%s.svg", i), path = outputdir, dpi = 300, width = 10, height = 10)
  }
  
  #####----------------------------------------------------------------#####
  ##### PCA plot
  #####----------------------------------------------------------------#####
  input.df <- countdf[meta.data$sample]
  pca.object <- prcomp(t(input.df), rank. = 2, scale. = FALSE)
  
  pcadf <- data.frame(pca.object$x)
  
  row.names(pcadf) <- meta.data$sample
  
  pcadf <- merge(pcadf, meta.data, by.x = "row.names", by.y = "sample")
  
  for (c in c("organ", "age", "condition", "full.condition")){
    pca.plot <- ggplot(data.frame(pcadf), aes_string(x="PC1", y="PC2", color= c, text = "Row.names", 
                                                     label = "SampleName")) +
      geom_point(size = 4) +
      theme_bw() + 
      theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
      geom_label_repel() + 
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0)
    ggsave(plot = pca.plot, filename = sprintf("PCA.%s.svg", c), path = outputdir, dpi = 300, width = 10, height = 10)
  }
  

  
  #####----------------------------------------------------------------------#####
  ##### Update 04.01.2024
  #####----------------------------------------------------------------------#####
  # Check gene expression of house keeping genes. 
  
  all.hk.genes <- c("Rpl13a", "Tbp", "Hprt1", "Actb", "Gapdh")
  hkdf <- countdf %>% subset(select = -c(gene_id)) %>%
    subset(gene_name %in% all.hk.genes) 
  
  row.names(hkdf) <- NULL
  hkdf <- hkdf %>%
    column_to_rownames("gene_name")
  hkdf.raw <- hkdf
  hkdf <- (hkdf - rowMeans(hkdf))/rowSds(as.matrix(hkdf))
  
  hkdf.pivot <- hkdf %>% rownames_to_column("gene_name") %>% pivot_longer(!gene_name, names_to = "SampleID", values_to = "exprs")
  
  hkdf.pivot <- merge(hkdf.pivot, meta.data, by.x = "SampleID", by.y = "sample")
  
  compare.plots <- list()
  for (c in c("organ", "age", "condition")){
    p.compare.all.hk.genes <- data.frame(hkdf.pivot) %>% ggplot(aes_string(x = c, y = "exprs", fill = c)) + geom_boxplot() + 
      facet_wrap(~gene_name) 
    compare.plots[[c]] <- p.compare.all.hk.genes
    ggsave(plot = p.compare.all.hk.genes, 
           filename = sprintf("compare_HK_genes.%s.svg", c), 
           path = outputdir,
           dpi = 300, 
           width = 10, 
           height = 10)
  }
}

#####----------------------------------------------------------------------#####
##### plot from deseq output
#####----------------------------------------------------------------------#####
plot_result_from_deseq_output <- function(input.deseq.output, path.to.save.output, input.metadata){
  sigdf <- input.deseq.output$resultdf.sig
  nonsigdf <- input.deseq.output$resultdf.nonsig
  
  writexl::write_xlsx(sigdf, file.path(path.to.save.output, "sig_DGE.xlsx"))
  writexl::write_xlsx(nonsigdf, file.path(path.to.save.output, "nonsig_DGE.xlsx"))
  
  ##### PCA plot
  input.df <- input.deseq.output$norm.count[input.metadata$sample]
  pca.object <- prcomp(t(input.df), rank. = 2, scale. = FALSE)
  
  pcadf <- data.frame(pca.object$x)
  
  row.names(pcadf) <- input.metadata$sample
  
  pcadf <- merge(pcadf, input.metadata, by.x = "row.names", by.y = "sample")
  
  pca.plot <- ggplot(pcadf, aes(x=PC1, y=PC2, color= full.condition, text = Row.names, label = SampleName)) +
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
  
  input.df <- input.deseq.output$all.resultdf
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
                    levels(input.metadata$full.condition)[[2]],
                    levels(input.metadata$full.condition)[[1]]
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
  
  input.to.heatmap <- subset(sig.genes.with.highlogFC, select = c("gene_id", "gene_name", input.metadata$sample))
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
  saveRDS(input.deseq.output, file.path(path.to.save.output, "deseq_output.rds"))
  write.csv(data.frame(status = c("finished")), 
            file.path(path.to.save.output, "finished.csv"))
}
