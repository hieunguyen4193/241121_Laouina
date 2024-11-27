gc()
rm(list = ls())

path.to.main.src <- "/home/hieunguyen/CRC1382/src_2023/241121_Laouina_Pabst_MolMedizin_3mRNAseq"

if ("svglite" %in% installed.packages() == FALSE){
  install.packages("svglite")
}

outdir <- "/home/hieunguyen/CRC1382/outdir"
PROJECT <- "ALaouina_241121"
  
path.to.main.output <- file.path(outdir, PROJECT)
dir.create(path.to.main.output, showWarnings = FALSE, recursive = TRUE)
path.to.01.output <- file.path(path.to.main.output, "01_output")
dir.create(path.to.01.output, showWarnings = FALSE, recursive = TRUE)

source(file.path(path.to.main.src, "preparation.R"))
source(file.path(path.to.main.src, "helper_functions_v2.R"))

meta.data <- readxl::read_excel(file.path(path.to.main.src, "Sample sheet Hieu.xlsx")) %>%
  rowwise() %>%
  mutate(SampleID = str_split(SampleName, "_")[[1]][[1]] ) %>%
  mutate(condition = str_split(SampleName, "_")[[1]][[2]]) %>%
  mutate(age = str_split(SampleName, "_")[[1]][[3]]) %>%
  mutate(organ = ifelse(str_split(SampleName, "")[[1]][[1]] == "D", "distal", "proximal")) %>%
  mutate(sample = str_split(FileName, "_R1")[[1]][[1]])%>%
  mutate(full.condition = sprintf("%s_%s_%s", condition, age, organ))

#####----------------------------------------------------------------#####
##### RUN DESEQ2
#####----------------------------------------------------------------#####
salmon.output <- "/media/hieunguyen/HD01/outdir/CRC1382/241121_Laouina_Pabst_MolMedizin_3mRNAseq/star_salmon"
path.to.tx2gene <- file.path(salmon.output, "tx2gene.tsv")
deseq.dataset <- generate_DESeq2_dataset(path.to.salmon.quants = salmon.output,
                                         path.to.tx2gene = path.to.tx2gene,
                                         metadata = meta.data)
tx2gene <- read_tsv(path.to.tx2gene, 
                    col_names = c("transcript_id", "gene_id", "gene_name"), 
                    show_col_types = FALSE)
deseq.output <- run_DESeq2_and_preprocess(deseq.dataset, tx2gene, thresh.pval = 0.05)

#####----------------------------------------------------------------#####
##### HEATMAP OF SELECTED GENES
#####----------------------------------------------------------------#####
countdf <- deseq.output$norm.count 
gene.list <- c("Lct", "Lyz1", "Tff3", "Zg16", "Chga", "Olfm4", "Dclk1")

countdf.filtered <- subset(countdf, gene_name %in% gene.list) %>%
  subset(select = -c(gene_id)) 
row.names(countdf.filtered) <- NULL
heatmapdf <- countdf.filtered %>% column_to_rownames("gene_name") %>% t() %>% as.data.frame()

heatmapdf.scaled <- (heatmapdf - rowMeans(heatmapdf))/rowSds(as.matrix(heatmapdf))
heatmap.plot <- heatmapdf.scaled %>% rownames_to_column("SampleID") %>% 
  pivot_longer(!SampleID, names_to = "signature", values_to = "z_score") %>%
  ggplot(aes(y = SampleID, x = signature, fill = z_score)) + geom_tile() + 
  scale_fill_distiller(palette = "RdBu") + 
  theme(axis.text = element_text(size = 22))
ggsave(plot = heatmap.plot, filename = "heatmap.svg", path = path.to.01.output, dpi = 300, width = 10, height = 10)

#####----------------------------------------------------------------#####
##### PCA plot
#####----------------------------------------------------------------#####
input.df <- deseq.output$norm.count[meta.data$sample]
pca.object <- prcomp(t(input.df), rank. = 2, scale. = FALSE)

pcadf <- data.frame(pca.object$x)

row.names(pcadf) <- meta.data$sample

pcadf <- merge(pcadf, meta.data, by.x = "row.names", by.y = "sample")

pca.plot <- ggplot(pcadf, aes(x=PC1, y=PC2, color= condition, text = Row.names, label = Row.names)) +
  geom_point(size = 4) +
  theme_bw() + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 12)) +
  geom_label_repel() + 
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)
ggsave(plot = pca.plot, filename = "PCA.svg", path = path.to.01.output, dpi = 300, width = 10, height = 10)
