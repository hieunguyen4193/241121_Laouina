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
source(file.path(path.to.main.src, "helper_functions.R"))

meta.data <- readxl::read_excel(file.path(path.to.main.src, "Sample sheet Hieu.xlsx")) %>%
  rowwise() %>%
  mutate(SampleID = str_split(SampleName, "_")[[1]][[1]] ) %>%
  mutate(condition = str_split(SampleName, "_")[[1]][[2]]) %>%
  mutate(age = str_split(SampleName, "_")[[1]][[3]]) %>%
  mutate(organ = ifelse(str_split(SampleName, "")[[1]][[1]] == "D", "distal", "proximal")) %>%
  mutate(sample = str_split(FileName, "_R1")[[1]][[1]])%>%
  mutate(full.condition = sprintf("%s_%s_%s", condition, age, organ))

#####----------------------------------------------------------------#####
##### Generate raw and normalized count data without HOUSE KEEPING GENES
#####----------------------------------------------------------------#####
salmon.output <- "/media/hieunguyen/HD01/outdir/CRC1382/241121_Laouina_Pabst_MolMedizin_3mRNAseq/star_salmon"
path.to.tx2gene <- file.path(salmon.output, "tx2gene.tsv")
deseq.dataset <- generate_DESeq2_dataset(path.to.salmon.quants = salmon.output,
                                         path.to.tx2gene = path.to.tx2gene,
                                         metadata = meta.data)
tx2gene <- read_tsv(path.to.tx2gene, 
                    show_col_types = FALSE)
tmp.tx2gene <- tx2gene[, c("gene_id", "gene_name")]
tmp.tx2gene <- tmp.tx2gene[!duplicated(tmp.tx2gene$gene_name),]

deseq.dataset <- generate_DESeq2_dataset(path.to.salmon.quants = salmon.output,
                                         path.to.tx2gene = path.to.tx2gene,
                                         metadata = meta.data)

deseq.dataset <- estimateSizeFactors(deseq.dataset)
deseq.dataset <- estimateDispersions(deseq.dataset)

norm.count <- counts(deseq.dataset, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("gene_id")
norm.count <- merge(norm.count, tmp.tx2gene, by.x = "gene_id", by.y = "gene_id")

raw.count <- counts(deseq.dataset, normalized = FALSE) %>% as.data.frame() %>% rownames_to_column("gene_id")
raw.count <- merge(raw.count, tmp.tx2gene, by.x = "gene_id", by.y = "gene_id")

preliminary_plot_from_deseq_output(norm.count, file.path(path.to.01.output, "norm_count"))
preliminary_plot_from_deseq_output(raw.count, file.path(path.to.01.output, "raw_count"))

##### estimate TRUE size factor:
meta.data <- colData(deseq.dataset)
sample.list <- list()
for (o in unique(meta.data$organ)){
  sample.list[[o]] <- subset(meta.data, meta.data$organ == o)$sample
}

all.hk.genes <- c("Rpl13a", "Tbp", "Hprt1", "Actb", "Gapdh")
hk.gene <- "Gapdh"

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

##### manually assign size factor for proximal and distal sample

sizeFactors(deseq.dataset) <- to_vec(
  for (item in colData(deseq.dataset)$organ){
    if (item == "proximal"){
      fc
    } else {
      1
    }
  }
)
deseq.dataset.corrected <- estimateDispersions(deseq.dataset)

norm.count.corrected <- counts(deseq.dataset.corrected, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("gene_id")
norm.count.corrected <- merge(norm.count.corrected, tmp.tx2gene, by.x = "gene_id", by.y = "gene_id")
preliminary_plot_from_deseq_output(norm.count.corrected, file.path(path.to.01.output, "norm_count_corrected_Gapdh"))
