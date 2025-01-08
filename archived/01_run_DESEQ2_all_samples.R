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
                    show_col_types = FALSE)

deseq.output <- run_DESeq2_and_preprocess(deseq.dataset, tx2gene, thresh.pval = 0.05, explicit.size.factor = FALSE)

#####----------------------------------------------------------------------#####
##### get Fold change of HK genes w.r.t ORGAN
#####----------------------------------------------------------------------#####
sample.list <- list()
for (o in unique(meta.data$organ)){
  sample.list[[o]] <- subset(meta.data, meta.data$organ == o)$sample
}
hk.genes <- c("Rpl13a", "Tbp", "Hprt1", "Actb", "Gapdh")
hkdf <- counts(deseq.dataset) %>%  as.data.frame() %>% rownames_to_column("gene_id")

genedf <- subset(tx2gene, select = c(gene_id, gene_name)) 
genedf <- genedf[!duplicated(genedf$gene_id), ]

hkdf <- merge(hkdf, genedf, by.x = "gene_id", by.y = "gene_id")
hkdf <- subset(hkdf, select = -c(gene_id)) %>% subset(gene_name %in% hk.genes)
row.names(hkdf) <- NULL
hkdf <- hkdf %>%
  column_to_rownames("gene_name")
hkdf.raw <- hkdf

# Fold change = distal / proximal
mean.hkdf.1 <- hkdf.raw["Gapdh", sample.list$distal] %>% as.numeric() %>% mean()
mean.hkdf.2 <- hkdf.raw["Gapdh", sample.list$proximal] %>% as.numeric() %>% mean()
fc <- mean.hkdf.1/mean.hkdf.2

preliminary_plot_from_deseq_output(input.deseq.output = deseq.output, 
                                   outputdir = file.path(path.to.01.output, "raw_deseq_output"))

###### see https://support.bioconductor.org/p/130564/

#####
explicit.size.factor <- TRUE
hk.gene <- "Gapdh"
deseq.output.corrected <- run_pipeline_DESEQ2(salmon.output = salmon.output,
                                    path.to.save.output = file.path(path.to.01.output, "corrected"),
                                    filtered.metadata = meta.data,
                                    thresh.pval = 0.05,
                                    explicit.size.factor = explicit.size.factor,
                                    hk.gene = hk.gene)
preliminary_plot_from_deseq_output(input.deseq.output = deseq.output.corrected, 
                                   outputdir = file.path(path.to.01.output, "corrected"))