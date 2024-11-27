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
source(file.path(path.to.main.src, "preparation.R"))
source(file.path(path.to.main.src, "helper_functions_v2.R"))

meta.data <- readxl::read_excel(file.path(path.to.main.src, "Sample sheet Hieu.xlsx")) %>%
  rowwise() %>%
  mutate(SampleID = str_split(SampleName, "_")[[1]][[1]] ) %>%
  mutate(condition = str_split(SampleName, "_")[[1]][[2]]) %>%
  mutate(age = str_split(SampleName, "_")[[1]][[3]]) %>%
  mutate(organ = ifelse(str_split(SampleName, "")[[1]][[1]] == "D", "distal", "proximal")) %>%
  mutate(sample = str_split(FileName, "_R1")[[1]][[1]])

salmon.output <- "/media/hieunguyen/HD01/outdir/CRC1382/241121_Laouina_Pabst_MolMedizin_3mRNAseq/star_salmon"

path.to.tx2gene <- file.path(salmon.output, "tx2gene.tsv")
deseq.dataset <- generate_DESeq2_dataset(path.to.salmon.quants = salmon.output,
                                         path.to.tx2gene = path.to.tx2gene,
                                         metadata = meta.data)
tx2gene <- read_tsv(path.to.tx2gene, 
                    col_names = c("transcript_id", "gene_id", "gene_name"), 
                    show_col_types = FALSE)
deseq.output <- run_DESeq2_and_preprocess(deseq.dataset, tx2gene, thresh.pval = 0.05)