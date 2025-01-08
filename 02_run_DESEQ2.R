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
path.to.02.output <- file.path(path.to.main.output, "02_output")
dir.create(path.to.02.output, showWarnings = FALSE, recursive = TRUE)

source(file.path(path.to.main.src, "preparation.R"))
source(file.path(path.to.main.src, "helper_functions.R"))

meta.data <- readxl::read_excel(file.path(path.to.main.src, "Sample sheet Hieu.xlsx")) %>%
  rowwise() %>%
  mutate(SampleID = str_split(SampleName, "_")[[1]][[1]] ) %>%
  mutate(condition = str_split(SampleName, "_")[[1]][[2]]) %>%
  mutate(age = str_split(SampleName, "_")[[1]][[3]]) %>%
  mutate(organ = ifelse(str_split(SampleName, "")[[1]][[1]] == "D", "distal", "proximal")) %>%
  mutate(sample = str_split(FileName, "_R1")[[1]][[1]]) %>%
  mutate(full.condition = sprintf("%s_%s_%s", condition, age, organ))

all.comp <- list(
  set1 = data.frame(
    condition1 = c("SPF_young_proximal", "SPF_young_distal", "SPF_old_proximal", "SPF_old_distal"),
    condition2 = c("GF_young_proximal", "GF_young_distal", "GF_old_proximal", "GF_old_distal")
  ),
  set2 = data.frame(
    condition1 = c("SPF_young_proximal", "SPF_old_proximal", "GF_young_proximal", "GF_old_proximal"),
    condition2 = c("SPF_young_distal", "SPF_old_distal", "GF_young_distal", "GF_old_distal")
  ),
  set3 = data.frame(
    condition1 = c("SPF_young_proximal", "SPF_young_distal", "GF_young_proximal", "GF_young_distal"),
    condition2 = c("SPF_old_proximal", "SPF_old_distal", "GF_old_proximal", "GF_old_distal")
  )
)


salmon.output <- "/media/hieunguyen/HD01/outdir/CRC1382/241121_Laouina_Pabst_MolMedizin_3mRNAseq/star_salmon"
path.to.tx2gene <- file.path(salmon.output, "tx2gene.tsv")
tx2gene <- read_tsv(path.to.tx2gene, show_col_types = FALSE)

hk.gene <- "Gapdh"

for (explicit.size.factor in c(TRUE, FALSE)){
  for (input.set in names(all.comp)){
    input.samples <- all.comp[[input.set]]
    for (i in seq(1, nrow(input.samples))){
      condition1 <- input.samples[i, "condition1"]
      condition2 <- input.samples[i, "condition2"]
      input.metadata <- subset(meta.data, meta.data$full.condition %in% c(condition1, condition2))
      
      ##### FC = condition2/condition1
      input.metadata$full.condition <- factor(input.metadata$full.condition, levels = c(condition1, condition2))
      
      deseq.dataset <- generate_DESeq2_dataset(path.to.salmon.quants = salmon.output,
                                               path.to.tx2gene = path.to.tx2gene,
                                               metadata = input.metadata)
      
      path.to.save.output <- file.path(path.to.02.output, 
                                       input.set, 
                                       sprintf("%s_vs_%s", condition1, condition2), 
                                       sprintf("sizeFactor_adjust_%s_%s", hk.gene, explicit.size.factor))
      dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)
      
      deseq.output <- run_DESeq2_and_preprocess(deseq.dataset,
                                                tx2gene,
                                                thresh.pval = 0.05,
                                                geneid.to.remove = NULL,
                                                explicit.size.factor = explicit.size.factor,
                                                hk.gene = hk.gene)
      
      plot_result_from_deseq_output(deseq.output, path.to.save.output, input.metadata)
    }
  }
}

