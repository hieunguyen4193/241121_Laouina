# $ aws s3 cp --no-sign-request s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/README.txt .
# download: s3://ngi-igenomes/igenomes/Homo_sapiens/Ensembl/GRCh37/Annotation/README.txt to ./README.txt
 
# $ cat README.txt
# The contents of the annotation directories were downloaded from Ensembl on: July 17, 2015.
 
# Gene annotation files were downloaded from Ensembl release 75. SmallRNA annotation files were downloaded from miRBase release 21.

# #!/bin/bash
 
# VERSION=108
# wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
# wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/homo_sapiens/Homo_sapiens.GRCh38.$VERSION.gtf.gz

# VERSION=113
# wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
# wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz
nextflow="/home/hieunguyen/CRC1382/src_2023/241121_Laouina_Pabst_MolMedizin_3mRNAseq/nextflow"
main_nf="/home/hieunguyen/CRC1382/src_2023/241121_Laouina_Pabst_MolMedizin_3mRNAseq/rnaseq/main.nf";
outdir="/home/hieunguyen/CRC1382/outdir/41121_Laouina_Pabst_MolMedizin_3mRNAseq";
path_to_ref="/home/hieunguyen/CRC1382/outdir/ref";
work=${outdir}/work;
inputdir="/home/hieunguyen/Downloads/FASTQ/241121_Laouina_Pabst_MolMedizin_3mRNAseq";
mkdir -p ${work}
mkdir -p $outdir;
path_to_sample_sheet="/home/hieunguyen/CRC1382/src_2023/241121_Laouina_Pabst_MolMedizin_3mRNAseq/SampleSheet.csv"
gtf=${path_to_ref}/Mus_musculus.GRCm39.113.gtf;
fasta=${path_to_ref}/Mus_musculus.GRCm39.dna.primary_assembly.fa;
max_mem="75GB";
max_cpus=20;

${nextflow} run ${main_nf} \
--input ${path_to_sample_sheet} \
--outdir ${outdir} \
--fasta $fasta \
--gtf $gtf \
-profile docker \
-resume \
--max_cpus $max_cpus  \
--max_memory $max_mem \
--publish_dir_mode copy \
-w ${work} 
