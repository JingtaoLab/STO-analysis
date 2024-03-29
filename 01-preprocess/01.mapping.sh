###download gtf and fasta from ensembl
wget https://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-93/gtf/homo_sapiens/Homo_sapiens.GRCh38.93.gtf.gz 

###install singularit and download SAW from DockerHub
singularity build SAW_5.4.sif docker://stomics/saw:05.4.0

###build index
singularity exec SAW_5.4.sif mapping --runMode genomeGenerate \
    --genomeDir reference/STAR_SJ100 \
    --genomeFastaFiles reference/genome.fa \
    --sjdbGTFfile reference/genes.gtf \
    --sjdbOverhang 99 \
    --runThreadN 12


####mapping
bash ./SAW/stereoPipeline.sh \
    -sif SAW_5.4.sif \
    -splitCount 1 \  ## 16 or 64 for Q4, 1 for Q40
    -maskFile SN.h5 \
    -fq1 $dataDir/reads/lane1_read_1.fq.gz,...,$dataDir/reads/laneN_read_1.fq.gz  \
    -fq2 $dataDir/reads/lane1_read_2.fq.gz,...,$dataDir/reads/laneN_read_2.fq.gz \ 
    -speciesName <speciesName> \ ###human mouse
    -tissueType <tissueName> \   ####embryo
    -refIndex reference/STAR_SJ100 \
    -annotationFile reference/genes.gtf \  ## GFF or GTF
    -threads 16 \
    -outDir $outDir/result \
    -imageRecordFile $dataDir/image/<SN_date_time_version>.ipr \ # [optional] when image is given and has passed QC
    -imageCompressedFile $dataDir/image/<SN_date_time_version>tar.gz \ # [optional] when image is given and has passed QC
    -doCellBin Y  # [optional] when you want to do the cell segmentation and get cell gene expression data
