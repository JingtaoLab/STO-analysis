#This text is for helping readers reproduce regulon analysis

#step1: Get gene expression matrix file and 3D coordinate file
python step1_data.py $in_h5ad $sample_name $outdir

#step1: pyscenic docker version is installed

#step2: Download config files for human from 
https://resources.aertslab.org/cistarget/

#step3: Open .sh file, editing the parmeters
#step4: Run pyscenic
sh step2_run_pyscenic.sh