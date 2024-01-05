#This text help reader to reproduce Hotspot gene module analysis
#The input data is .h5ad files which generate from our building pipeline

#step1: Get regulon matrix file and 3D coordinate file

The file with 3D coordinate was from SCENIC step1_get_data

#After running step1, you can get regulon matrix and 3D coordinate matrix

#step2: Run Hotspot
python step2_hotspot.py $regulon_matrix $coordinate_file $outdir $prefix

#step3: Transfer hotspot object to .csv files and plot module scores in 3D space
python $hotspot_object $outdir
