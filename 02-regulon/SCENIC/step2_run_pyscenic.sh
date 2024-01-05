echo "Start: `date`"
tf=${tf_file}
exprDat=${expression_matrix}

pfx=$prefix
grnOut=${pfx}.adjacencies.tsv
ctxOut=${pfx}.regulons.tsv
aucellOut=${pfx}.auc_mtx.tsv
dbDir=${config_file_outdir}
cores=6

# run grn
echo "Run GRNBoost"
pyscenic grn -o $grnOut --num_workers $cores --transpose -m grnboost2  $exprDat $tf

# run ctx
echo "Run RcisTarget"
pyscenic ctx --expression_mtx_fname $exprDat \
                         --transpose \
                         --num_workers $cores \
                         --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
                         --output $ctxOut \
                         $grnOut \
                         $dbDirhg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather,hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
                         

# run aucell
echo "Run AUCell"
pyscenic aucell --output=$aucellOut --num_workers=$cores --transpose $exprDat $ctxOut

echo "Done: `date`"