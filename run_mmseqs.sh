
OUTDIR=linclust_out
FASTAS=../processed_pOAS_fastas
NAME=pOAS_hlcombined_linclust_seqid95c80
cd $OUTDIR; mmseqs easy-linclust \
$FASTAS/*.fasta $NAME  /tmp/ --min-seq-id 0.95 -c 0.8 --cov-mode 1

OUTDIR=linclust_out
FASTAS=../processed_pOAS_fastas_novac_nodis
NAME=pOAS_hlcombined_novac_nodis_linclust_seqid95c80
cd $OUTDIR; mmseqs easy-linclust \
$FASTAS/*.fasta $NAME  /tmp/ --min-seq-id 0.95 -c 0.8 --cov-mode 1
