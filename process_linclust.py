from process_oas import linclust_to_df


#linclust_to_df('linclust_out/pOAS_hlcombined_novac_nodis_linclust_seqid95c80_rep_seq.fasta',
#               "s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
#               "/gstore/scratch/u/mahajs17/repositories/antibody_database/linclust_naive")
#linclust_to_df('linclust_out/pOAS_hlcombined_linclust_seqid95c80_rep_seq.fasta',
#               "s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
#               "/gstore/scratch/u/mahajs17/repositories/antibody_database/linclust_all")
#linclust_to_df('linclust_out/pOAS_hlcombined_human_novac_nodis_linclust_seqid95c80_rep_seq.fasta',
#               "s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
#               "/gstore/scratch/u/mahajs17/repositories/antibody_database/linclust_naive_human")
linclust_to_df('linclust_out/pOAS_hlcombined_human_rabbit_linclust_seqid95c80_rep_seq.fasta',
                   "s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
                   "/gstore/scratch/u/mahajs17/repositories/antibody_database/linclust_naive_human_rabbit") 
