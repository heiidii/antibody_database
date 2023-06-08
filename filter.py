import os
from process_oas import process_files


#process_files(outpath='processed_pOAS_fastas', basename='pairedseqs')
#filters = {'Vaccine': 'None','Disease': 'None'}
#process_files(outpath='processed_pOAS_fastas_novac_nodis', basename='pairedseqs_novac_nodis',
#            filter_by=filters)
#filters = {'Vaccine': 'None','Disease': 'None','Species':'human'}
#process_files(outpath='processed_pOAS_fastas_human_novac_nodis', basename='pairedseqs_human_novac_nodis',
#            filter_by=filters)
#filters = {'Species':['human', 'rabbit']}
#process_files(outpath='processed_pOAS_fastas_human_rabbit', basename='pairedseqs_human_rabbit',
#            filter_by=filters)
#others 'HIV'
filters = {'Disease': 'SARS-COV-2'}
outpath='processed_pOAS_fastas_sarscov2'
process_files(outpath=outpath, basename='pairedseqs_sarscov2',
            filter_by=filters)
filters = {'Disease': 'HIV'}
outpath='processed_pOAS_fastas_hiv'
process_files(outpath=outpath, basename='pairedseqs_hiv',
            filter_by=filters)