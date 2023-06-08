import os
import sys
import json
import pandas as pd
import s3fs
from multiprocessing import Pool
import gzip
from pyarrow import Table, dataset
from pyarrow.parquet import ParquetDataset
from pyarrow.util import guid


HEAVY_ID_COL = 'sequence_id_heavy'
LIGHT_ID_COL = 'sequence_id_light'
HEAVY_AA_COL = 'sequence_alignment_aa_heavy'
LIGHT_AA_COL = 'sequence_alignment_aa_light'
COMBINED_AA_COL = 'sequence_heavy_and_light'
COMBINED_ID_COL = 'sequence_id_combined'
ID_SEP = ':'
OAS_columns = [HEAVY_ID_COL, LIGHT_ID_COL, HEAVY_AA_COL, LIGHT_AA_COL]


def _clean_string(bad_string):
    return ''.join([s.rstrip() for s in str(bad_string).rstrip() if not s in 
                ['}', '{', '\"', '\\', "'"]])


def _read_metadata(metadata_line):
    metadata_line = _clean_string(metadata_line)
    meta_dict = {_clean_string(dat.split(':')[0]):_clean_string(dat.split(':')[1]) 
                 for dat in metadata_line.split(',')
                    if len(dat.split(':'))==2}

    print(meta_dict)
    return meta_dict

def _get_combined_id(df, file_id):
    return file_id + ID_SEP + df[HEAVY_ID_COL].astype(str) + ID_SEP  \
        + df[LIGHT_ID_COL].astype(str)

def _split_combined_id(combined_id):
    dat = combined_id.split(ID_SEP)
    return {'file_id': dat[0], HEAVY_ID_COL: dat[1], LIGHT_ID_COL: dat[2]}


def csvs_to_fastas(csv_path, filter_by=None):

    file_id=os.path.basename(csv_path).split('_')[0]
    if not filter_by is None:
        fs = s3fs.S3FileSystem()
        with fs.open(csv_path,'rb') as f:
            g = gzip.GzipFile(fileobj=f)
            filter_row = g.readline()
            print(filter_row)
        meta_dict = ','.join(pd.read_csv(csv_path, nrows=0).columns)
        meta_dict = json.loads(meta_dict)
        #meta_dict = _read_metadata(filter_row)
        for key, value in filter_by.items():
            if key in meta_dict:
                if isinstance(value,str): 
                	if meta_dict[key].find(value)==-1:
                    		#skip this file
                            return None
                elif isinstance(value,list):
                    found=False
                    for v in value:
                     found = found or (meta_dict[key].find(v)!=-1)
                    if not found: 
                     return None
                else:
                     sys.exit('Not implemented: ', value)
    
    df = pd.read_csv(csv_path, skiprows=[0], usecols=OAS_columns, low_memory=True)
    print(len(df[OAS_columns[0]]))
    if df.empty:
        return None
    df[COMBINED_AA_COL] = df[HEAVY_AA_COL]+ df[LIGHT_AA_COL]
    df[COMBINED_ID_COL] = _get_combined_id(df, file_id)
    fasta_lines = [f'>{seq_id}\n{seq}\n' for seq_id, seq in zip(df[COMBINED_ID_COL].tolist(),
                                                          df[COMBINED_AA_COL].tolist())]
    print(csv_path, len(fasta_lines))
    return ''.join(fasta_lines), file_id


def process_files(s3_path="s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
                 max_files=None,
                 outpath='./',
                 basename='pairedseqs',
                 filter_by=None):
    os.makedirs(outpath, exist_ok=True)
    s3 = s3fs.S3FileSystem()
    files = s3.ls(s3_path)
    files.sort()
    files_process = files if max_files is None else files[:max_files]
    j=0
    for i, file in enumerate(files_process):
        output = csvs_to_fastas(csv_path=f"s3://{file}", filter_by=filter_by)
        if output is None:
            continue
        data, file_id = output
        split_fasta_name = f"{outpath}/{basename}_{file_id}_{j:04d}.fasta"
        print('Writing ', i, len(files_process))
        with open(split_fasta_name, 'w') as f:
            for lines in data:
                f.write(lines)
        j+=1


def process_files_parallel(s3_path="s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
                            max_files=8):
    s3 = s3fs.S3FileSystem()
    files = [f's3://{path}' for path in s3.ls(s3_path)[:max_files]]
    combined_fasta = 'out.fasta'
    print(files)

    with Pool(4) as pool_run:
        data, file_id = pool_run.map(csvs_to_fastas, files)

    with open(combined_fasta, 'w') as f:
        for lines in data:
            f.write(lines)


def linclust_to_df(linclust_fasta, s3_path, out_path):
    fasta_lines = open(linclust_fasta).readlines()
    seq_ids = [line.rstrip()[1:] for line in fasta_lines
                    if line.startswith('>')]
    file_to_heavy_map, file_to_light_map = {}, {}
    for seq_id in seq_ids:
        split_id = _split_combined_id(seq_id)
        file_id = split_id['file_id']
        if file_id in file_to_heavy_map:
            file_to_heavy_map[file_id].append(split_id[HEAVY_ID_COL])
            file_to_light_map[file_id].append(split_id[LIGHT_ID_COL])
        else:
            file_to_heavy_map[file_id] = [split_id[HEAVY_ID_COL]]
            file_to_light_map[file_id] = [split_id[LIGHT_ID_COL]]

    s3 = s3fs.S3FileSystem()
    files = s3.ls(s3_path)
    files.sort()
    print(out_path)
    os.makedirs(out_path, exist_ok=True)
    for file in files:
        print(file)
        csv_file = f"s3://{file}"
        file_id = os.path.basename(csv_file).split('_')[0]
        if file_id in file_to_heavy_map:
            meta_data = ','.join(pd.read_csv(csv_file, nrows=0).columns)
            meta_data = json.loads(meta_data)
            df = pd.read_csv(csv_file, skiprows=[0], low_memory=True)
            df_filter =  df[df[HEAVY_ID_COL].isin(file_to_heavy_map[file_id]) &
                                df[LIGHT_ID_COL].isin(file_to_light_map[file_id])]
            df_filter['Species'] = [meta_data['Species'] for _ in df_filter[HEAVY_ID_COL]]
            df_filter['Disease'] = [meta_data['Disease'] for _ in df_filter[HEAVY_ID_COL]]
            df_filter['Vaccine'] = [meta_data['Vaccine'] for _ in df_filter[HEAVY_ID_COL]]
            outfile = os.path.basename(file).split('.')[0] + '_filtered.csv.gz'
            df_filter.to_csv(path_or_buf = f'{out_path}/{outfile}', compression="gzip", index = None, sep = ",",
                                 header=True, encoding='utf-8-sig')

#adapted from walk-jump-sampler code
SALT_DEFAULT = guid()
def pandas_to_parquet(
    df: pd.DataFrame,
    path: str,
    existing_data_behavior: str = "delete_matching",
    partition_cols=None,
    salt: str = SALT_DEFAULT,
):
    """Write the dataframe as parquet to the specified path.

    The data will be split into multiple files, each one containing at most 65k rows.
    """
    table = Table.from_pandas(df)

    partitioning = None
    if partition_cols:
        part_schema = table.select(partition_cols).schema
        partitioning = dataset.partitioning(part_schema, flavor="hive")

    basename_template = salt + "-{i}.parquet"

    parquet_format = dataset.ParquetFileFormat()

    dataset.write_dataset(
        table,
        path,
        partitioning=partitioning,
        use_threads=True,
        format=parquet_format,
        basename_template=basename_template,
        existing_data_behavior=existing_data_behavior,
        max_rows_per_file=1 << 16,
        max_rows_per_group=1 << 16,
    )


def parquet_to_pandas(path: str, columns=None, **kwargs) -> pd.DataFrame:
    return (
        ParquetDataset(path, use_legacy_dataset=False, **kwargs)
        .read(columns=columns)
        .to_pandas(self_destruct=True)
    )

def annotate_df_by_file_and_seq_ids(df_path, s3_path, out_path):
    print('Untested')
    df = parquet_to_pandas(df_path)
    file_id = os.path.basename(df_path).split('_')[0]
    print(file_id)
    s3 = s3fs.S3FileSystem()
    files = s3.glob(s3_path + f'/{file_id}*')
    print(s3_path + f'/{file_id}_*')
    #print(files)
    assert len(files) == 1
    s3_file = files[0]
    meta_dict = ','.join(pd.read_csv(s3_file, nrows=0).columns)
    meta_dict = json.loads(meta_dict)
    df['Species'] = meta_dict['Species']
    os.makedirs(out_path, exist_ok=True)
    outfile = out_path + '/' +os.path.basename(df_path).split('.')[0] + '_with_species.csv'
    #df.to_csv(path_or_buf = f'{out_path}/{outfile}', compression="gzip", index = None, sep = ",",
    #         header=True, encoding='utf-8-sig')
    SPLIT_COLUMN = ['partition']
    pandas_to_parquet(df, outfile, partition_cols=SPLIT_COLUMN)


if __name__ == '__main__':
    print('Not implemented')
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
    #filters = {'Disease': 'SARS-COV-2'}
    #linclust_to_df('linclust_out/pOAS_hlcombined_novac_nodis_linclust_seqid95c80_rep_seq.fasta',
    #               "s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
    #               "/gstore/scratch/u/mahajs17/repositories/antibody_database/linclust_naive")
    #linclust_to_df('linclust_out/pOAS_hlcombined_linclust_seqid95c80_rep_seq.fasta',
    #               "s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
    #               "/gstore/scratch/u/mahajs17/repositories/antibody_database/linclust_all")
    #linclust_to_df('linclust_out/pOAS_hlcombined_human_novac_nodis_linclust_seqid95c80_rep_seq.fasta',
    #               "s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
    #               "/gstore/scratch/u/mahajs17/repositories/antibody_database/linclust_naive_human") 
    #linclust_to_df('linclust_out/pOAS_hlcombined_human_rabbit_linclust_seqid95c80_rep_seq.fasta',
    #               "s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
    #               "/gstore/scratch/u/mahajs17/repositories/antibody_database/linclust_naive_human_rabbit")
    #annotate_df_by_file_and_seq_ids("s3://prescient-data-dev/sandbox/mahajs17/walk-jump/data/parquet/large_molecules/linclust_poas_naive/",
    #                "s3://prescient-data-dev/sandbox/freyn6/raw/poas/",
    #              "/gstore/scratch/u/mahajs17/repositories/antibody_database/linclust_naive_human_annotated")
