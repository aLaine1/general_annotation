import sys
import os
import csv
import gzip
import datetime
import numpy as np
import pandas as pd

DataFrame = pd.core.frame.DataFrame

def loadTable(path: str, dtype=str, **kwargs) -> DataFrame:
    options = ""
    if ".tsv" in path:
        df = pd.read_csv(path, dtype=dtype, sep="\\t",engine='python',**kwargs)
    elif ".csv" in path:
        df = pd.read_csv(path, dtype=dtype, sep=",",engine='python',**kwargs)
    else:
        sys.exit("Invalid format of input file. Should be a TSV or CSV file (can be gzipped)")
    return df

def mergeAll(input_table,bam_and_gff,blast,chim,unique_id_col,column_to_keep,output_file):
    merged = mergeBAM_GFF(input_table,bam_and_gff,unique_id_col,column_to_keep,output_file)
    chim_loaded = loadTable(chim)
    merged=pd.merge(left=merged,right=chim_loaded,how="outer",on=unique_id_col)
    merged['mapped_to'] = np.where(pd.isna(merged['mapped_to']),merged["reference_chim"],merged["mapped_to"])
    merged.drop('reference_chim',axis=1,inplace=True)
    for blast_ali in blast:
        blast_loaded = loadTable(blast_ali)
        merged = pd.merge(left = merged,right = blast_loaded, how = "outer", on = unique_id_col)
    merged.to_csv(output_file,sep="\t", index = False,na_rep="NA")

def mergeBAM_GFF(input_table,bam_and_gff,unique_id_col,column_to_keep,output_file):
    if column_to_keep == "all":
        base_file = loadTable(input_table)
    else:
        base_file = loadTable(input_table,usecols=column_to_keep)
    bam_gff = loadTable(bam_and_gff)
    merged = pd.merge(left = base_file,right = bam_gff, how = "outer", on = unique_id_col)
    merged.fillna(value={"mapped_to":"None"})
    return merged
