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
    df = pd.read_csv(path, dtype=dtype, sep="\\t",engine='python',header=None, **kwargs)
    return df

def extractChimFromFile(input_file,output_file,unique_id_col,reference):
    base_chim = loadTable(input_file)
    output_index = [unique_id_col,"type","seg1_cj","seg2_cj","reference_chim"]
    CHIM_data_raw = []
    for index, row in base_chim.iterrows():
        tag = row[9]
        type = "chimeric_junction"
        chr1 = row[0]
        strand1 = row[2]
        chr2 = row[3]
        strand2 = row[5]
        tag = row[9]
        start1 = row[10]
        start2 = row[12]
        if chr1 == chr2 and strand1 == strand2:
            if strand1 == "+":
                if row[10] > row[12]:
                    type = "circRNA"
            elif strand1 == "-":
                if row[12] > row[10]:
                    type = "circRNA"
        seg_1 = ":".join([chr1,start1,strand1])
        seg_2 = ":".join([chr2,start2,strand2])
        new_line = [tag,type,seg_1,seg_2,reference+"(chimeric)"]
        CHIM_data_raw.append(new_line)
    CHIM_data = pd.DataFrame(CHIM_data_raw,columns = output_index)
    CHIM_data.to_csv(output_file,sep="\t",index = False)
