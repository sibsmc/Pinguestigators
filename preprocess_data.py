#!/bin/python


import pandas as pd
import os
import numpy as np

data_dir = "Data"


data_paths = [
    os.path.join(data_dir, f)
    for f in os.listdir(data_dir)
    if f.endswith(".xlsx")
]


def fix_age(data):
    old_colname = "Age.Bracket"
    new_colname = "~Age"
    data = data.rename(columns={old_colname: new_colname})

    def fix_age_cell(cell):
        V = [int(b) for b in cell.split("-")]
        return np.ceil(np.mean(V))

    data[new_colname] = data[new_colname].apply(fix_age_cell)

    return data


def fix_sex(data):
    sex_map = {
        "male": 1,
        "female": 0
    }

    data["Sex"] = data["Sex"].apply(lambda v: sex_map[v])

    return data


for file_path in data_paths:
    data = pd.read_excel(file_path, index=False)
    output_path = os.path.splitext(file_path)[0] + ".csv"

    data = fix_age(data)
    data = fix_sex(data)

    data.to_csv(output_path, index=False)
