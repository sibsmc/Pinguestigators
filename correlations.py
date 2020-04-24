#!/bin/python


import pandas as pd
import numpy as np

DATA_WITHSEQ = pd.read_csv("Data/GTEx_pancreas_liver_images_liverfat_pancreasfat_seq.csv")
DATA_ALL = pd.read_csv("Data/GTEx_pancreas_liver_images_liverfat_pancreasfat.csv")



fats = ["Fat.Percentage_pancreas", "Fat.Percentage_liver"]


dsall_correlations = DATA_ALL.corr()
dseq_coorelations = DATA_WITHSEQ.corr()

print("Full dataset correlations:")
print(dsall_correlations)
print()

print("SEQ dataset correlations:")
print(dseq_coorelations)
