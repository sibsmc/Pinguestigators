#!/bin/python


import pandas as pd
import numpy as np
import copy
from Models import DiseaseCATMean, RNA_SVR

from sklearn import metrics

# Load datasets;
DATA_WITHSEQ = pd.read_csv("Data/GTEx_pancreas_liver_images_liverfat_pancreasfat_seq.csv")
DATA_ALL = pd.read_csv("Data/GTEx_pancreas_liver_images_liverfat_pancreasfat.csv")


CFAT_LIVER = "Fat,Percentage_liver"
CFAT_PANC = "Fat,Percentage_pancreas"

Models = [DiseaseCATMean, RNA_SVR]

# Select which dataset is going to be used.
DATASET = DATA_WITHSEQ


def SplitDataset(dataset, pct_test=10):
    """

    Splits the dataset randomly in a training and a testing part.

    """

    SIZE = dataset.shape[0]
    TEST_SIZE = pct_test / SIZE * 100

    all_rows = range(SIZE)
    test_mask = np.random.choice(all_rows, size=round(TEST_SIZE))

    TRAIN = dataset.drop(test_mask, axis=0)
    TEST = dataset.iloc[test_mask]

    return TRAIN, TEST


DATASET_TRAIN, DATASET_TEST = SplitDataset(DATASET)
print()
print("-- Pinguestigators FAT predictors --")
print("Total dataset size: %i" % DATASET.shape[0])
print("Training dataset size: %i" % DATASET_TRAIN.shape[0])
print("Testing dataset size: %i" % DATASET_TEST.shape[0])
print()

for Model in Models:
    # Initialize model;
    model = Model.Model(DATASET_TRAIN)

    REAL = []
    PREDICTED = []

    for i, ROW in DATASET_TEST.iterrows():

        _ROW = ROW.copy()

        # Fetch real fat % from patient record.
        real = [_ROW[CFAT_LIVER], _ROW[CFAT_PANC]]
        REAL.append(real)

        # Delete result values... the model should not know them.
        del _ROW[CFAT_LIVER]
        del _ROW[CFAT_PANC]

        pred = model.Predict(_ROW)
        PREDICTED.append(pred)

        print("Real %s\nPred %s\n\n" % (real, pred))


    MSEs = metrics.mean_squared_error(REAL, PREDICTED, multioutput="raw_values")

    print("Model %s MSE:" % (model.name))
    print("Liver: %.2f" % MSEs[0])
    print("Pancreas: %.2f" % MSEs[1])

    print()
