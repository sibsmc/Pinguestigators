#!/bin/python


import pandas as pd
import numpy as np
# from Models import DiseaseCAT

# Load datasets;
DATA_WITHSEQ = pd.read_csv("Data/GTEx_pancreas_liver_images_liverfat_pancreasfat_seq.csv")
DATA_ALL = pd.read_csv("Data/GTEx_pancreas_liver_images_liverfat_pancreasfat.csv")

CFAT_LIVER = "Fat.Percent_liver"
CFAT_PANC = "Fat.Percent_pancreas"

Models = []

# Select which dataset is going to be used.
DATASET = DATA_ALL


def SplitDataset(dataset, pct_test=0.1):
    """

    Splits the dataset randomly in a training and a testing part.

    """

    SIZE = dataset.shape[0]
    all_rows = range(SIZE)
    test_mask = np.random.choice(all_rows, size=round(pct_test * SIZE))

    TRAIN = dataset.iloc[~test_mask]
    TEST = dataset.iloc[test_mask]

    return TRAIN, TEST


def CalculateMSE(REAL, PRED):
    """

    Each of REAL and PRED come in as [(liver%, panc%)]...

    Output: [liver, pancMSE]
    """
    nb_cats = len(REAL[0])
    result = []

    for i in range(nb_cats):
        cat_results = []
        for r, p in zip(REAL, PRED):
            e = (REAL - PRED) ** 2
            cat_results.append(e)

        result.append(np.mean(cat_results))

    return result


DATASET_TRAIN, DATASET_TEST = SplitDataset(DATASET)


for Model in Models:
    # Initialize model;
    model = Model.Model(DATASET_TRAIN)

    REAL = []
    PREDICTED = []
    for ROW in DATASET_TEST.iloc():
        _ROW = ROW.copy

        # Fetch real fat % from patient record.
        real = (_ROW[CFAT_LIVER], _ROW[CFAT_PANC])
        REAL.append(real)

        # Delete result values... the model should not know them.
        del _ROW[CFAT_LIVER]
        del _ROW[CFAT_PANC]

        pred = model.Predict(_ROW)
        PREDICTED.append(pred)

        print("Real %s\nPred %s" % (real, pred))

    MSEs = CalculateMSE(REAL, PREDICTED)
    print("Model %s MSE:" % (model.name))

    print("Liver: %.2f" % MSEs[0])
    print("Pancreas: %.2f" % MSEs[1])

    print()
