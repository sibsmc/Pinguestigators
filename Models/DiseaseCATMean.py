#!/bin/python


import pandas as pd
import numpy as np

from sklearn.linear_model import LinearRegression

PCAT_PAN = "Pathology.Categories_pancreas"
PCAT_LIV = "Pathology.Categories_liver"

FPCT_PAN = "Fat,Percentage_pancreas"
FPCT_LIV = "Fat,Percentage_liver"


class Model():
    """

    This model computes fat percents for liver and pancreas (independently)
    based on the mean fat percents observed for each single disease present in the
    "Pathology.Categories_XX" rows. If such categories is not available for
    given record, it will fall back to age alone as the criteria.

    """
    name = "Disease CAT Mean"

    def __init__(self, training_data):
        self.BANKS, self.LRMODELS = BuildModel(training_data)

    def Predict(self, Individue):

        age = Individue["~Age"]
        pancreas_fat = self._PredictEach(Individue[PCAT_PAN],
                                         age,
                                         self.BANKS["DISEASE_PANCREAS"],
                                         self.LRMODELS["PANCREAS"])

        liver_fat = self._PredictEach(Individue[PCAT_LIV],
                                      age,
                                      self.BANKS["DISEASE_LIVER"],
                                      self.LRMODELS["LIVER"])

        return [liver_fat, pancreas_fat]

    def _PredictEach(self, disease, age, BANK, fallbacklr_model):
        Diseases = SplitDiseaseField(disease)
        values = []
        if Diseases:
            values = [BANK[d] for d in Diseases if d in BANK.keys()]

        if not values:
            values = fallbacklr_model.predict([[age]])

        return np.mean(values)


def SplitDiseaseField(field):
    try:
        return [f.strip() for f in field.split(",")]
    except:
        return []


def CalculateFatForDisease(disease_values):
    disease_probs = {}
    for key in disease_values.keys():
        disease_probs[key] = np.mean(disease_values[key])

    return disease_probs


def UpdateBank(bank, diseases, fat_percent):

    for d in diseases:
        if d not in bank.keys():
            bank[d] = []

        bank[d].append(fat_percent)


def BuildModel(training_data):
    """

    Builds data structures used by the model to predict:
    FAT probabilities associated with each single disease.

    """

    # A dict of disease:[float]

    BANKS = {
        "DISEASE_PANCREAS": {},
        "DISEASE_LIVER": {},
        "AGE_PANCREAS": {},
        "AGE_LIVER": {}
    }

    MODELS = {
        "LIVER": [[], []],
        "PANCREAS": [[], []]
    }

    def model_fit(X, Y):
        m = LinearRegression()
        m.fit(np.array([X]).reshape(-1, 1), np.array([Y]).reshape(-1, 1))
        return m

    for i, row in training_data.iterrows():

        dp = SplitDiseaseField(row[PCAT_PAN])
        dl = SplitDiseaseField(row[PCAT_LIV])

        fat_pancreas = row[FPCT_PAN]
        fat_liver = row[FPCT_LIV]
        age = row["~Age"]
        UpdateBank(BANKS["DISEASE_PANCREAS"], dp, fat_pancreas)
        UpdateBank(BANKS["DISEASE_LIVER"], dl, fat_liver)


        MODELS["LIVER"][0].append(age)
        MODELS["LIVER"][1].append(fat_liver)
        MODELS["PANCREAS"][0].append(age)
        MODELS["PANCREAS"][1].append(fat_pancreas)

        # UpdateBank(BANKS["AGE_PANCREAS"], [age], fat_pancreas)
        # UpdateBank(BANKS["AGE_LIVER"], [age], fat_liver)

    return {k: CalculateFatForDisease(BANKS[k]) for k in BANKS.keys()}, {k: model_fit(*MODELS[k]) for k in MODELS.keys()}


