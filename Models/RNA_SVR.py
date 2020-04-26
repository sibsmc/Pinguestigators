#!/bin/python

import pandas as pd
import sklearn

TissueSamples = pd.read_csv("Data/count_matrix_target_subset.tsv", sep="\t")


def GetGeneMask(Descriptions):
    return [True for x in Descriptions]

GeneMask = GetGeneMask(TissueSamples["Description"])


def TissueSampleRow(kind):
    return "Tissue.Sample.ID_%s" % kind


def FatPercentRow(kind):
    return "Fat,Percentage_%s" % kind


class Model():
    name = "RNAseq SVR"
    def __init__(self, testing_dataset):

        self.modelPancreas = BuildModel(testing_dataset, "pancreas")
        self.modelLiver = BuildModel(testing_dataset, "liver")

    def Predict(self, Input):
        try:
            pancreas_fat = self.modelPancreas.predict([RowToModelInput(Input, "pancreas")])
            liver_fat = self.modelLiver.predict([RowToModelInput(Input, "liver")])

            return [liver_fat[0], pancreas_fat[0]]
        except (KeyError, ValueError):
            print("Prediction error (missing RNAseq data for sample)....")
            return [20, 20]


def RowToModelInput(row, kind):

    try:
        SampleID = row[TissueSampleRow(kind)]

        TrueSampleIDs = [r for r in TissueSamples.columns
                         if r.startswith(SampleID)]

        if not TrueSampleIDs:
            return None

        TrueSampleID = TrueSampleIDs[0]

        sample = TissueSamples[[TrueSampleID]]

        return sample[GeneMask].values.reshape(-1,)
    except KeyError:
        print("Key error: %s" % SampleID)
        return None


def BuildModel(dataset, kind):
    def GetOutput(row):
        if kind in ["pancreas", "liver"]:
            return row[FatPercentRow(kind)]

    model = sklearn.svm.SVR(kernel="linear")
    print(dataset.shape[0])
    pre_inputs = [RowToModelInput(row, kind) for i, row in dataset.iterrows()]

    inputs = [i for i in pre_inputs if i is not None]

    outputs = [
        GetOutput(row)
        for i, (I, row) in enumerate(dataset.iterrows())
        if pre_inputs[i] is not None
    ]

    model.fit(inputs, outputs)

    return model
