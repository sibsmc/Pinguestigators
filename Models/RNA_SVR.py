#!/bin/python

import pandas as pd
import sklearn
from sklearn.preprocessing import StandardScaler, MinMaxScaler

TissueSamples = pd.read_csv("Data/count_matrix_target_subset.tsv", sep="\t")


def GetGeneMask(kind, Descriptions):

    Pathways = pd.read_csv("Data/%s_pathways_genes.csv" % kind, sep=" ")

    def GetGenes (PathwaysRecord):
        allGenes = []
        SplitCell = lambda v: v.split(";")
        for i, row in PathwaysRecord.iterrows():
            allGenes += SplitCell(row["Genes"])
        return allGenes

    def GetMaskPiece(Gene):
        if Gene in AllGenes:
            return True
        else:
            return False

    AllGenes = GetGenes(Pathways)

    Mask = [GetMaskPiece(Gene) for Gene in Descriptions]

    return Mask


GeneMask = {
    "pancreas": GetGeneMask("pancreas", TissueSamples["Description"]),
    "liver": GetGeneMask("liver", TissueSamples["Description"])
}


def TissueSampleRow(kind):
    return "Tissue.Sample.ID_%s" % kind


def FatPercentRow(kind):
    return "Fat,Percentage_%s" % kind


class Model():
    name = "RNAseq SVR"
    def __init__(self, testing_dataset):

        self.modelPancreas, self.pscaler = BuildModel(testing_dataset, "pancreas")
        self.modelLiver, self.lscaler = BuildModel(testing_dataset, "liver")

    def Predict(self, Input):
        try:
            pancreas_base = [RowToModelInput(Input, "pancreas")]
            pancreas_fat = self.modelPancreas.predict(pancreas_base)

            liver_base = [RowToModelInput(Input, "liver")]
            liver_fat = self.modelLiver.predict(liver_base)

            return [liver_fat[0], pancreas_fat[0]]
        except (KeyError, ValueError) as e:
            print("Prediction error (missing RNAseq data for sample)....")
            return [20, 20]


def RowToModelInput(row, kind):


    SampleID = row[TissueSampleRow(kind)]

    TrueSampleIDs = [r for r in TissueSamples.columns
                     if r.startswith(SampleID)]

    if not TrueSampleIDs:
        return None

    TrueSampleID = TrueSampleIDs[0]
    try:
        sample = TissueSamples[[TrueSampleID]]

        Masked = sample[GeneMask[kind]]

        return Masked.values.reshape(-1,)
    except KeyError as e:
        print("Key error: %s" % SampleID)
        return None


def BuildModel(dataset, kind):
    def GetOutput(row):
        if kind in ["pancreas", "liver"]:
            return row[FatPercentRow(kind)]

    model = sklearn.svm.SVR(degree=2, kernel="poly", gamma="auto", tol=1e-7, max_iter=10e6)

    pre_inputs = [RowToModelInput(row, kind) for i, row in dataset.iterrows()]

    inputs = [i for i in pre_inputs if i is not None]

    outputs = [
        GetOutput(row)
        for i, (I, row) in enumerate(dataset.iterrows())
        if pre_inputs[i] is not None
    ]
    assert inputs

    # SCALER DISABLED!
    scaler = MinMaxScaler()
    scaler.fit(inputs)
    # inputs = scaler.transform(inputs)
    model.fit(inputs, outputs)

    return model, scaler
