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


def CheckResult(R):
    if 0 <= R <= 100:
        return R
    else:
        print("Bizarre result: %s" % R)
        return 10


class Model():
    name = "RNAseq SVR"
    def __init__(self, testing_dataset):

        self.modelPancreas, self.pscaler = BuildModel(testing_dataset, "pancreas")
        self.modelLiver, self.lscaler = BuildModel(testing_dataset, "liver")

    def Predict(self, Input):
        try:
            pancreas_base = [RowToModelInput(Input, "pancreas")]
            pancreas_fat = self.modelPancreas.predict(pancreas_base)[0]

            liver_base = [RowToModelInput(Input, "liver")]
            liver_fat = self.modelLiver.predict(liver_base)[0]

            return [CheckResult(liver_fat), CheckResult(pancreas_fat)]

        except (KeyError, ValueError) as e:
            print("Prediction error (missing RNAseq data for sample).... \n%s" % e)
            return [20, 20]


def RowToModelInput(row, kind):
    """

    This converts a patient row into inputs for the SVR.
    In this model we use RNAseq values as inputs.

    """
    SampleID = row[TissueSampleRow(kind)]

    TrueSampleIDs = [r for r in TissueSamples.columns
                     if r.startswith(SampleID)]

    if not TrueSampleIDs:
        return None

    TrueSampleID = TrueSampleIDs[0]
    assert len(TrueSampleIDs) <= 1

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

    model = sklearn.svm.SVR(degree=2, kernel="poly",
                            gamma="auto", tol=1e-7, max_iter=-1)

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
