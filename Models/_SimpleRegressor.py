#!/bin/python


import numpy as np
import pandas as pd
from itertools import combinations
import sklearn
import sklearn.preprocessing
from sklearn.metrics import r2_score
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression


class Model():
    name = "SimpleRegressor"

    def __init__(self, training_data):
        self.LiverModel = BuildModel("liver", training_data)
        self.PancreasModel = BuildModel("pancreas", training_data)

    def Predict(self, Individue):
        LIVER_X = [ExtractValues("liver", Individue)]
        PANCR_X = [ExtractValues("pancreas", Individue)]
        liver_fat = self.LiverModel.predict(LIVER_X)
        pancreas_fat = self.PancreasModel.predict(PANCR_X)

        return liver_fat, pancreas_fat


def SplitDiseasesField(field):
    # account for Nan categories
    try:
        return [f.strip() for f in field.split(",")]
    except:
        return ["None"]


def FlattenDiseases(kind, dataSet):
    """
    Some categories appear in groups of 2 or more. FlattenDiseases will separate the
    individual diseases that appear in the Pathology.Categories_... column into single,
    unique values

    Param:
    kind (str): either "pancreas" or "liver"
    dataSet (DataFrame): either the sequence file or the non sequence file

    Returns:
    (set) single_categories: a set of the diseases as unique values
    """

    disease_col = dataSet["Pathology.Categories_{}".format(
        kind)].dropna().unique()

    single_categories = ["None"]

    # assume dropped NaN values
    for rawCategory in disease_col:
        rawCategory = rawCategory.strip()
        categories = [rawCategory.split(",")]

        for diseases in categories:
            for d in diseases:
                single_categories.append(d.strip())

    return set(single_categories)


def ProcessDiseases(kind, dataSet):
    """
    Create a mapping between the `kind` fat percentage and every instance a disease occurs.
    Return this mapping as a dictionary
    {"disease": [fat percent at instance of "disease"]}

    Param:
    kind (str): either "pancreas" or "liver"
    dataSet (DataFrame): either the sequence file or the non sequence file

    Returns:
    (dictionary) prediction_struct: dictionary of the form {"disease": [fat percent at instance of "disease"]}

    Dependency: FlattenDiseases()
    """

    flat = FlattenDiseases(kind, dataSet)

    try:
        # select diseases and percentages
        categories = "Pathology.Categories_{}".format(kind)
        pcts = "Fat,Percentage_{}".format(kind)
        # create the dataframe with the specified columns
        df = dataSet[[categories, pcts]]
    except:
        # select diseases and percentages
        categories = "Pathology.Categories_{}".format(kind)
        pcts = "Fat.Percentage_{}".format(kind)
        # create the dataframe with the specified columns
        df = dataSet[[categories, pcts]]

    # fill any nan values with a string "None"
    #values = {categories: "None"}
    #df = df.fillna(value = values)

    # mapping between fat percents and records with each disease
    prediction_struct = {}

    for disease in flat:
        for i in range(0, len(df)):
            if disease in df.iloc[i][0]:
                prediction_struct.setdefault(disease, [df.iloc[i][1]])
                prediction_struct[disease].append(df.iloc[i][1])

    return prediction_struct


def CalculateFatForDisease(disease_values):
    disease_probs = {}
    for key in disease_values.keys():
        disease_probs[key] = np.mean(disease_values[key])

    return disease_probs


### Some final processing to feed into the model ###
def miniBatch(kind, dataSet, batchSize):
    """
    generate mini batches of all possible diseases for a data set
    """

    # all possible diseases for the dataset
    possibleDiseases = list(FlattenDiseases(kind, dataSet))

    batches = []
    for i in range(0, len(possibleDiseases)-batchSize, batchSize):
        batches.append(possibleDiseases[i:i+batchSize])

    return batches


def contains(a, b):
    """
    return if a contains b
    """
    return any(item in a for item in b)

# simple example with contains
#print(contains([1,2,3,4], [5, 1]))

# encode the pathology categories to numbers


def encodePathology(kind, dataSet, verbose=True):
    """
    modify dataSet in place to encode diseases

    Param:
    kind (str): either "pancreas" or "liver"
    dataSet (DataFrame): either the sequence file or the non sequence file
    verbose (boolean) [optional]: default to `True`. At the end of the method call,
    print a preview of the columns.

    Returns:
    None

    """


    # column names for
    pathology_categories_kind = 'Pathology.Categories_{}'.format(kind)

    # get the data for each new column
    mixed_categories_kind = dataSet[pathology_categories_kind]
    flat_categories_kind = list(FlattenDiseases(kind, dataSet))

    for name in flat_categories_kind:
        dataSet[name] = dataSet[pathology_categories_kind].apply(
            lambda x: 1 if name in flat_categories_kind else 0)

    if verbose:
        print(dataSet.columns)
        dataSet[flat_categories_kind].head(6)


def ExtractValues(kind, dataSet, featureNames=["~Age", "Sex", 'saponification', 'atrophy', 'pancreatitis',
                                               'congestion', 'hemorrhage', 'scarring', 'inflammation', 'nodularity',
                                               'metaplasia', 'cyst', 'no_abnormalities', 'diabetic', 'desquamation',
                                               'clean_specimens', 'necrosis', 'fibrosis', 'calcification', 'sclerotic',
                                               'None']):


    # Some features are not present, so we need this to avoid runtime IndexErrors;
    presentFeatureNames = [x for x in featureNames if x in dataSet.keys()]
    x = np.asarray(dataSet[presentFeatureNames])
    try:  # Maybe we are predicting... so we won't have this column available.
      y = np.asarray(dataSet["Fat.Percentage_{}".format(kind)])
    except:
      y = []

    return x, y


def BuildModel(kind, dataSet):
    # train test split
    df = encodePathology(kind, dataSet)
    x, y = ExtractValues(kind, dataSet)

    model = LinearRegression()
    model.fit(x, y)

    return model
