#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sklearn
import sklearn.preprocessing
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression


noSeq = pd.read_csv("../Data/GTEx_pancreas_liver_images_liverfat_pancreasfat.csv")
"""
Focus only on pancreas data
"""


#collect the features and labels
featureNames = ["~Age", "Sex", 'saponification', 'atrophy', 'pancreatitis',
       'congestion', 'hemorrhage', 'scarring', 'inflammation', 'nodularity',
       'metaplasia', 'cyst', 'no_abnormalities', 'diabetic', 'desquamation',
       'clean_specimens', 'necrosis', 'fibrosis', 'calcification', 'sclerotic',
       'None']

x  = np.asarray(noSeq[featureNames])
y = np.asarray(noSeq["Fat.Percentage_pancreas"])

#train test split
x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=0.8, test_size=0.2)

model = LinearRegression()
model.fit(x_train, y_train)
predictions = model.predict(x_test)
predictions = predictions.reshape(-1,1)
y_test = y_test.reshape(-1,1)
assert predictions.shape == y_test.shape

pancreas_RMSE = (mean_squared_error(y_test, predictions))**.5
print(predictions, "\n\n\n", pancreas_RMSE)
#run in my notebook pancreas_MSE = 18.841485307059244
#on first pass without any tuning

"""
Same manual process, but now focusing on liver data
"""

featureNames = ["~Age", "Sex", 'necrosis', 'diabetic', 'pancreatitis',
       'nodularity', 'fibrosis', 'None', 'metaplasia', 'sclerotic',
       'saponification', 'scarring', 'congestion', 'atrophy',
       'no_abnormalities', 'cyst', 'hemorrhage', 'inflammation',
       'clean_specimens', 'calcification', 'desquamation']

x  = np.asarray(noSeq[featureNames])
y = np.asarray(noSeq["Fat.Percentage_liver"])

#train test split
x_train, x_test, y_train, y_test = train_test_split(x, y, train_size=0.8, test_size=0.2)

model = LinearRegression()
model.fit(x_train, y_train)
predictions = model.predict(x_test)

predictions = predictions.reshape(-1,1)
y_test = y_test.reshape(-1,1)
assert predictions.shape == y_test.shape

liver_RMSE = (mean_squared_error(y_test, predictions))**0.5
print(predictions, "\n\n\n", liver_RMSE)
#run in my notebook liver_MSE = 11.881679131462677
#on first pass without any tuning 
