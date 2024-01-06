# -*- coding: utf-8 -*-
"""
Created on Tue May 23 23:34:27 2023

@author: kevin
"""

import os
import pickle
import numpy as np
from xgboost.sklearn import XGBClassifier
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedGroupKFold
import pandas as pd
import sklearn
import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc
from sklearn.model_selection import GridSearchCV
from skopt import BayesSearchCV


#set tube and training file
tube = "tube1"
train_file = "training_file.csv" #name of file containing labeled training data
train_dir = "/input_path/" #insert path to folder with training files
output_dir = "/output_path/"

os.chdir(output_dir)

df_train = pd.read_csv(train_dir + train_file, index_col = 0)

#Training file should contain a column with Sample labels (named Sample.x here), Class labels, and feature columns

df_train["Sample.x"] = df_train["Sample.x"].astype('str') 

Y = df_train[["Class"]].values
Y[Y == "Normal"] = 0
Y[Y == "Blast"] = 1
Y = Y.astype(int)
X = df_train.iloc[:, 6:] #include only the desired feature columns for model training

group_k_fold= StratifiedGroupKFold(10)
gkf = group_k_fold.split(X, Y, groups = df_train["Sample.x"].astype("category").cat.codes)
gkf = [x for x in gkf]


XGB = XGBClassifier(tree_method = "gpu_hist", scale_pos_weight = 1)
optimization_dict = {'max_depth': [6,8,10], 'gamma' : [1,5], 'colsample_bytree' :[0.6,0.8,1], 
                 'eta' : [0.05, 0.3], 'n_estimators' : [75, 200], 
                 'min_child_weight': [1,5]}

model = BayesSearchCV(XGB, optimization_dict, 
                  scoring='roc_auc', cv=gkf, verbose=1)

model.fit(X,Y)

filename = tube + "_xgb.sav"
pickle.dump(model, open(filename, 'wb'))


