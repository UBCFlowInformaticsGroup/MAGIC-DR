# -*- coding: utf-8 -*-
"""
Created on Wed Feb  1 17:00:09 2023

@author: kevin
"""
import shap
import os
import pickle
import numpy as np
from xgboost.sklearn import XGBClassifier
from sklearn.metrics import confusion_matrix
import pandas as pd
import sklearn
import matplotlib.pyplot as plt
import matplotlib as mpl
import scanpy as sc
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

np.random.seed(999)

sc.set_figure_params(scanpy=True, dpi=300, dpi_save=300, frameon=True, vector_friendly=True, fontsize=14, figsize=None, color_map=None, format='pdf', facecolor=None, transparent=False, ipython_format='png2x')

# Get the jet colormap
jet = plt.cm.get_cmap('jet')

# Define the new colormap
newcolors = jet(np.linspace(0, 0.9, 256))  # We select until 0.9 to get to bright red
newcmp = ListedColormap(newcolors)

#set tube and load files
tube = "tube1"
query_file = "query_case.csv" #file containing compensated flow data of query case(s) after extended blast gating
query_dir = "/input_dir/"
train_file = "training_data.csv" #file containing training data (will be used to anchor UMAP)

output_dir = "/output_dir/" + tube + "/"
os.chdir(output_dir)
df_query = pd.read_csv(output_dir + query_file, index_col = 0)

df_train = pd.read_csv("path_to_train_file" + train_file, index_col = 0)
df_train = df_train.loc[df_train["Class"] == "Normal", :].sample(100000) #subsample 100000 normal events from training data to anchor UMAP

model_file = "/path_to_model/" + tube + "_xgb.sav"
model = pickle.load(open(model_file, 'rb'))

#Parameters for tube 1 model: SSC.A, CD7, CD56, CD34, CD117, CD33, CD38, HLA.DR, CD13, CD45, CD133, CD2
#Parameters for tube 2 model: SSC.A, CD34, CD117, CD33, CD38, HLA.DR, CD13, CD45, CD19, CD22, CD15, CD10
#Parameters for tube 3 model: SSC.A, CD36, CD14, CD11b, CD34, CD117, CD64, HLA.DR, CD13, CD45, CD4, CD33, CD123


#Adjust column names as required depending on input format to match the formatting below


if tube == "tube1":
    df_query = df_query[['SSC.A', 'CD7', 'CD56', 'CD34', 'CD117', 'CD33', 'CD38',
                         'HLA.DR', 'CD13', 'CD45', 'CD133', 'CD2']]

if tube == "tube2":
    df_query = df_query[['SSC.A', 'CD15', 'CD22', 'CD34', 'CD117', 'CD19', 'CD38',
                         'HLA.DR', 'CD13', 'CD45', 'CD10', 'CD33']]

if tube == "tube3":
    df_query = df_query[['SSC.A', 'CD36', 'CD14', 'CD11b', 'CD34', 'CD117', 'CD64', 
                    'HLA.DR', 'CD13', 'CD45', 'CD4', 'CD33', 'CD123']]




df_query = df_query.drop("FSC-A", axis = 1)





#If query file contains multiple samples, this will loop through each sample
samples = df_query.index.str.split(pat = "_").str.get(1) #assumes each row is labeled SampleX_event#

for sample in samples.unique():
    #Create train and test sets
    
 
    df_test = df_query.loc[samples == sample, :]
   
    
    #Compute metrics
    preds = model.predict(df_test)
    probs = model.predict_proba(df_test)
    blast_prob = probs[:,1]
    df_test["Blast.Prob"] = blast_prob

    
    
    if tube == "tube1":
       channels = ["CD7", "CD56", "CD34", "CD117", "CD33", "CD38", "HLA.DR", 
                            "CD13", "CD45", "CD133", "CD2",'SSC.A']
    if tube == "tube2":
       channels = ["CD34", "CD117", "CD33", "CD38", "HLA.DR", 
                            "CD13", "CD45", "CD19", "CD22", "CD15", "CD10", "SSC.A"]
    if tube == "tube3":
       channels = ["CD34", "CD117", "CD33", "CD14", "HLA.DR", "CD11b",
                             "CD13", "CD45", "CD64", "CD36", "CD4", "CD123", "SSC.A"]
   
   
   
   
    X_test = df_test[channels]
    
    #output SHAP plots for high proability predicted blast events
    explainer = shap.TreeExplainer(model = model.best_estimator_)
    shap_values = explainer(X_test)
    top_shap = shap_values[blast_prob >= 0.85]
   
    
    plt.clf()
    shap.summary_plot(top_shap,show=False)
    plt.savefig(sample+"_shap.png", dpi = 300)
    
    #Prpeare data for UMAP
    adata_test = sc.AnnData(X_test)
    adata_test.obsm["transformed"] = np.asanyarray(X_test)

    
    adata_test.obs["Blast.Prob"] = df_test["Blast.Prob"].values
   
    
    X_normal = df_train[channels]
    adata_normal = sc.AnnData(X_normal)
    adata_normal.obsm["transformed"] = np.asanyarray(X_normal)
    
    
    
    adata_merge = adata_normal.concatenate(adata_test)
    
  
    
    sc.pp.neighbors(adata_merge, n_neighbors=20, n_pcs=12, use_rep = "transformed", random_state = 33 )
    sc.tl.umap(adata_merge,min_dist = 0.5, random_state = 23, maxiter = 500, )
    
    file = sample+"_prob"
    sc.pl.umap( adata_merge, color=["Blast.Prob"], add_outline=False,
           legend_fontsize=10,frameon=False, sort_order = False,
           title=' ',alpha = .8, color_map = newcmp, vmin = 0, vmax = 1, s = 5, save = file)

 
    file = sample + "_heat"    
    sc.pl.umap(adata_merge, color = channels,
           palette ='Set2',color_map=mpl.cm.inferno, frameon = False, size = 1, s = 6, ncols = 3, vmin = -0.2, vmax = 4, save = file)
    
    
    X_test["UMAP.x"] = adata_merge.obsm["X_umap"][len(X_normal):,0]
    X_test["UMAP.y"] = adata_merge.obsm["X_umap"][len(X_normal):,1]
    X_test["Blast.Prob"] = df_test["Blast.Prob"]
    
    
    X_test.to_csv(sample + "_output.csv")
    
   
    

