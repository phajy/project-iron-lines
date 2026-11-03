import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO

def ParseQDP(filename):

    with open(filename) as f:
        fileString = f.read()
        
        fileString = fileString[fileString.find("!")+2:]
        numcols = fileString[0:fileString.find("\n")].count(" ")+1
        splitter = "NO "*numcols + "\n"

        strings = fileString.split(splitter)

        dfs = []

        for (i, string) in enumerate(strings):
            dfs[i] = pd.read_csv(StringIO(string), sep=" ")
            
            
        
        

ParseQDP("/home/brad/OneDrive/UniWork/Year4/Dissertation/0693781301_ratio_plot_lp_ph_pl_ga.qdp")