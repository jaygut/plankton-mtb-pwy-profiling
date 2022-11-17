
'''Plotting a Venn diagram to visualize differences between e.g. set of active enzymes (E.C. numbers) in an environmental sample (i.e. MTX/MTG dataset) 
   mapped through the humann3 profiler using the default DBs (chocophlan + UniRef90) versus a custom protein DB assembled from a bunch of genomic assemblies 
   from a great variety of plankton species annotated using BUSCO + Augustus + Diamond (against UniRef90) 
'''

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns

from subprocess import Popen, call, STDOUT, PIPE
import os
import pandas as pd
import numpy as np
import matplotlib
import sys

from matplotlib_venn import venn2
from matplotlib_venn import venn3


if __name__ == "__main__":
    
    matplotlib.rcParams['savefig.dpi'] = 500
    matplotlib.rcParams['figure.dpi'] = 400
    sns.set_style("whitegrid", {'axes.grid' : False})
    sns.set_context("paper")
    sns.set(font='serif')
    # sns.set_style("white", {
    #         "font.family": "serif",
    #         "font.serif": ["Times", "Palatino", "serif"]
    #     })
    sns.set_style('ticks')

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    ec_set_merged = pd.read_csv(input_file)

    fig, ax = plt.subplots(figsize=(5,5))
    ax.set_title("Unique E.C. numbers mapped through HUMAnN 3.0")
    ec_set0 = set(ec_set_merged.query('DB == "Default"')["EC"].values)
    ec_set1 = set(ec_set_merged.query('DB == "Custom"')["EC"].values)
    venn2([ec_set0, ec_set1], set_labels = ('Default DBs', 'Custom Marine \nPlankton DB'),  ax = ax)
    fig.savefig(output_file,dpi=1200,bbox_inches="tight", pad_inches=0.1)