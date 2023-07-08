import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from collections import defaultdict
import seaborn as sns
from sklearn.metrics import ConfusionMatrixDisplay
import pickle

from Data_Processing.EAG_DataProcessing_Library import open_wave


def name_con(f):
    n = os.path.basename(f)
    tn = n.split("_")
    if tn[1].lower() == '100':
        x = '1:100'
    elif tn[1].lower() == '1k':
        x = '1:1k'
    elif tn[1].lower() == '10k':
        x = '1:10k'
    else:
        x = ''
    name = x + ' ' + tn[2]
    return name


def csv_plot(FILE, NAME):
    t = open_wave(FILE)
    plt.title(label=NAME, size=10)
    plt.plot(t)

def Open_Pickle_Jar(filepath, perms='rb', array=True):
    file = open(filepath, perms)
    data = pickle.load(file)
    file.close()
    if array == True:
        data = np.asarray(data, dtype=object)  # this turns the current list of lists into a nested array (easier to index)
    return data



def pickle_to_DF(FILE):
    f=Open_Pickle_Jar(FILE)
    dd = defaultdict(list)

    for d in (f[0]): # you can list as many input dicts as you want here
        for key, value in d.items():
            dd[key].append(value)
    return pd.DataFrame.from_dict(dd)

def ClassifierResults(BASE, BUTTER, PROCESS,CH, QCTHRESH,Feature, MODEL):

    df=pickle_to_DF(f'{BASE}Butter{BUTTER}/{PROCESS}/{CH}/ClassifierResults/{Feature}/_QC_T_{QCTHRESH}/{MODEL}_Results.pickle')
    return df

def extract_CM(DF, cumulative=False, mean=True):
    if cumulative == True:
        sumCM=DF['confusion_matrix'].sum()
    else:
        sumCM=None
    if mean == True:
        meanCM=DF['confusion_matrix'].sum()/DF['confusion_matrix'].sum().sum()
    else:
        meanCM=None
    return sumCM, meanCM

def plot_CM(CM,LABELS,title):
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.xticks(fontsize=14,weight='bold', rotation=45)
    plt.tick_params(axis='both', rotation=45)
    plt.yticks(fontsize=14,weight='bold')
    plt.title(title ,fontsize=20, y=1.05, weight='bold')
    disp=ConfusionMatrixDisplay(CM,display_labels=LABELS)
    CMDISP=disp.plot(cmap=plt.cm.Reds,values_format='g',ax=ax)
    plt.xticks(rotation=45)
    plt.imshow(CM, cmap=plt.cm.Blues, vmin=0, vmax=1)
    return CMDISP

def ViPlot(data, title, num_odors, savedir):
    plt.figure(figsize=(18, 9))
    plt.xticks(rotation=0, fontsize=12, weight='bold')
    plt.yticks(fontsize=14, weight='bold')
    plt.xlabel('Dataset', fontsize=20, weight='bold')
    plt.ylabel('Mean Accuracy', fontsize=20, weight='bold')
    plt.title(title, fontsize=20)
    plt.ylim(0, 1)
    plt.axhline(y=(1 / num_odors), linestyle='--', color='black')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    # Create the violin plot
    ax = sns.violinplot(data=data, color="skyblue")

    # Calculate the mean for each column
    means = data.mean()

    # Add annotations to the plot
    for i, mean in enumerate(means):
        ax.text(i, mean, f'Mean: {mean:.2f}', ha='center', va='top', fontsize=12)

    plt.ylim(0, 1.5)
    plt.tight_layout()
    #plt.savefig(savedir)
    plt.show()
