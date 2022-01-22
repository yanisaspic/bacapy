"""
Librairie pour aggreger des donnees quantitatives.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

import math
import pandas as pd
from scipy.stats import zscore

def getDataFrameRowByStd(dataFrame, ascend):
    """Renvoie la ligne d'un dataframe d'apres sa valeur d'ecart-type."""
    selectedRow = dataFrame.std(axis=1).sort_values(ascending=ascend).index[0]
    return dataFrame.loc[selectedRow]
def getDataFrameMinStdRow(dataFrame):
    """Renvoie la ligne d'un dataframe avec le plus petit ecart-type."""
    return getDataFrameRowByStd(dataFrame, True)
def getDataFrameMaxStdRow(dataFrame):
    """Renvoie la ligne d'un dataframe avec le plus grand ecart-type."""
    return getDataFrameRowByStd(dataFrame, False)

def countPosAndNeg(series):
    """Renvoie le nombre de valeurs positives et negatives dans une Series pandas."""
    return len(series[series>0]), len(series[series<0])
def getUpDownZScore(nUp, nDown):
    """Renvoie un score Z centre autour de 0.
    
    @ Marini, F., Ludt, A., Linke, J. et al. 
    GeneTonic: an R/Bioconductor package for streamlining the interpretation of RNA-seq data. 
    BMC Bioinformatics 22, 610 (2021). 
    https://doi.org/10.1186/s12859-021-04461-5
    """
    # ZeroDivisionError si nUp==0 et nDown==0 => utilise max(1).
    return (nUp-nDown) / math.sqrt(max(1, nUp+nDown))
def getDataFrameUpDownZAggregate(dataFrame):
    """
    Calcule un score Z pour chaque colonne d'un dataFrame tel que :
    z = (nUp - nDown) / sqrt(nUp + nDown)
    où nUp et nDown sont le nombre de cellules strictement positives et negatives d'une colonne.
    """
    return dataFrame.apply(lambda col: getUpDownZScore(*countPosAndNeg(col)))

def getDataFrameMeanAggregate(dataFrame):
    """Renvoie la moyenne des lignes d'un dataframe."""
    return dataFrame.mean()
def getDataFrameNormalZAggregate(dataFrame):
    """
    Applique une normalisation Z pour chaque ligne d'un dataFrame tel que :
    Xnorm = (X - mu) / sigma
    où mu et sigma sont les moyennes et ecart-types respectifs d'une ligne.

    Renvoie la moyenne des lignes normalisees.
    """
    return getDataFrameMeanAggregate(zscore(dataFrame, axis=1))

def getDataFrameAggregate(dataFrame, method):
    """Methodes d'aggregation compatibles avec des donnees de type ratio."""
    algorithms = {
        'mean': getDataFrameMeanAggregate,
        'minStd': getDataFrameMinStdRow,
        'maxStd': getDataFrameMaxStdRow,
        'upDownZ': getDataFrameUpDownZAggregate
    }
    return algorithms[method](dataFrame)