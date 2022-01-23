"""
Library to aggregate quantitative data.

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
    """
    Return a dataframe row depending on its standard deviation' value.

    Parameters
    ----------
    dataFrame : pandas.DataFrame
    ascend : Boolean

    ----------
    Returns : pandas.DataFrame
    """
    # Prevent error message if empty dataframe
    selectedRow = dataFrame.std(axis=1).sort_values(ascending=ascend).index[0]
    return dataFrame.loc[selectedRow]
def getDataFrameMinStdRow(dataFrame):
    """
    Return Dataframe row with smallest standard deviation
    
    Parameters
    ----------
    dataFrame : pandas.DataFrame

    ----------
    Returns : pandas.DataFrame
    """
    return getDataFrameRowByStd(dataFrame, True)
def getDataFrameMaxStdRow(dataFrame):
    """
    Return Dataframe row with largest standard deviation
    
    Parameters
    ----------
    dataFrame : pandas.DataFrame

    ----------
    Returns : pandas.DataFrame
    """
    return getDataFrameRowByStd(dataFrame, False)

def countPosAndNeg(series):
    """
    Return the number of positive and negative values in a pandas Series.
    
    Parameters
    ----------
    series : pandas.Series

    ----------
    Returns : tuple
    """
    return len(series[series>0]), len(series[series<0])
def getUpDownZScore(nUp, nDown):
    """
    Return a Z-score centered around 0.
    
    @ Marini, F., Ludt, A., Linke, J. et al. 
    GeneTonic: an R/Bioconductor package for streamlining the interpretation of RNA-seq data. 
    BMC Bioinformatics 22, 610 (2021). 
    https://doi.org/10.1186/s12859-021-04461-5
    
    Parameters
    ----------
    nUp : double
    nDown : double

    ----------
    Returns : double
    """
    # ZeroDivisionError si nUp==0 et nDown==0 => utilise max(1).
    return (nUp-nDown) / math.sqrt(max(1, nUp+nDown))
def getDataFrameUpDownZAggregate(dataFrame):
    """
    Compute Z-score for each column of a Dataframe such as :
    z = (nUp - nDown) / sqrt(nUp + nDown)
    Where nUp and nDown are the numbers of cells strictly positives and negatives of the column.
    
    Parameters
    ----------
    dataFrame : pandas.DataFrame

    ----------
    Returns : pandas.DataFrame
    """
    return dataFrame.apply(lambda col: getUpDownZScore(*countPosAndNeg(col)))

def getDataFrameMeanAggregate(dataFrame):
    """
    Return average of dataFrame' rows.

    Parameters
    ----------
    dataFrame : pandas.DataFrame

    ----------
    Returns : Series or DataFrame (if level specified)
    """
    return dataFrame.mean()
def getDataFrameNormalZAggregate(dataFrame):
    """
    Apply Z-score normalisation for each row of dataFrame such as :
    Xnorm = (X - mu) / sigma
    Where mu and sigma are the mean and std of the row.

    Parameters
    ----------
    dataFrame : pandas.DataFrame

    ----------
    Returns : Series or DataFrame (if level specified)
    """
    return getDataFrameMeanAggregate(zscore(dataFrame, axis=1))

def getDataFrameAggregate(dataFrame, method):
    """
    Aggregation methods compatibles with ratio data

    Parameters
    ----------
    dataFrame : pandas.DataFrame
    """
    algorithms = {
        'mean': getDataFrameMeanAggregate,
        'minStd': getDataFrameMinStdRow,
        'maxStd': getDataFrameMaxStdRow,
        'upDownZ': getDataFrameUpDownZAggregate
    }
    return algorithms[method](dataFrame)
