"""
Librairie pour aggreger les donnees d'expression genique.

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
from handle_files import getNodeData

def getDataFrameRowByStd(dataFrame, ascend):
    selectedRow = dataFrame.std(axis=1).sort_values(ascending=ascend).index[0]
    return dataFrame.loc[selectedRow]
def getDataFrameMinStdRow(dataFrame):
    return getDataFrameRowByStd(dataFrame, True)
def getDataFrameMaxStdRow(dataFrame):
    return getDataFrameRowByStd(dataFrame, False)

def countPosAndNeg(series):
    return len(series[series>0]), len(series[series<0])
def getUpDownZScore(nUp, nDown):
    # ZeroDivisionError si nUp==0 et nDown==0 => utilise max(1).
    return (nUp-nDown) / math.sqrt(max(1, nUp+nDown))
def getDataFrameUpDownZRow(dataFrame):
    """
    Calcule un score Z pour chaque colonne d'un dataFrame tel que :
    z = (nUp - nDown) / sqrt(nUp + nDown)
    où nUp et nDown sont le nombre de cellules strictement positives et negatives d'une colonne.
    
    @ Marini, F., Ludt, A., Linke, J. et al. 
    GeneTonic: an R/Bioconductor package for streamlining the interpretation of RNA-seq data. 
    BMC Bioinformatics 22, 610 (2021). 
    https://doi.org/10.1186/s12859-021-04461-5
    """
    return dataFrame.apply(lambda col: getUpDownZScore(*countPosAndNeg(col)))

def getDataFrameMeanRow(dataFrame):
    return dataFrame.mean()
def getDataFrameNormalZRow(dataFrame):
    """
    Applique une normalisation Z pour chaque ligne d'un dataFrame tel que :
    Xnorm = (X - mu) / sigma
    où mu et sigma sont les moyennes et ecart-types respectifs d'une ligne.

    Renvoie la moyenne des lignes normalisees.
    """
    return getDataFrameMeanRow(zscore(dataFrame, axis=1))

def getDataFrameRow(dataFrame, method):
    """Methodes d'aggregation compatibles avec des donnees de type ratio."""
    algorithms = {
        'mean': getDataFrameMeanRow,
        'minStd': getDataFrameMinStdRow,
        'maxStd': getDataFrameMaxStdRow,
        'upDownZ': getDataFrameUpDownZRow
    }
    return algorithms[method](dataFrame)

def getNodeExpression(nodeId, nodeIdToGenes, method):
    """
    Applique un algorithme d'aggregation des valeurs d'expression pour un element BioCyc.

    Parametres
    ----------
    nodeId : string, l'id BioCyc d'un element (substrat, produit ou reaction)
    nodeIdToGenes : dictionnaire associant un nodeId (cle) a des genes (valeur)
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']
    """
    nodeData = getNodeData(nodeId, nodeIdToGenes)
    if len(nodeData['level'])==0:    # substract or product
        return []
    if method=='normalZ':
        return getDataFrameNormalZRow(nodeData['level']).tolist()
    return getDataFrameRow(nodeData['ratio'], method).tolist()

def setExpressionProperty(graph, nodeIdToGenes, method):
    """
    Ajoute ou met a jour la propriete 'expression' du graph.
    Elle contient les valeurs d'expression apres aggregation des elements BioCyc.

    Parametres
    ----------
    graph : objet tlp.Graph()
    nodeIdToGenes : dictionnaire associant l'id BioCyc d'un element (cle) a des genes (valeur)
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']
    """
    graph.getDoubleVectorProperty('expression')
    for node in graph.getNodes():
        nodeId = graph['id'][node]
        graph['expression'][node] = getNodeExpression(nodeId, nodeIdToGenes, method)