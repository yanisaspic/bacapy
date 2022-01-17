import pandas as pd
import math
from scipy.stats import zscore
import tulip as tlp

"""
V A R I A B L E S
"""
genesFilename = "mapGeneLocus.csv"
levelsFilename = "ecoliK12_levels.csv"
ratiosFilename = "ecoliK12_ratio.csv"

"""
F U N C T I O N S
"""
def loadDataFiles():
    genes = pd.read_csv(genesFilename, sep=';')
    levels = pd.read_csv(levelsFilename, sep=';')
    ratios = pd.read_csv(ratiosFilename, sep=';')
    return genes[genes['locus'].isin(levels['locus'])], levels, ratios

def getGeneLocus(geneName):
    return genes.loc[genes['gene name']==geneName, 'locus'].values[0]
def getLocusLevel(locus):
    return levels.loc[levels['locus']==locus].drop('locus', axis=1).to_numpy()
def getLocusRatio(locus):
    return ratios.loc[ratios['locus']==locus].drop('locus', axis=1).to_numpy()
def getLocusData(locus):
    return {'level': getLocusLevel(locus), 'ratio': getLocusRatio(locus)}
def getGeneData(geneName):
    return getLocusData( getGeneLocus(geneName) )

def getNodeData(nodeId, nodeIdToGenes):
    nodeData = {'level': {}, 'ratio': {}}
    for dataType in ['level', 'ratio']:
        for geneName in nodeIdToGenes[nodeId]:
            geneData = getGeneData(geneName)
            nodeData[dataType][geneName] = geneData[dataType][0]
        nodeData[dataType] = pd.DataFrame.from_dict(nodeData[dataType], orient="index")
    return nodeData

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

def isGeneWithData(geneName):
    try:
        getGeneData(geneName)
        return True
    except IndexError:
        return False
def parseGeneAssociation(graph):
    """
    Renvoie un dictionnaire associant l'id BioCyc d'un element (cle) a des genes (valeur).
    Seuls les genes dont l'expression est connue sont conserves.

    Paramètres
    ----------
    graph : objet tlp.graph() 
    """
    nodeIdToGenes = {}
    for node in graph.getNodes():
        nodeId = graph['id'][node]
        nodeIdToGenes[nodeId] = []
        for mot in graph['geneAssociation'][node].split(' or '):
            if isGeneWithData(mot[2:-2]):
                nodeIdToGenes[nodeId].append(mot[2:-2])
    return nodeIdToGenes
    
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
    
"""
M A I N
"""
genes, levels, ratios = loadDataFiles()
expressionMethod = 'upDownZ' # options: 'mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ'

def main(graph):
    print(genes)
    nodeIdToGenes = parseGeneAssociation(graph)
    setExpressionProperty(graph, nodeIdToGenes, expressionMethod)
