"""
Librairie pour interagir avec les structures de donnees d'expression genique.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

import pandas as pd

genesFilename = "mapGeneLocus.csv"
levelsFilename = "ecoliK12_levels.csv"
ratiosFilename = "ecoliK12_ratio.csv"

def loadDataFiles():
    # genes, levels et ratios directement appeles par d'autres fonctions : utiliser global
    global genes, levels, ratios
    genes = pd.read_csv(genesFilename, sep=';')
    levels = pd.read_csv(levelsFilename, sep=';')
    ratios = pd.read_csv(ratiosFilename, sep=';')
    genes = genes[genes['locus'].isin(levels['locus'])]

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

def getNodeData(nodeId, nodeIdToGenes):
    nodeData = {'level': {}, 'ratio': {}}
    for dataType in ['level', 'ratio']:
        for geneName in nodeIdToGenes[nodeId]:
            geneData = getGeneData(geneName)
            nodeData[dataType][geneName] = geneData[dataType][0]
        nodeData[dataType] = pd.DataFrame.from_dict(nodeData[dataType], orient="index")
    return nodeData