"""
Librairie pour interagir avec les elements metaboliques micro (reaction, substrat ou produit)

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

import pandas as pd
from handle_genes import isGeneWithData, getGeneData
from aggregate_data import getDataFrameAggregate, getDataFrameNormalZAggregate

def getReactionIds(graph):
    """Renvoie l'ensemble des ids BioCyc correspondant a une reaction."""
    ids=graph.getStringProperty('id')
    isReaction=graph.getBooleanProperty('reaction')
    reactionIds=[]
    for n in graph.getNodes():
        if isReaction:
            reactionIds.append(ids[n])
    return reactionIds

def parseGeneAssociation(graph):
    """
    Renvoie un dictionnaire associant l'id BioCyc d'un element (cle) a des genes (valeur).
    Les elements sont des substrats, des produits ou des reactions.
    Seuls les genes dont l'expression a ete mesuree sont conserves.

    Parametres
    ----------
    graph : objet tlp.graph() 
    """
    biocycIdToGenes = {}
    for node in graph.getNodes():
        biocycId = graph['id'][node]
        biocycIdToGenes[biocycId] = []
        for mot in graph['geneAssociation'][node].split(' or '):
            if isGeneWithData(mot[2:-2]):
                biocycIdToGenes[biocycId].append(mot[2:-2])
    return biocycIdToGenes

def getReactionData(biocycId, biocycIdToGenes):
    """Renvoie les mesures d'expression des genes impliques dans une reaction.
    Si biocycId correspond a un substrat, les dataframes renvoyes sont vides."""
    reactionData = {'level': {}, 'ratio': {}}
    for dataType in ['level', 'ratio']:
        for geneName in biocycIdToGenes[biocycId]:
            geneData = getGeneData(geneName)
            reactionData[dataType][geneName] = geneData[dataType][0]
        reactionData[dataType] = pd.DataFrame.from_dict(reactionData[dataType], orient="index")
    return reactionData

def getReactionExpression(reactionId, reactionIdToGenes, method):
    """
    Applique un algorithme d'aggregation des valeurs d'expression pour un element BioCyc.

    Parametres
    ----------
    reactionId : string, l'id BioCyc d'un element (substrat, produit ou reaction)
    reactionIdToComponents : dictionnaire associant un reactionId (cle) a ses genes (valeur)
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']
    """
    nodeData = getReactionData(reactionId, reactionIdToGenes)
    if len(nodeData['level'])==0:    # substract or product
        return []
    if method=='normalZ':
        return getDataFrameNormalZAggregate(nodeData['level']).tolist()
    return getDataFrameAggregate(nodeData['ratio'], method).tolist()

def setReactionExpressionProperty(graph, reactionIdToGenes, method):
    """
    Ajoute ou met a jour la propriete 'expression' du graph.
    Elle contient les valeurs d'expression apres aggregation des elements BioCyc (reaction, substrat ou produit).

    Parametres
    ----------
    graph : objet tlp.Graph()
    reactionIdToComponents : dictionnaire associant un reactionId (cle) a ses genes (valeur)
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']
    """
    graph.getDoubleVectorProperty('expression')
    for node in graph.getNodes():
        reactionId = graph['id'][node]
        graph['expression'][node] = getReactionExpression(reactionId, reactionIdToGenes, method)