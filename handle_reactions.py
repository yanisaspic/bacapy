"""
Library to interact with reactions-level data.

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

def parseGeneAssociation(graph):
    """
    Return dictionary with the BioCyc ID of an element as key,
    and the list of genes as value. The elements can be substrates,
    products, or reactions. Only the genes with an expression are kept.

    Parameters
    ----------
    graph: tlp.Graph()

    Returns 
    -------
    biocycIdToGenes: dict
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
    """
    Return expression data of the genes involved in a reaction.
    If biocycId = ID of a substrate, an empty dataframe

    Parameters
    ----------
    biocycId: str
    biocycIdToGenes: dict

    Returns
    -------
    reactionData: pandas.DataFrame
    """
    reactionData = {'level': {}, 'ratio': {}}
    for dataType in ['level', 'ratio']:
        for geneName in biocycIdToGenes[biocycId]:
            geneData = getGeneData(geneName)
            reactionData[dataType][geneName] = geneData[dataType][0]
        reactionData[dataType] = pd.DataFrame.from_dict(reactionData[dataType], orient="index")
    return reactionData

def getReactionExpression(reactionId, reactionIdToGenes, method):
    """
    Aggregate expression values for a BioCyc element.

    Parameters
    ----------
    reactionId: str
        substrate, product, or reaction' ID from BioCyc
    reactionIdToGenes: dict 
        reactionId as key and genes as values
    method: str
        from ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']
    
    Returns
    -------
    pandas.DataFrame
    """
    nodeData = getReactionData(reactionId, reactionIdToGenes)
    if len(nodeData['level']) == 0:    # substract or product
        return []
    if method=='normalZ':
        return getDataFrameNormalZAggregate(nodeData['level']).tolist()
    return getDataFrameAggregate(nodeData['ratio'], method).tolist()

def setReactionExpressionProperty(graph, reactionIdToGenes, method):
    """
    Add or Update the 'expression' property of graph
    this property contain Expression' values after aggregation of
    BioCyc elements (reaction, substrate, or product).

    Parameters
    ----------
    graph: tlp.Graph
    reactionIdToComponents: dict
        reactionId as key and genes as values
    method: str
        from ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']
    """
    graph.getDoubleVectorProperty('expression')
    for node in graph.getNodes():
        reactionId = graph['id'][node]
        graph['expression'][node] = getReactionExpression(reactionId, reactionIdToGenes, method)
