"""
This library is dedicated to handling pathways data.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

import handle_graphs as hg
from tulip import tlp
from aggregate_data import getDataFrameAggregate

def getPathwayNodes(graph, pathwayId, pathwayIdsToReactions):
    """
    Return list of nodes involved in a pathway (reaction, substrate, and product).
    
    Parameters
    ----------
    graph : tlp.Graph
    pathwayId : str, BioCyc ID of a pathway
    pathwayIdsToReactions : python Dictionary, pathwayId as key and reactions as values
    """
    reactionNodes = hg.filterNodesWithIds(graph, pathwayIdsToReactions[pathwayId])
    pathwayNodes = set(reactionNodes)
    for n in reactionNodes:
        pathwayNodes.update(hg.getNeighborNodes(graph, n))
    return list(pathwayNodes)
def drawPathwaySubGraphs(graph, pathwayIdsToReactions):
    """
    Draw subgraph for each pathway
    
    Parameters
    ----------
    graph : tlp.Graph
    pathwayIdsToReactions : python Dictionary, pathwayId as key and reactions as values
    """
    for pathwayId in pathwayIdsToReactions:
        pathwayNodes = getPathwayNodes(graph, pathwayId, pathwayIdsToReactions)
        graph.inducedSubGraph(pathwayNodes, name=pathwayId)

def getOnePathwayExpression(pathwaySubGraph, method):
    """
    Return aggregated expression of a pathway from expression values of its reactions
    
    Parameters
    ----------
    pathwaySubGraph : tlp.Graph
    method : string, from ['mean', 'maxStd', 'minStd', 'upDownZ']

    ----------
    return : list
    """
    reactionsExpression = hg.getExpression(pathwaySubGraph)
    # erreur avec minStd et maxStd quand les expressions de reactions n'ont pas ete mesurees
    if len(reactionsExpression)==0:
        return []
    return getDataFrameAggregate(reactionsExpression, method).to_list()

def getAllPathwaysExpression(graph, method):
    """
    Return Dictionary with pathway as keys and list of aggregated expressions as values
    (computed from expression of involved reactions).
    The expression is computed from sub-graphs each corresponding to independent pathways
    
    Parameters
    ----------
    graph : tlp.Graph
    method : str, from ['mean', 'maxStd', 'minStd', 'upDownZ']

    ----------
    Return : python Dictionary
    """
    pathwayIdToExpression = {}
    for pathwaySubGraph in graph.getSubGraphs():
        pathwayName = pathwaySubGraph.getName()
        pathwayExpression = getOnePathwayExpression(pathwaySubGraph, method)
        pathwayIdToExpression[pathwayName] = pathwayExpression
    return pathwayIdToExpression

def setPathwayExpressionProperty(quotientGraph, pathwaysExpression):
    """
    Add or Update the Expression Property of graph.
    With the aggregated expression values of BioCyc Elements (reaction, substrate, or product)
   
    Parameters
    ----------
    quotientGraph : tlp.Graph
    pathwaysExpression : python Dictionary, pathwayId as key and aggregated expression as value
    """
    quotientGraph.getDoubleVectorProperty('expression')
    for n in quotientGraph.getNodes():
        pathwayLabel = quotientGraph['viewLabel'][n]
        quotientGraph['expression'][n] = pathwaysExpression[pathwayLabel]

def customizeGraph(graph, timestamp):
    """
    Adapt size and color of graph' nodes depending on expression level at timestamp

    Parameters
    ----------
    graph : tlp.Graph
    timestamp : int
    """
    for n in graph.getNodes():
        pathwayExpression = graph['expression'][n]
        # on n'affiche pas les pathways dont l'expression n'a pas ete mesuree :
        if len(pathwayExpression) == 0:
            graph['viewSize'][n] = (0, 0, 0)
        else:
            tpExpression = pathwayExpression[timestamp]
            size = abs(tpExpression)*400
            graph['viewSize'][n] = (size, size, 0)
            graph['viewColor'][n] = hg.getColorScale().getColorAtPos(tpExpression)
        graph['viewShape'][n] = tlp.NodeShape.Circle

def drawQuotientGraphs(graph, pathwaysExpression, timestamps):
    """
    Draw Quotient Graph of pathways for each timestamp

    Parameters
    ----------
    graph : tlp.Graph
    timestamps : list of int
    """
    for t in timestamps:
        quotientGraphName = f'tp{t+1} quotient graph'
        qg = hg.getQuotientGraph(graph, quotientGraphName)
        setPathwayExpressionProperty(qg, pathwaysExpression)
        customizeGraph(qg, t)
