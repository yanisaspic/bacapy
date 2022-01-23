"""
This library is dedicated to miscellaneous Tulip graph functions.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

from tulip import tlp
import pandas as pd

def getColorScale():
    colorScale = tlp.ColorScale([])
    colorScale.setColorAtPos(-1.0, tlp.Color.Blue)
    colorScale.setColorAtPos(0.0, tlp.Color.White)
    colorScale.setColorAtPos(1.0, tlp.Color.Red)
    return colorScale

def getRootGraph():
    for graph in tlp.getRootGraphs():
        return graph
def newSubGraph(graph, subGraphName):
    """
    Create a new sub-graph
    If a graph of the same name already exist, it is delete beforehand.
    
    Parametres
    ----------
    graph : tlp.Graph
    subGraphName : str

    ----------
    Return : tlp.Graph
    """
    sg = graph.getSubGraph(subGraphName)
    if sg != None:
        graph.delAllSubGraphs(sg)
    return graph.addCloneSubGraph(subGraphName)
def getWorkingGraph(graph, workingGraphName='working graph', originalGraphName='original graph'):
    wg, og = newSubGraph(graph, workingGraphName), newSubGraph(graph, originalGraphName)
    return wg
def renameSubGraph(graph, subGraphName, subGraphNewName):
    """
    Rename a sub-graph
    
    Parameters
    ----------
    graph : tlp.Graph
    subGraphName : str
    subGraphNewName : str
    """
    sg = graph.getSubGraph(subGraphName)
    sg.setAttribute('name', subGraphNewName)

def renameLabelsWithProperty(graph, propertyName):
    """
    Rename the labels of a graph' nodes using one of the graph' property.
    
    Parameters
    ----------
    graph : tlp.Graph
    propertyName : str
    """
    currentLabels = graph.getStringProperty('viewLabel')
    newLabels = graph.getStringProperty(propertyName)
    for n in graph.getNodes():
        currentLabels[n] = newLabels[n]

def getIdsFromNodes(graph, nodes):
    """
    Return an array of nodes' IDs from the graph.

    Parameters
    ----------
    graph : tlp.Graph
    nodes : array, tlp.Node

    ----------
    Return : array, str
    """
    nodeIds = []
    ids = graph.getStringProperty('id')
    for n in nodes:
        nodeIds.append(ids[n])
    return nodeIds
def getNeighborNodes(graph, node):
    """
    For a given node, return the list of all its direct neighbors
    
    Parameters
    ----------
    graph : tlp.Graph
    nodes : array, tlp.Node

    ----------
    Return : array, tlp.Node
    """
    return tlp.reachableNodes(graph, node, 1, direction=tlp.UNDIRECTED)
def deleteNodes(graph, nodes):
    """
    Delete nodes from graph
    
    Parameters
    ----------
    graph : tlp.Graph
    nodes : array, tlp.Node
    """
    for n in nodes:
        graph.delNode(n)

def splitNodesWithIds(graph, nodeIds):
    """
    Separate nodes depending on where they belong in an array of nodes' ids
    
    Parameters
    ----------
    graph : tlp.Graph
    nodeIds : list, str

    ----------
    Return : tuple, list, tlp.Node
    """
    includedNodes, excludedNodes = [], []
    ids = graph.getStringProperty('id')
    for n in graph.getNodes():
        if ids[n] in nodeIds:
            includedNodes.append(n)
        else:
            excludedNodes.append(n)
    return includedNodes, excludedNodes
def filterNodesWithIds(graph, nodeIds, excluded=False):
    """
    Return nodes depending on their belonging in an array of nodes' IDs.
    
    Parameters
    ----------
    graph : tlp.Graph
    nodeIds : list, str
    excluded : boolean, (False) keep nodes in the list

    ----------
    Return : tuple, list, tlp.Node
    """
    return splitNodesWithIds(graph, nodeIds)[excluded]

def getQuotientGraph(graph, quotientGraphName='quotient graph'):
    """
    Create Quotient Graph from a graph with sub-graph
    
    Parameters
    ----------
    graph : tlp.Graph
    quotientGraphName = str
    
    ----------
    Return : tlp.Graph
    """
    params = tlp.getDefaultPluginParameters("Quotient Clustering", graph)
    params["use name of subgraph"] = True
    #params["layout quotient graph(s)"] = True
    graph.applyAlgorithm("Quotient Clustering", params)
    quotientGraphDefaultName = 'quotient of ' + graph.getAttribute('name')
    renameSubGraph(getRootGraph(), quotientGraphDefaultName, quotientGraphName)
    return getRootGraph().getSubGraph(quotientGraphName)
    
def getExpression(graph):
    """
    Return a dataFrame with a graph' nodes' expressions levels.
    
    Parameters
    ----------
    graph : tlp.Graph

    ----------
    Return : pandas.DataFrame
    """
    expression = graph['expression']
    expressionList = []
    for n in graph.getNodes():
        expressionList.append(expression[n])
    return pd.DataFrame.from_records(expressionList)
