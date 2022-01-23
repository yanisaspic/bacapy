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
    """
    Create a new color scale.
    
    Returns
    -------
    colorScale: tlp.ColorScale
        a color scale
    """
    colorScale = tlp.ColorScale([])
    colorScale.setColorAtPos(0.0, tlp.Color.Blue)
    colorScale.setColorAtPos(0.05, tlp.Color.Azure)
    colorScale.setColorAtPos(0.33, tlp.Color.White)
    colorScale.setColorAtPos(0.95, tlp.Color.OrangeRed)
    colorScale.setColorAtPos(1.0, tlp.Color.Red)
    return colorScale

def getRootGraph():
    """
    Return the root graph.
    
    Returns
    -------
    graph: tlp.Graph
    """
    for graph in tlp.getRootGraphs():
        return graph
    
def newSubGraph(graph, subGraphName):
    """
    Create a new sub-graph. If a graph of the same name
    already exist, it is delete beforehand.
    
    Parameters
    ----------
    graph: tlp.Graph
    subGraphName: str

    Returns
    -------
    tlp.Graph
    """
    sg = graph.getSubGraph(subGraphName)
    if sg != None:
        graph.delAllSubGraphs(sg)
    return graph.addCloneSubGraph(subGraphName)


def getWorkingGraph(graph, workingGraphName = 'working graph',
                    originalGraphName = 'original graph'):
    """
    Get the working graph.
    
    Parameters
    ----------
    graph:  tlp.Graph
    workingGraphName: str
    originalGraphName: str
    
    Returns
    -------
    wg: tlp.Graph
    """
    og = newSubGraph(graph, originalGraphName)
    wg = newSubGraph(graph, workingGraphName)
    return wg


def renameSubGraph(graph, subGraphName, subGraphNewName):
    """
    Rename a sub-graph.
    
    Parameters
    ----------
    graph: tlp.Graph
    subGraphName: str
    subGraphNewName: str
    """
    sg = graph.getSubGraph(subGraphName)
    sg.setAttribute('name', subGraphNewName)

def renameLabelsWithProperty(graph, propertyName):
    """
    Rename the labels of a graph' nodes using one of the graph' property.
    
    Parameters
    ----------
    graph: tlp.Graph
    propertyName: str
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
    graph: tlp.Graph
    nodes: list
        the list of Tlp.Node 

    Returns
    -------
    nodeIds: list
        the list of nodes id (str)

    """
    nodeIds = []
    ids = graph.getStringProperty('id')
    for n in nodes:
        nodeIds.append(ids[n])
    return nodeIds

def getNeighborNodes(graph, node):
    """
    For a given node, return the list of all its 
    direct neighbors.
    
    Parameters
    ----------
    graph: tlp.Graph
    node: tlp.Node

    Returns
    -------
    nn: list
        the list of neighor nodes (Tlp.Node) 
    
    """
    nn = tlp.reachableNodes(graph, node, 1, direction = tlp.UNDIRECTED)
    return nn

def deleteNodes(graph, nodes):
    """
    Delete nodes from graph.
    
    Parameters
    ----------
    graph : tlp.Graph
    nodes: list
        the list of Tlp.Node 
    """
    for n in nodes:
        graph.delNode(n)

def splitNodesWithIds(graph, nodeIds):
    """
    Separate nodes depending on where they belong in an
    array of nodes' ids.
    
    Parameters
    ----------
    graph : tlp.Graph
    nodeIds : list, str

    Returns
    -------
    includedNodes: list
    excludedNodes: list
    """
    includedNodes, excludedNodes = [], []
    ids = graph.getStringProperty('id')
    for n in graph.getNodes():
        if ids[n] in nodeIds:
            includedNodes.append(n)
        else:
            excludedNodes.append(n)
    return includedNodes, excludedNodes

def filterNodesWithIds(graph, nodeIds, excluded = False):
    """
    Return nodes depending on their belonging in an
    array of nodes' IDs.
    
    Parameters
    ----------
    graph: tlp.Graph
    nodeIds: list
        the list of nodes id (str)
    excluded: boolean
        if false, keep nodes in the list

    Returns
    -------
    nodeIdList: list
        list of excluded or included nodes id (str)
    """
    nodesIdList = splitNodesWithIds(graph, nodeIds)[excluded]

def getQuotientGraph(graph, quotientGraphName = 'quotient graph'):
    """
    Create Quotient Graph from a graph with sub-graph.
    
    Parameters
    ----------
    graph : tlp.Graph
    quotientGraphName: str
    
    Returns
    -------
    tlp.Graph
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

    Returns
    -------
    pandas.DataFrame
    """
    expression = graph['expression']
    expressionList = []
    for n in graph.getNodes():
        expressionList.append(expression[n])
    return pd.DataFrame.from_records(expressionList)
