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
    """Cree un nouveau sous-graphe. 
    Si un graphe du meme nom existe deja, il est prealablement supprime."""
    sg = graph.getSubGraph(subGraphName)
    if sg != None:
        graph.delAllSubGraphs(sg)
    return graph.addCloneSubGraph(subGraphName)
def getWorkingGraph(graph, workingGraphName='working graph', originalGraphName='original graph'):
    wg, og = newSubGraph(graph, workingGraphName), newSubGraph(graph, originalGraphName)
    return wg
def renameSubGraph(graph, subGraphName, subGraphNewName):
    """Renomme un sous-graphe."""
    sg = graph.getSubGraph(subGraphName)
    sg.setAttribute('name', subGraphNewName)

def renameLabelsWithProperty(graph, propertyName):
    """Renomme les labels des noeuds d'un graphe d'apres une propriete du graphe."""
    currentLabels = graph.getStringProperty('viewLabel')
    newLabels = graph.getStringProperty(propertyName)
    for n in graph.getNodes():
        currentLabels[n] = newLabels[n]

def getIdsFromNodes(graph, nodes):
    """Renvoie une liste d'ids correspondant aux noeuds du graphe indiques."""
    nodeIds = []
    ids = graph.getStringProperty('id')
    for n in nodes:
        nodeIds.append(ids[n])
    return nodeIds
def getNeighborNodes(graph, node):
    """Pour un noeud donne, renvoie la liste de ses voisins directs sur le graphe."""
    return tlp.reachableNodes(graph, node, 1, direction=tlp.UNDIRECTED)
def deleteNodes(graph, nodes):
    """Supprime les noeuds indiques du graphe."""
    for n in nodes:
        graph.delNode(n)

def splitNodesWithIds(graph, nodeIds):
    """Separe les noeuds selon leur appartenance a une liste d'ids de noeuds."""
    includedNodes, excludedNodes = [], []
    ids = graph.getStringProperty('id')
    for n in graph.getNodes():
        if ids[n] in nodeIds:
            includedNodes.append(n)
        else:
            excludedNodes.append(n)
    return includedNodes, excludedNodes
def filterNodesWithIds(graph, nodeIds, excluded=False):
    """Renvoie des noeuds selon leur appartenance a une liste d'ids de noeuds.
    
    Parametres
    ----------
    graph : objet tlp.Graph()
    nodeIds : list de strings, les ids des noeuds
    excluded : bool, les noeuds conserves sont ceux inclus (False) ou exclus (True) de la liste
    """
    # False==0 et True==1
    return splitNodesWithIds(graph, nodeIds)[excluded]

def getQuotientGraph(graph, quotientGraphName='quotient graph'):
    """Construit un graphe quotient depuis un graphe avec des sous-graphes."""
    params = tlp.getDefaultPluginParameters("Quotient Clustering", graph)
    params["use name of subgraph"] = True
    #params["layout quotient graph(s)"] = True
    graph.applyAlgorithm("Quotient Clustering", params)
    quotientGraphDefaultName = 'quotient of ' + graph.getAttribute('name')
    renameSubGraph(getRootGraph(), quotientGraphDefaultName, quotientGraphName)
    return getRootGraph().getSubGraph(quotientGraphName)
    
def getExpression(graph):
    """
    Renvoie un dataFrame contenant les expressions des noeuds d'un graphe.
    
    Parametres
    ----------
    graph : objet tlp.Graph
    """
    expression = graph['expression']
    expressionList = []
    for n in graph.getNodes():
        expressionList.append(expression[n])
    return pd.DataFrame.from_records(expressionList)
