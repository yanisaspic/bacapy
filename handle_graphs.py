"""
Librairie pour la gestion des graphes Tulip.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

def newSubGraph(graph, subGraphName):
    """Cree un nouveau sous-graphe. 
    Si un graphe du meme nom existe deja, il est prealablement supprime."""
    sg = graph.getSubGraph(subGraphName)
    if sg != None:
        graph.delAllSubGraphs(sg)
    return graph.addCloneSubGraph(subGraphName)

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
    neighborNodes = []
    for nn in graph.getInOutNodes(node):
        neighborNodes.append(nn)
    return neighborNodes
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