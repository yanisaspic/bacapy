"""
Librairie pour interagir avec les elements metaboliques macro (pathways)

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

from handle_graphs import filterNodesWithIds, getNeighborNodes, getQuotientGraph

def getPathwayNodes(graph, pathwayId, pathwayIdsToReactions):
    """Renvoie la liste des noeuds impliques dans une voie metabolique (reaction, substrat et produit).
    
    Parametres
    ----------
    graph : objet tlp.Graph()
    pathwayId : string, l'id BioCyc d'un pathway
    pathwayIdsToReactions : dictionnaire associant un pathwayId (cle) a ses reactions (valeur)
    """
    reactionNodes = filterNodesWithIds(graph, pathwayIdsToReactions[pathwayId])
    pathwayNodes = set(reactionNodes)
    for n in reactionNodes:
        pathwayNodes.update(getNeighborNodes(graph, n))
    return list(pathwayNodes)
def drawPathways(graph, pathwayIdsToReactions):
    """Dessine des sous-graphes correspondant a une voie metabolique unique."""
    for pathwayId in pathwayIdsToReactions:
        pathwayNodes = getPathwayNodes(graph, pathwayId, pathwayIdsToReactions)
        graph.inducedSubGraph(pathwayNodes, name=pathwayId)
    return getQuotientGraph(graph)