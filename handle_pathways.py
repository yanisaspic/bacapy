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
    """Renvoie la liste des noeuds impliques dans une voie metabolique (reaction, substrat et produit).
    
    Parametres
    ----------
    graph : objet tlp.Graph()
    pathwayId : string, l'id BioCyc d'un pathway
    pathwayIdsToReactions : dictionnaire associant un pathwayId (cle) a ses reactions (valeur)
    """
    reactionNodes = hg.filterNodesWithIds(graph, pathwayIdsToReactions[pathwayId])
    pathwayNodes = set(reactionNodes)
    for n in reactionNodes:
        pathwayNodes.update(hg.getNeighborNodes(graph, n))
    return list(pathwayNodes)
def drawPathwaySubGraphs(graph, pathwayIdsToReactions):
    """Dessine des sous-graphes correspondant a une voie metabolique unique."""
    for pathwayId in pathwayIdsToReactions:
        pathwayNodes = getPathwayNodes(graph, pathwayId, pathwayIdsToReactions)
        graph.inducedSubGraph(pathwayNodes, name=pathwayId)

def getOnePathwayExpression(pathwaySubGraph, method):
    """
    Renvoie l'expression agregee d'une voie metabolique d'apres les valeurs
    d'expression des reactions associees.
    
    Parametres
    ----------
    pathwaySubGraph : tlp.Graph
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ']
    """
    reactionsExpression = hg.getExpression(pathwaySubGraph)
    # erreur avec minStd et maxStd quand les expressions de reactions n'ont pas ete mesurees
    if len(reactionsExpression)==0:
        return []
    return getDataFrameAggregate(reactionsExpression, method).to_list()

def getAllPathwaysExpression(graph, method):
    """
    Renvoie un dictionnaire associant une voie metabolique (cle) a une liste de 
    donnees d'expression agregees (valeur) calculee depuis l'expression des reactions associees.
    L'expression est calculee depuis des sous-graphes correspondant chacun a
    une voie metabolique unique. 
    
    Parametres
    ----------
    graph : tlp.Graph
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ']

    ----------
    Return : dictionnaire
    """
    pathwayIdToExpression = {}
    for pathwaySubGraph in graph.getSubGraphs():
        pathwayName = pathwaySubGraph.getName()
        pathwayExpression = getOnePathwayExpression(pathwaySubGraph, method)
        pathwayIdToExpression[pathwayName] = pathwayExpression
    return pathwayIdToExpression

def setPathwayExpressionProperty(quotientGraph, pathwaysExpression):
    """Ajoute ou met a jour la propriete expression du graph.
    Elle contient les valeurs d'expression apres aggregation des elements BioCyc (reaction, substrat ou produit).

    Parametres
    ----------
    quotientGraph : objet tlp.Graph
    pathwaysExpression : dictionnaire associant un pathwayId (cle) a son expression agregee (valeur)
    timestamp : estampille correspondant au graph quotient"""
    quotientGraph.getDoubleVectorProperty('expression')
    for n in quotientGraph.getNodes():
        pathwayLabel = quotientGraph['viewLabel'][n]
        quotientGraph['expression'][n] = pathwaysExpression[pathwayLabel]

def customizeGraph(graph, timestamp):
    """Adapte la taille et la couleur des noeuds d'un graphe selon le niveau d'expression au timestamp donne."""
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
    """Dessine un graphe quotient de pathways pour chaque timestamp indique."""
    for t in timestamps:
        quotientGraphName = f'tp{t+1} quotient graph'
        qg = hg.getQuotientGraph(graph, quotientGraphName)
        setPathwayExpressionProperty(qg, pathwaysExpression)
        customizeGraph(qg, t)
