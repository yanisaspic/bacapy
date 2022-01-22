"""
Unique script a executer pour repondre a la problematique.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

#from handle_requests
from handle_graphs import renameLabelsWithProperty, filterNodesWithIds, getNeighborNodes, getIdsFromNodes, deleteNodes
from handle_genes import loadDataFiles
from handle_reactions import parseGeneAssociation, setReactionExpressionProperty

expressionMethod = 'upDownZ' # options: 'mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ'

def main(graph):
    loadDataFiles()

    # utilisation des ids sur le graphe
    renameLabelsWithProperty(graph, 'id')

    # calcul de l'expression des reactions:
    reactionIdToGenes = parseGeneAssociation(graph)
    setReactionExpressionProperty(graph, reactionIdToGenes, expressionMethod)
