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
from handle_graphs import getWorkingGraph, renameLabelsWithProperty
from handle_genes import loadDataFiles
from handle_reactions import parseGeneAssociation, setReactionExpressionProperty
from handle_pathways import drawPathways
from handle_requests import mainBioCyc


expressionMethod = 'upDownZ' # options: 'mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ'

def main(graph):
    loadDataFiles()
    wg = getWorkingGraph(graph)
    renameLabelsWithProperty(wg, 'id')
    pathways = mainBioCyc(graph)    # takes several minutes  
    
    # calcul de l'expression des reactions:
#    reactionIdToGenes = parseGeneAssociation(graph)
#    setReactionExpressionProperty(graph, reactionIdToGenes, expressionMethod)
    
#
    
    quotientGraphNames = []
    for i in range(1, 18):
        quotientGraphNames.append(f"tp {i}") 
    drawPathways(wg, pathways, quotientGraphNames)
