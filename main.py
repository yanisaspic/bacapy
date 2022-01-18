"""
Unique script a executer pour repondre a la problematique.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

from handle_files import loadDataFiles, parseGeneAssociation
from aggregate_expression import setExpressionProperty

expressionMethod = 'normalZ' # options: 'mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ'

def main(graph):    
    loadDataFiles()
    nodeIdToGenes = parseGeneAssociation(graph)
    setExpressionProperty(graph, nodeIdToGenes, expressionMethod)
    # expression = graph.getDoubleVectorProperty('expression')
    # test = []
    # for n in graph.getNodes():
    #     test.append(n)
    # node936 = test[936]
    # print(expression[node936])
