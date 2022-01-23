"""
Main script to run to complete the Tulip project.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

from handle_requests import filterBiocycPathways
from handle_graphs import getWorkingGraph, renameLabelsWithProperty
from handle_genes import loadGeneFiles
from handle_reactions import parseGeneAssociation, setReactionExpressionProperty
from handle_pathways import drawPathwaySubGraphs, drawQuotientGraphs, getAllPathwaysExpression
from handle_heatmap import getHeatmap, convertToDataFrame

nTimestamps = 17
reactionExpressionMethod = 'normalZ' # options: 'mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ'
pathwayExpressionMethod = 'upDownZ' # options: 'mean', 'maxStd', 'minStd', 'upDownZ'

def main(graph):
    
    wg = getWorkingGraph(graph)
    
    # query BioCyc to get pathways and remove nodes without pathways
    pathways = filterBiocycPathways(wg)
    print(f"{len(pathways)} pathways found.")

    loadGeneFiles()
    renameLabelsWithProperty(wg, 'id')
     
    # compute the expression score of reactions
    reactionIdToGenes = parseGeneAssociation(wg)
    setReactionExpressionProperty(wg, reactionIdToGenes, reactionExpressionMethod)    

    # split the pathways into subgraphs
    drawPathwaySubGraphs(wg, pathways)
    pathwaysExpression = getAllPathwaysExpression(wg, pathwayExpressionMethod)
    # draw quotient graphs for each timestamp
    drawQuotientGraphs(wg, pathwaysExpression, timestamps=range(nTimestamps))
    
    # draw the heatmap
    pathwaysExpressionDataFrame = convertToDataFrame(pathwaysExpression, nTimestamps)
    getHeatmap(pathwaysExpressionDataFrame, clusterize=True)
