
from tulip import tlp
import pandas as pd
from aggregate_expression import *
from heatmap import *

def main(graph):
    wg=graph.getSubGraph("Workingraph")
    quotient=graph.getSubGraph("quotient of Workingraph")
    
    
    
    path_expr = {}
    for subgraph in wg.getSubGraphs():
        pathway = subgraph.getName() 
        
        expr_list = []
        for node in subgraph.getNodes():
            expr = subgraph["expression"][node]   
            if expr != []  :
                expr_list.append(expr)
        
        expr_df = pd.DataFrame.from_records(expr_list)
        aggregated_expr = getDataFrameRow(expr_df, "mean")

        path_expr[pathway] = aggregated_expr.to_list()
    
    
    None_list = [0]*17
    for i in path_expr : 
        if len(path_expr[i]) == 0:
             path_expr[i] = None_list
    
    pathway_expression = pd.DataFrame.from_dict(path_expr).T
    tp = "tp"
    col = []
    for i in range(1,18):
        col.append( tp+str(i) )
    pathway_expression.columns = col
    

    
    heat = graph.getSubGraph("heatmap")
    if heat == None :
        heat = graph.addSubGraph("heatmap")
    else : 
        heat.clear()

    #heat = tlp.loadGraph("C:/Users/Antonin Colajanni/Desktop/M2/R_BOURQUI/bacapy/Empty_graph.tlpx")
    #print(heat)
    
    heatmap(heat, pathway_expression, "expr")
    
    
    #print(graph.getRoot())
 
    """
    for path in path_expr :
        print("\n Pathway :",path, "nombre de gènes :", len(path_expr[path]) )
        for i in path_expr[path]:
            print(i)
    """
