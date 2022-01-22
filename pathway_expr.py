from tulipgui import tlpgui
from tulip import tlp
import pandas as pd
from aggregate_expression import *
from heatmap import *


def get_expr_pathway(graph, method = "mean"):
    """
    Créer un dictionnaire rassemblant les données d'expression 
    des réactions contenus dans chaque pathway
    
    Paramètres
    ----------
    graph : tlp.Graph
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']
    
    """
    path_expr = {}
    for subgraph in graph.getSubGraphs():
        pathway = subgraph.getName() 
        
        expr_list = []
        for node in subgraph.getNodes():
            expr = subgraph["expression"][node]   
            if expr != []  :
                expr_list.append(expr)
        
        expr_df = pd.DataFrame.from_records(expr_list)
        aggregated_expr = getDataFrameRow(expr_df, method)
        path_expr[pathway] = aggregated_expr.to_list()
        
    return path_expr

def fill_na(path_expr, tp = 17):
    """
    Rempli les pathways d'une liste de 0 si ils contiennent aucune expression
    
    Paramètres
    ----------
    graph : tlp.Graph
    tp : int, nombre de timestamp
    
    """
    None_list = [0]*tp
    for i in path_expr : 
        if len(path_expr[i]) == 0:
             path_expr[i] = None_list
    return path_expr
             
def get_aggregated_expr_pathway(graph, tp = 17, method="mean"):
    """
    Créer un dictionnaire rassemblant les données d'expression 
    des réactions contenus dans chaque pathway et les aggrègent
    
    Paramètres
    ----------
    graph : tlp.Graph
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']
    
    """
    
    path_expr = get_expr_pathway(graph, method)
    path_expr = fill_na(path_expr, tp)
    
    
    pathway_expression = pd.DataFrame.from_dict(path_expr).T
    tp = "tp"
    col = []
    #for i in range(1,18):
    [col.append( tp+str(i) ) for i in range(1,18)]
    pathway_expression.columns = col     
    return(pathway_expression)
         
             
def main(graph):
    wg=graph.getSubGraph("Workingraph")
    quotient=graph.getSubGraph("quotient of Workingraph")
    
    pathway_expression = get_aggregated_expr_pathway(wg)
    
    heat = graph.getSubGraph("heatmap")
    if heat == None :
        heat = graph.addSubGraph("heatmap")
    else : 
        heat.clear()
        
    heat = tlp.newGraph()
    heat.setName("heatmap")
    
   
    
    print(tlpgui.getViewsOfGraph(heat))

    heatmap(heat, pathway_expression, "expr")  
    prop = tlp.DataSet()
    view = tlpgui.getViewsOfGraph(heat)    
    print(tlpgui.getViewsOfGraph(graph))
    heat = tlpgui.createView("Node Link Diagram view",heat,prop,True)  
    
    
    
    
    #heat.delLocalProperty("viewFont")
    #print(font)
    #print(pathway_expression)
    
    """
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
    

    """
    

    

    """
    for path in path_expr :
        print("\n Pathway :",path, "nombre de gènes :", len(path_expr[path]) )
        for i in path_expr[path]:
            print(i)
    """
