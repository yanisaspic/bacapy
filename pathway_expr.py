from tulipgui import tlpgui
from tulip import tlp
import pandas as pd
from aggregate_expression import *
from heatmap import *




def get_expr_pathway(graph, method = "mean"):
    """
    Créer un dictionnaire rassemblant les données d'expression en liste 
    de liste des réactions contenus dans chaque pathway
    
    Paramètres
    ----------
    graph : tlp.Graph
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']

    ----------
    Return : dictionnaire python
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
    Rempli les pathways d'une liste de "NA" si ils ne contiennent aucune expression
    
    Paramètres
    ----------
    path_expr : tlp.Graph
    tp : int, nombre de timestamp
    
    ----------
    Return : pandas.Dataframe
    """
    None_list = ["NA"]*tp
    for i in path_expr : 
        if len(path_expr[i]) == 0:
             path_expr[i] = None_list
    return path_expr
             
def get_aggregated_expr_pathway(graph, tp = 17, method="mean"):
    """
    Créer un dataframe qui aggrègent les données d'expression de réaction
    par pathway
    
    Paramètres
    ----------
    graph : tlp.Graph
    tp : int, nombre de timestamp
    method : string parmi ['mean', 'maxStd', 'minStd', 'upDownZ', 'normalZ']
    
    ----------
    Return : pandas.Dataframe
    """
    
    path_expr = get_expr_pathway(graph, method)
    path_expr = fill_na(path_expr, tp)
    
    pathway_expression = pd.DataFrame.from_dict(path_expr).T
    tp = "tp"
    col = []
    [col.append( tp+str(i) ) for i in range(0,17)]
    pathway_expression.columns = col     
    return(pathway_expression)

             
def main(graph):
    wg=graph.getSubGraph("Workingraph")
    quotient=graph.getSubGraph("quotient of Workingraph")
    
    pathway_expression = get_aggregated_expr_pathway(wg, method = "mean")
    
    roots = tlp.getRootGraphs()
    heat = None
    has_a_view = False
    for gr in roots :
       if gr.getName() == "heatmap" :
           heat=gr
           heat.clear()
           has_a_view = True

    if heat == None : 
        heat = tlp.newGraph()
        heat.setName("heatmap")

    df = cluster_row(pathway_expression)

    heatmap(heat, df, "expr", True) 
    
    if not has_a_view : 
        prop = tlp.DataSet()
        #view = tlpgui.getViewsOfGraph(heat)    
        heat = tlpgui.createView("Node Link Diagram view",heat,prop,True)  
