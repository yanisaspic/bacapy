import pandas as pd
import numpy as np
from tulip import tlp

import scipy.cluster.hierarchy as spc
from handle_files import *


def get_nodes(gr, lignes=2, colonnes=2, nodes = {}):
    """
    Créer autant de nodes que nécessaire pour la heatmap
    
    Paramètres
    ----------
    gr : tlp.Graph
    lignes : int 
    colonnes : int
    nodes : python dictionnaire 
    
    """
    for i in range(lignes) : 
        nodes[i] = {}
        for j in range(colonnes):
            node = gr.addNode()
            nodes[i][j] = node
    
    return nodes
   
def grid_map(nodes, lignes, colonnes, graph  ):
    """
    Ordonne les noeuds pour créer une heatmap
    
    Paramètres
    ----------
    nodes : dictionnaire python contenant la liste des noeuds
    lignes : int 
    colonnes : int
    graph : tlp.Graph
    
    """
    graph['viewShape'].setAllNodeValue(0)
    graph['viewSize'].setAllNodeValue( (0.9 ,1 ,0) )
    for i in range(colonnes):   
        for j in range(lignes):

            # On place les noeud de la première colonne plus à gauche : nom des lignes
            #if i == 0 :
             #   graph['viewLayout'][nodes[j][i]] = (i-3,-j,0)
            #else :
            graph['viewLayout'][nodes[j][i]] = (i,-j,0)

def add_property(graph,df, nodes, property_name="expr") :
    """
    Ajoute une propriété (type int/double) 
    permettant de colorer les cases de la heatmap
    
    
    Paramètres
    ----------
    gr : tlp.Graph
    df : pandas.DataFrame
    property_name : str, nom de la propriété à créer
    nodes : python dictionnaire 
    
    """
    prop = graph.getDoubleProperty(property_name)
    
    ligne = 0
    col = 0
    while col < len(nodes[0]) :
        for row in range( len(nodes) ):
            value = df.iloc[row,col]               
            graph[property_name][ nodes[row][col] ] = value
        col += 1


def color_heatmap(graph, prop, nodes):
    """
    Colore la heatmap en fonction de la propriété donnée en paramètres
    
    Paramètres
    ----------
    graph : tlp.Graph
    prop : nom de la propriété du graph à mapper sur les noeuds 
    nodes : python dictionnaire 
    
    """
    # dict qui récup les paramètres
    params = tlp.getDefaultPluginParameters("Color Mapping", graph)
    
    # On map la propriété 'prop' dans à la case ('property')
    params['property'] = graph[prop]
    
    colorScale = tlp.ColorScale([])    
    colorScale.setColorAtPos(0.0, tlp.Color.Blue)
    colorScale.setColorAtPos(0.5, tlp.Color.White)
    colorScale.setColorAtPos(1.0, tlp.Color.Red)
    params["color scale"] = colorScale
    
    graph.applyColorAlgorithm("Color Mapping", params)


def add_scale(graph, df, prop):
    """
    Ajoute une échelle de couleur
    
    Paramètres
    ----------
    graph : tlp.Graph
    df : pandas.DataFrame
    prop : nom de la propriété du graph à mapper sur les noeuds 

    """
    col = df.shape[1]
    
    global_max = int(max(df.max().tolist()[1:]))
    global_min = int(min(df.min().tolist()[1:]))
    
    for i in range(global_min, global_max+1, 1):
        #print(i)        
        node = graph.addNode()
        graph['viewLayout'][node] = (col+2, i-4, 0)
        graph[prop][node] = int(i)
        graph['viewLabel'][node] = str(i)
        # position du label : gauche
        graph['viewLabelPosition'][node] = 3
        
def add_features(graph, df, feature) :
    """
    Ajoute les nomes des colonnes ou des lignes à la heatmap
    
    Paramètres
    ----------
    graph : tlp.Graph
    df : pandas.DataFrame
    feature : str ["header","index"] 

    """
    
    if feature == "header" :
        ft = df.columns
        
    if feature == "index" :
        ft = df.index
        
    count = 0
    for col in ft :
        
        node = graph.addNode()
    #    graph['viewFont'][node] = "C:/Users/Antonin Colajanni/Desktop/M2/R_BOURQUI/bacapy/OpenSans-Regular.ttf"
        #graph['viewFont'][node] =  
        graph['viewLabel'][node] = col
        graph['viewFontSize'][node] = 18
        #graph['viewLabelPosition'][node] = 4

        #graph['viewLabelPosition'][node] = 1
        graph['viewColor'][node] = tlp.Color.White
        
                  
        if feature == "header" :
            graph['viewLayout'][node] = (count, 1, 0)
            graph['viewRotation'][node] = 90
                    
        if feature == "index" :
            graph['viewLayout'][node] = (-3, -count, 0)
            graph['viewSize'][node] = ( (0.9 ,4 ,0) )

        count +=1


def heatmap(graph, df, property_name) :
    """
    Affiche une heatmap à partir d'un pandas.DataFrame
    
    Paramètres
    ----------
    graph : tlp.Graph
    df : pandas.DataFrame
    prop : nom de la propriété du graph à mapper sur les noeuds 

    """
    
    lignes = df.shape[0]
    col = df.shape[1]
    nodes = get_nodes(graph, lignes, col)
    add_property(graph, df, nodes,  property_name)
    grid_map(nodes, lignes, col, graph)
    add_scale(graph, df, property_name)
    color_heatmap(graph, property_name, nodes)
    add_features(graph, df, "header")
    add_features(graph, df, "index")
    
#def cluster_row():
    

### Test
def main(graph):
    
    #graph.clear()
    


    genes, levels, ratios = loadDataFiles()

    ratio = ratios.iloc[30:130, :]
    lignes = ratio.shape[0]
    col = ratio.shape[1]

    
    """
    genes = ratio['locus']  
    rownames = ratio.locus.to_list()
    ratio = ratio.drop(axis=1, columns= "locus")
    ratio.index = rownames  
    
    
    transposed = ratio.T
    transposed.columns = genes
    transposed = transposed.drop(axis=0, index = 'locus')
    transposed = transposed.astype("float64")    
    corr_matrix = transposed.corr()
    
    pdist = spc.distance.pdist(corr_matrix)
    
    linkage = spc.linkage(pdist, method='ward' )
    idx = spc.fcluster(linkage, 0.5 * pdist.max(), 'distance')
    
    ratio['clust'] = idx
    ratio_clust = ratio.sort_values(by=['clust'], )
    ratio = ratio_clust.drop(columns = 'clust', axis = 1)
    

    heatmap(graph, ratio,  "aaa")
    """


      
