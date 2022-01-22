import pandas as pd
import numpy as np
from tulip import tlp
from tulipgui import tlpgui

import scipy.cluster.hierarchy as spc
from handle_files import *


def get_nodes(gr, lignes=2, colonnes=2, nodes = {}):
    """
    Créer autant de nodes que nécessaire pour la heatmap
    Les stocke dans un dictionnaire indicé par leur 
    numéro de ligne et colonne
    
    Paramètres
    ----------
    gr : tlp.Graph
    lignes : int 
    colonnes : int
    nodes : dictionnaire python
    
    ----------
    Return : dictionnaire python
    
    """
    for i in range(lignes) : 
        nodes[i] = {}
        for j in range(colonnes):
            node = gr.addNode()
            nodes[i][j] = node
    
    return nodes
   
def grid_map(nodes, lignes, colonnes, graph):
    """
    Ordonne les noeuds pour créer une grille
    
    Paramètres
    ----------
    nodes : dictionnaire python contenant la liste des noeuds
    lignes : int 
    colonnes : int
    graph : tlp.Graph
    
    """
    graph['viewShape'].setAllNodeValue(0)
    graph['viewSize'].setAllNodeValue( (0.95 ,1 ,0) )
    for i in range(colonnes):   
        for j in range(lignes):

            graph['viewLayout'][nodes[j][i]] = (i,-j,0)

def add_property(graph,df, nodes, property_name="expr") :
    """
    Ajoute une propriété (type double) 
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
    Colore la heatmap en fonction de la propriété donnée en paramètre
    
    Paramètres
    ----------
    graph : tlp.Graph
    prop : nom de la propriété du graph à mapper sur les noeuds 
    nodes : dictionnaire python
    
    """
    params = tlp.getDefaultPluginParameters("Color Mapping", graph)

    params['property'] = graph[prop]
    #params["override minimum value"] = True
    #params["override maximum value"] = True
    params["type"] = "uniform"
    
    colorScale = tlp.ColorScale([])    
    colorScale.setColorAtPos(0.0, tlp.Color.Blue)
    colorScale.setColorAtPos(0.5, tlp.Color.White)
    colorScale.setColorAtPos(1.0, tlp.Color.Red)
    params["color scale"] = colorScale
    
    graph.applyColorAlgorithm("Color Mapping", params)


def add_scale(graph, df, prop):
    """
    Ajoute une échelle de couleur à la heatmap
    
    Paramètres
    ----------
    graph : tlp.Graph
    df : pandas.DataFrame
    prop : nom de la propriété du graph à mapper sur les noeuds 

    """
    col = df.shape[1]
    
    global_max = int(max(df.max().tolist())*10)
    global_min = int(min(df.min().tolist())*10)
        
    for i in range(global_min, global_max+1, 10):        
        node = graph.addNode()
        graph['viewLayout'][node] = (col+2, (i-40)/10, 0)
        graph[prop][node] = i/10
        graph['viewLabel'][node] = str(i/10)
        graph['viewLabelPosition'][node] = tlp.LabelPosition.Left
        graph['viewFontSize'][node] = 14
        
def add_features(graph, df, feature) :
    """
    Ajoute les noms des colonnes ou des lignes à la heatmap
    
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
        graph['viewLabel'][node] = col
        graph['viewFontSize'][node] = 18
        graph['viewShape'][node] = tlp.NodeShape.CubeOutlinedTransparent
        graph['viewColor'][node] = tlp.Color.White
                  
        if feature == "header" :
            graph['viewLayout'][node] = (count, 1, 0)
            graph['viewSize'][node] = ( (0.2 ,0.2 ,0) )
            graph['viewRotation'][node] = 90
            graph['viewLabelPosition'][node] = tlp.LabelPosition.Right
                    
        if feature == "index" :
            graph['viewLayout'][node] = (-0.5, -count, 0)
            graph['viewSize'][node] = ( (0.1 ,0.9 ,0) )
            graph['viewLabelPosition'][node] = tlp.LabelPosition.Left

        count +=1

def remove_NA(df,col="tp1") : 
    """
    Enlève les lignes d'un dataframe contenant "NA"
    
    ----------
    Paramètres
    df : pandas.DataFrame
    col : str, nom de colonne
    
    ----------
    Return : pandas.DataFrame
    """
    return df.loc[df[col] != "NA"]
    
def get_corr_mat(df):
    """
    Créer une matrice de corrélation par 
    paire entre les lignes du dataframe
    
    ----------
    Paramètres
    df : pandas.DataFrame
    col : str, nom de colonne
    
    ----------
    Return : pandas.DataFrame
    """
    df = remove_NA(df)
    transposed = df.T
    transposed.columns = df.index
    transposed = transposed.astype("float64")    
    return transposed.corr()
    
def cluster_row(df):
    """
    Ordonne les colonne d'un Dataframe en fonction de la 
    corrélation entre les lignes
    
    
    ----------
    Paramètres
    df : pandas.DataFrame
    
    ----------
    Return : pandas.DataFrame
    """

    df = remove_NA(df)
    corr_matrix = get_corr_mat(df)
    pdist = spc.distance.pdist(corr_matrix)
    linkage = spc.linkage(pdist, method='ward' )
    idx = spc.fcluster(linkage, 0.5 * pdist.max(), 'distance')
    df['clust'] = idx
    df_clust = df.sort_values(by=['clust'], )
    df = df_clust.drop(columns = 'clust', axis = 1)
    return df 


def heatmap(graph, df, property_name="expr", cluster = True) :
    """
    Affiche une heatmap à partir d'un pandas.DataFrame
    
    Paramètres
    ----------
    graph : tlp.Graph
    df : pandas.DataFrame
    prop : nom de la propriété du graph à mapper sur les noeuds 
    cluster : Booleen, True : Ordonne les lignes de la heatmap

    """
    
    lignes = df.shape[0]
    col = df.shape[1]
    nodes = get_nodes(graph, lignes, col)
    
    if cluster :
        df = cluster_row(df)
    
    add_property(graph, df, nodes,  property_name)
    grid_map(nodes, lignes, col, graph)
    add_scale(graph, df, property_name)
    color_heatmap(graph, property_name, nodes)
    add_features(graph, df, "header")
    add_features(graph, df, "index")
