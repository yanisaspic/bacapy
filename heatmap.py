import pandas as pd
from scipy.stats import zscore
import math
from statistics import mean
from tulip import tlp
#from tulipgui import tlpgui
from scipy.stats import zscore
import scipy.cluster.hierarchy as spc

"""
V A R I A B L E S
"""
genesFilename = "mapGeneLocus.csv"
levelsFilename = "ecoliK12_levels.csv"
ratiosFilename = "ecoliK12_ratio.csv"

"""
F U N C T I O N S
"""
def loadDataFiles(genesFilename="mapGeneLocus.csv",levelsFilename="ecoliK12_levels.csv",ratiosFilename="ecoliK12_ratio.csv" ):
    genes = pd.read_csv(genesFilename, sep=';')
    levels = pd.read_csv(levelsFilename, sep=';')    
    ratios = pd.read_csv(ratiosFilename, sep=';')
    return genes, levels, ratios

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
            if i == 0 :
                graph['viewLayout'][nodes[j][i]] = (i-3,-j,0)
            else :
                graph['viewLayout'][nodes[j][i]] = (i,-j,0)

def add_property(graph,df, nodes, property_name="expr") :
    """
    à refaire ?
    
    
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

            if is_integer(value):
                
                graph["expr"][ nodes[row][col] ] = value
            else : 
                graph["viewLabel"][ nodes[row][col] ] = value

        col += 1

def is_integer(value):
    """
    Renvoit True si value est un integer
    
    Paramètres
    ----------
    value : int or str
    
    """
    try :
        int(value)
        return True
    except ValueError :
        return False

def color_heatmap(graph, prop, nodes):
    """
    Colore la heatmap
    
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
    colorScale.setColorAtPos(1.0, tlp.Color.Brown)
    params["color scale"] = colorScale
    
    graph.applyColorAlgorithm("Color Mapping", params)
    
    lignes = len(nodes)
    colonnes = len(nodes[0])
    for i in range(colonnes):   
        for j in range(lignes):
            if i == 0 :
                graph['viewColor'][nodes[j][i]] = (255,255,255)
    
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
        
def add_header(graph, df) :
    """
    Ajoute les colonnes de la heatmap
    
    Paramètres
    ----------
    graph : tlp.Graph
    df : pandas.DataFrame

    """

    colonne = 0
    for col in df.columns[1:] :
        
        node = graph.addNode()
  
        graph['viewLayout'][node] = (colonne+1, 1, 0)
        graph['viewLabel'][node] = col
        graph['viewFontSize'][node] = 16
        graph['viewRotation'][node] = 90
        graph['viewLabelPosition'][node] = 4

        #graph['viewLabelPosition'][node] = 1
        graph['viewColor'][node] = (255,255,255)
        
        colonne +=1

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
    
    grid_map(nodes, lignes, col, graph)
    add_scale(graph, df, property_name)
    add_property(graph, df, nodes)
    color_heatmap(graph, property_name, nodes)
    add_header(graph, df)
    
#def cluster_row():
    
    
    
def main(graph):
    
    graph.clear()
    
    
    genes, levels, ratios = loadDataFiles()

    ratio = ratios.iloc[30:130, :]
    
    lignes = ratio.shape[0]
    col = ratio.shape[1]
    
    
    
    genes = ratio['locus']    
    transposed = ratio.T
    transposed.columns = genes
    transposed = transposed.drop(axis=0, index = 'locus')
    transposed = transposed.astype("float64")    
    corr_matrix = transposed.corr()
    
    pdist = spc.distance.pdist(corr_matrix)
    linkage = spc.linkage(pdist, method='ward')
    idx = spc.fcluster(linkage, 0.5 * pdist.max(), 'distance')
    ratio['clust'] = idx
    ratio_clust = ratio.sort_values(by=['clust'])
    ratio = ratio_clust.drop(columns = 'clust', axis = 1)
    
    heatmap(graph, ratio,  "expr")

        
