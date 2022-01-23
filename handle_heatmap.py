"""
This library is dedicated to heatmaps creation.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

from tulip import tlp
from tulipgui import tlpgui
import pandas as pd
from handle_graphs import getColorScale
import scipy.cluster.hierarchy as spch

def convertToDataFrame(pathwaysExpression, expectedSize):
    """
    Converts the pathways and expression values dictionary into a dataFrame.
    
    Parameters
    ----------
    pathwaysExpression : dict associating a pathwayId (key) to its aggregated expression (value)
    expectedSize :
    """
    # fills with NA to call pd.DataFrame.from_records()
    for pathway, expression in pathwaysExpression.items():
        nExpressionValues = len(expression)
        if nExpressionValues < expectedSize:
            expression += [None] * (expectedSize - nExpressionValues)
            pathwaysExpression[pathway] = expression
    dataFrame =  pd.DataFrame.from_records(pathwaysExpression).T
    dataFrame.columns = [f'tp {t+1}' for t in range(expectedSize)]
    return dataFrame.dropna().astype("float64")

def getNodes(graph, nRows, nCols):
    """
    Creates a dictionary of rows*cols nodes.
    The nodes are indexed by their row and column respectively.
    
    Parameters
    ----------
    graph : tlp.Graph
    nRows : int 
    nCols : int
    
    ----------
    Return : dictionary
    """
    nodes = {}
    for i in range(nRows) : 
        nodes[i] = {}
        for j in range(nCols):
            node = graph.addNode()
            nodes[i][j] = node
    return nodes
   
def getGridMap(graph, nodes):
    """
    Places nodes on a graph to create a gridmap.
    
    Parameters
    ----------
    nodes : dict of nodes
    graph : tlp.Graph
    """
    graph['viewShape'].setAllNodeValue(0)
    graph['viewSize'].setAllNodeValue( (0.95, 1, 0) )
    for iRow in nodes.keys():   
        for iCol in nodes[iRow].keys():
            graph['viewLayout'][nodes[iRow][iCol]] = (iCol, -iRow, 0)

def addDoubleProperty(graph, dataFrame, nodes, propertyName) :
    """
    Adds a double property to a graph.
    
    Parameters
    ----------
    gr : tlp.Graph
    dataFrame : pandas.DataFrame with nRows and nCols
    propertyName : str, name of the created property
    nodes : dict of nodes with nRows keys, and nCols keys for each row
    """
    graph.getDoubleProperty(propertyName)
    for iRow in nodes.keys():
        for iCol in nodes[iRow].keys():
            node = nodes[iRow][iCol]
            propertyValue = dataFrame.iloc[iRow, iCol]
            graph[propertyName][node] = propertyValue

def colorNodes(graph, propertyName):
    """
    Colorizes a graph's nodes according to the property values indicated.

    Parameters
    ----------
    graph : tlp.Graph
    propertyName : name of the property used for the nodes color mapping
    """
    params = tlp.getDefaultPluginParameters("Color Mapping", graph)
    params['property'] = graph[propertyName]
    params["type"] = "uniform"
    params["color scale"] = getColorScale()
    graph.applyColorAlgorithm("Color Mapping", params)

def addScale(graph, dataFrame, propertyName):
    """
    Adds a color scale to the heatmap.
    
    Parametres
    ----------
    graph : tlp.Graph
    dataFrame : pandas.DataFrame
    propertyName : name of the property used for the nodes color mapping
    """    
    globalMax = int(dataFrame.max().max() * 10)
    globalMin = int(dataFrame.min().min() * 10)
    cols = dataFrame.shape[1]

    for i in range(globalMin, globalMax+1, 10):        
        node = graph.addNode()
        graph['viewLayout'][node] = (cols+2, (i-40)/10, 0)
        graph[propertyName][node] = i/10
        graph['viewLabel'][node] = str(i/10)
        graph['viewLabelPosition'][node] = tlp.LabelPosition.Left
        graph['viewFontSize'][node] = 14

def addFeature(graph, feature):
    """
    Adds a template node for any feature (header or index) to the heatmap.
    Returns the corresponding node.

    Parameters
    ----------
    graph : tlp.Graph
    feature : str, the name of a header or an index of the heatmap
    """
    node = graph.addNode()
    graph['viewLabel'][node] = feature
    graph['viewFontSize'][node] = 18
    graph['viewColor'][node] = tlp.Color.White
    return node

def addHeader(graph, feature, iCol):
    """
    Adds the name of a dataFrame column to the heatmap.

    Parameters
    ----------
    graph : tlp.Graph
    feature : str, the name of a header for the heatmap
    iCol : int, the number of the corresponding column
    """
    node = addFeature(graph, feature)
    graph['viewLayout'][node] = (iCol, 1, 0)
    graph['viewSize'][node] = ( (0.2 ,0.2 ,0) )
    graph['viewRotation'][node] = 90
    graph['viewLabelPosition'][node] = tlp.LabelPosition.Right

def addIndex(graph, feature, iRow):
    """
    Adds the name of a dataFrame row to the heatmap.

    Parameters
    ----------
    graph : tlp.Graph
    feature : str, the name of an index for the heatmap
    iRow : int, the number of the corresponding row
    """
    node = addFeature(graph, feature)
    graph['viewLayout'][node] = (-1, -iRow, 0)
    graph['viewSize'][node] = ( (0.1 ,0.9 ,0) )
    graph['viewLabelPosition'][node] = tlp.LabelPosition.Left

def addHeadersAndIndexes(graph, dataFrame):
    """
    Adds the names of the rows and the columns of the dataFrame to the heatmap.

    Parameters
    ----------
    graph : tlp.Graph
    df : pandas.DataFrame
    """
    headers, indexes = dataFrame.columns, dataFrame.index
    for iCol in range(len(headers)):
        addHeader(graph, headers[iCol], iCol)
    for iRow in range(len(indexes)):
        addIndex(graph, indexes[iRow], iRow)
    
#def getCorrelationMatrix(dataFrame):
#    """
#    Creates a paired correlation matrix between the rows of the dataframe.
#    
#    ----------
#    Parameters
#    df : pandas.DataFrame
#    col : str, nom de colonne
#    
#    ----------
#    Return : pandas.DataFrame
#    """
#    transposedDataFrame = dataFrame
#    transposedDataFrame.columns = dataFrame.index
#    transposedDataFrame = transposedDataFrame.astype("float64")    
#    return transposedDataFrame.corr()
    
def clusterizeDataFrame(dataFrame):
    """
    Clusterizes a dataFrame using correlation and distance calculations.
    
    ----------
    Paramètres
    dataFrames : pandas.DataFrame
    
    ----------
    Return : pandas.DataFrame
    """
    correlationMatrix = dataFrame.T.corr()
    pdist = spch.distance.pdist(correlationMatrix)
    linkage = spch.linkage(pdist, method='ward')
    idx = spch.fcluster(linkage, 0.5 * pdist.max(), 'distance')
    dataFrame['cluster'] = idx
    clusterizedDataFrame = dataFrame.sort_values(by='cluster')
    return clusterizedDataFrame.drop(columns='cluster', axis=1)


def drawHeatmap(graph, dataFrame, propertyName, clusterize) :
    """
    Displays a heatmap using expression values stored in a dataFrame.
    
    Parametres
    ----------
    graph : tlp.Graph
    dataFrame : pandas.DataFrame
    propertyName : name of the property used for the nodes color mapping
    clusterize : bool, if True the heatmap rows are clusterized
    """
    nodes = getNodes(graph, *dataFrame.shape)
    if clusterize :
        dataFrame = clusterizeDataFrame(dataFrame)
    
    addDoubleProperty(graph, dataFrame, nodes, propertyName)
    getGridMap(graph, nodes)
    addScale(graph, dataFrame, propertyName)
    colorNodes(graph, propertyName)
    addHeadersAndIndexes(graph, dataFrame)
    
def getHeatmap(dataFrame, propertyName="expression", clusterize=True):
    """
    Creates a new graph and displays a heatmap using expression values stored in a dataFrame.
    
    Parametres
    ----------
    dataFrame : pandas.DataFrame
    propertyName : name of the property used for the nodes color mapping
    clusterize : bool, if True the heatmap rows are clusterized
    """
    heatmapGraph = tlp.newGraph()
    heatmapGraph.setName('heatmap')
    drawHeatmap(heatmapGraph, dataFrame, propertyName, clusterize)
    defaultProperties = tlp.DataSet()
    tlpgui.createView("Node Link Diagram view", heatmapGraph, defaultProperties, show=True)
