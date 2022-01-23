"""
Library to interact with our gene expression data files.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

import pandas as pd

genesFilename = "mapGeneLocus.csv"
levelsFilename = "ecoliK12_levels.csv"
ratiosFilename = "ecoliK12_ratio.csv"

def loadGeneFiles():
    """
    Load 3 Files : the locus, their expression levels,
    and their differential expression.
    """
    # genes, levels and ratios are directly called
    # by other functions : use global
    global genes, levels, ratios
    genes = pd.read_csv(genesFilename, sep = ';')
    levels = pd.read_csv(levelsFilename, sep = ';')
    ratios = pd.read_csv(ratiosFilename, sep = ';')
    genes = genes[genes['locus'].isin(levels['locus'])]

def getGeneLocus(geneName):
    """
    Return the locus corresponding to a gene.
    
    Parameters
    ----------
    geneName: str
        the gene name
        
    Returns
    -------
    locus: str
        the locus corresponding to the gene 
    """
    locus = genes.loc[genes['gene name'] == geneName, 'locus'].values[0]
    return locus

def getLocusLevel(locus):
    """
    Return the expression level of a locus.
    
    Parameters
    ----------
    locus: str
        the locus code
        
    Returns
    -------
    expr: numpy array
        array of expression levels
    """
    expr = levels.loc[levels['locus'] == locus].drop('locus', axis = 1).to_numpy()
    return expr

def getLocusRatio(locus):
    """
    Return the differential expression of a locus.
    
    Parameters
    ----------
    locus: str
        the locus code
        
    Returns
    -------
    ratio: numpy array
        array of differential expressions
    """
    ratio = ratios.loc[ratios['locus'] == locus].drop('locus', axis = 1).to_numpy()
    return ratio

def getLocusData(locus):
    """
    Return the expression level and the differential
    expression of a locus.
    
    Parameters
    ----------
    locus: str
        the locus code
        
    Returns
    -------
    lvl_ratio: dict
        dict of expressions data
    """
    lvl_ratio = {'level': getLocusLevel(locus),
                 'ratio': getLocusRatio(locus)}
    return lvl_ratio

def getGeneData(geneName):
    """
    Return the expression level and the differential
    expression of a gene.
    
    Parameters
    ----------
    geneName: str
        the gene name 
    
    Returns
    -------
    data: dict
        dict of expression data
    """
    data = getLocusData(getGeneLocus(geneName))
    return data

def isGeneWithData(geneName):
    """ 
    Check if the expression of a gene was measured.
    
    Parameters
    ----------
    geneName: str
        the gene name 
    
    Return
    ------
    boolean
    """
    try:
        getGeneData(geneName)
        return True
    except IndexError:
        return False
