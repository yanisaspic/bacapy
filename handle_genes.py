"""
Librairie pour interagir avec les structures de donnees d'expression genique.

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
    Load 3 Files : the locus, their expression levels, and their differential expression
    """
    # genes, levels and ratios are directly called by other functions : use global
    global genes, levels, ratios
    genes = pd.read_csv(genesFilename, sep=';')
    levels = pd.read_csv(levelsFilename, sep=';')
    ratios = pd.read_csv(ratiosFilename, sep=';')
    genes = genes[genes['locus'].isin(levels['locus'])]

def getGeneLocus(geneName):
    """Return the locus corresponding to a gene."""
    return genes.loc[genes['gene name']==geneName, 'locus'].values[0]
def getLocusLevel(locus):
    """Return the expression level of a locus."""
    return levels.loc[levels['locus']==locus].drop('locus', axis=1).to_numpy()
def getLocusRatio(locus):
    """Return the differential expression of a locus."""
    return ratios.loc[ratios['locus']==locus].drop('locus', axis=1).to_numpy()
def getLocusData(locus):
    """Return the expression level and the differential expression of a locus."""
    return {'level': getLocusLevel(locus), 'ratio': getLocusRatio(locus)}
def getGeneData(geneName):
    """Return the expression level and the differential expression of a gene."""
    return getLocusData( getGeneLocus(geneName) )

def isGeneWithData(geneName):
    """ 
    Check if the expression of a gene was measured.
    
    -------
    Return : Boolean
    """
    try:
        getGeneData(geneName)
        return True
    except IndexError:
        return False