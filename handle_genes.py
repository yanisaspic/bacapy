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
    """Charge 3 fichiers : les locus, leur niveau d'expression et leur expression differentielle."""
    # genes, levels et ratios directement appeles par d'autres fonctions : utiliser global
    global genes, levels, ratios
    genes = pd.read_csv(genesFilename, sep=';')
    levels = pd.read_csv(levelsFilename, sep=';')
    ratios = pd.read_csv(ratiosFilename, sep=';')
    genes = genes[genes['locus'].isin(levels['locus'])]

def getGeneLocus(geneName):
    """Renvoie le locus correspondant a un gene."""
    return genes.loc[genes['gene name']==geneName, 'locus'].values[0]
def getLocusLevel(locus):
    """Renvoie le niveau d'expression d'un locus."""
    return levels.loc[levels['locus']==locus].drop('locus', axis=1).to_numpy()
def getLocusRatio(locus):
    """Renvoie l'expression differentielle d'un locus."""
    return ratios.loc[ratios['locus']==locus].drop('locus', axis=1).to_numpy()
def getLocusData(locus):
    """Renvoie le niveau d'expression et l'expression differentielle d'un locus."""
    return {'level': getLocusLevel(locus), 'ratio': getLocusRatio(locus)}
def getGeneData(geneName):
    """Renvoie le niveau d'expression et l'expression differentielle d'un gene."""
    return getLocusData( getGeneLocus(geneName) )

def isGeneWithData(geneName):
    """Verifie si l'expression d'un gene a ete mesuree."""
    try:
        getGeneData(geneName)
        return True
    except IndexError:
        return False