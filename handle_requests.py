"""
Librairie pour interagir avec la base de donnees BioCyc.

@ ASLOUDJ Yanis
@ COLAJANNI Antonin
@ DUGUE Berenice
@ JACQUES Patrick
@ SAUVESTRE Clement
@ SIMON Arnaud
"""

import requests
from xml.etree import ElementTree as ET
import handle_graphs as hg
from handle_reactions import getReactionIds
import time

def requestBiocyc(ID):
    """
    Performs a request for an object (reaction, pathway ...)
    with its ID in BioCyc. 
    
    Parameters:
        ID (string): ID of the requested object
    
    Returns:
        doc (xml): The BioCyc XML text of our requested object 
        
    """
    URL = "https://websvc.biocyc.org/getxml?ECOLI:" + ID
    response = requests.get(URL, timeout = 5)
    if response.status_code == 200:
        doc = ET.fromstring(response.text)
    # in case of a temporary ban
    elif response.status_code == 429:
        time.sleep(60)
        response = requests.get(URL, timeout = 5)
        doc = ET.fromstring(response.text)
    return doc

def getPathwayID(pathways, reaction):
    """
    Adds the metabolic pathway associated with the reaction
    to the list of pathways. 
    
    Parameters:
        pathways   (list): Unique list of pathways associated
                           to our reactions in BioCyc
        reaction (string): ID of the reaction for which we search
                           its pathways in BioCyc
    
    """
    doc = requestBiocyc(reaction)
    for e in doc.findall(".//in-pathway/Pathway"):
        ID = e.attrib['frameid']
        if ID not in pathways:
            pathways.append(ID)
   
def getPathways(reactions):
    """
    Gets the list of pathways associated to reactions in
    Biocyc and the list of reactions to delete not found 
    in BioCyc.
    
    Parameters:
        reactions (list): List of reaction IDs retrieved 
                          from our Tulip graph
                          
    Returns:
        pathways  (list): Unique list of pathways associated
                          to our reactions in BioCyc
        to_delete (list): List of reactions to delete in our
                          Tulip graph

    """
    pathways = []
    to_delete = []
    for r in reactions:
        try:
            getPathwayID(pathways, r)      
        except :
            to_delete.append(r)
    return pathways, to_delete

def getReactionsID(pathway):
    """
    Retrieves in BioCyc all reactions associated to a pathway.
    
    Parameters:
        pathway (string): Pathway id
        
    Returns:
        reactions (list): List of associated reactions
        
    """
    reactions = []
    doc = requestBiocyc(pathway)
    for e in doc.findall(".//reaction-list/Reaction"):
        reactions.append(e.attrib['frameid'])
    return reactions

def getReactions(pathways):
    """
    Associate for each identified metabolic pathway the 
    reactions that are associated with it.
    
    Parameters:
        pathways (list): List of identified pathways
    
    Returns:
        data (dict): String of pathway id as key and its 
                     associated reaction list as values
        
    """
    data = {}
    for p in pathways:
        try :
            reactions_list = getReactionsID(data, p)
            data[p] = reactions_list
        except :
            pass
    return data

def filterBiocycReactions(graph, notBR):
    """
    Deletes all reactions not found on BioCyc. Substrates are 
    also deleted if they are not used by found reactions.
    
    Parameters:
        graph (tlp.Graph): Ecoli K12 substrates - reactions graph 
        notBR (list): List of reaction IDs not found on BioCyc
    
    """
    biocycReactionNodes = hg.filterNodesWithIds(graph, notBR, excluded = True)
    nodeIdsToKeep = set(hg.getIdsFromNodes(biocycReactionNodes))
    for n in biocycReactionNodes:
        neighborNodes = hg.getNeighborNodes(graph, n)
        nodeIdsToKeep.update(hg.getIdsFromNodes(neighborNodes))
    hg.deleteNodes(graph, hg.filterNodesWithIds(graph, nodeIdsToKeep, excluded = True))

def mainBiocyc(graph):
    reactions = getReactionIds(graph)
    pathways, to_del = getPathways(reactions)
    data = getReactions(pathways)
    filterBiocycReactions(graph, to_del)
    return data
    