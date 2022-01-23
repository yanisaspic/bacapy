"""
This library is dedicated to BioCyc interactions.

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
import time

def requestBiocyc(ID):
    """
    Performs a request for an object (reaction, pathway ...)
    with its ID in BioCyc. 
    
    Parameters
    ----------
    ID : str, requested object' ID
    
    ----------
    Return : doc (xml), BioCyc XML text of our requested object 
        
    """
    # monitor requests 
    URL = "https://websvc.biocyc.org/getxml?ECOLI:" + ID
    response = requests.get(URL, timeout = 5)
    print(f"{round(time.time() - start, 2)}: {ID} ({response.status_code})")
    if response.status_code == 200:
        doc = ET.fromstring(response.text)
    # in case of a temporary ban
    elif response.status_code == 429:
        time.sleep(60)
        doc = requestBiocyc(ID)
    return doc

def getPathwaysID(pathways, reaction):
    """
    Adds the metabolic pathways associated with the reaction
    to the list of pathways. 
    
    Parameters
    ----------
    pathways   (list): Unique list of pathways associated
                        to our reactions in BioCyc
    reaction : str, ID of the reaction for which we search
                        its pathways in BioCyc
    
    """
    doc = requestBiocyc(reaction)
    for e in doc.findall(".//in-pathway/Pathway"):
        ID = e.attrib['frameid']
        if ID not in pathways:
            pathways.append(ID)
    return pathways
   
def getPathways(reactions):
    """
    Gets the list of pathways associated to reactions in
    Biocyc and the list of reactions to delete not found 
    in BioCyc.
    
    Parameters
    ----------
    reactions (list): List of reaction IDs retrieved 
                      from our Tulip graph
    
    ----------             
    Return :
        pathways  (list): Unique list of pathways associated
                          to our reactions in BioCyc
        to_delete (list): List of reactions to delete in our
                          Tulip graph

    """
    pathways = []
    to_delete = []
    for r in reactions:
        try:
            pathways = getPathwaysID(pathways, r)
        except :
            to_delete.append(r)
    return pathways, to_delete

def getReactionIdsFromPathway(pathway):
    """
    Retrieves in BioCyc all reactions associated to a pathway.
    
    Parameters
    ----------
    pathway (string): Pathway id
        
    ----------
    Return : reactions (list): List of associated reactions ids
        
    """
    print(f"{time.time() - start}: {pathway}")
    reactions = []
    doc = requestBiocyc(pathway)
    for e in doc.findall(".//reaction-list/Reaction"):
        print(e.attrib['frameid'])
        reactions.append(e.attrib['frameid'])
    return reactions

def getPathwayIdsToReactions(pathways):
    """
    Associate for each identified metabolic pathway the 
    reactions that are associated with it.
    
    Parameters
    ----------
    pathways (list): List of identified pathways

    ----------
    Returns : data (dict): String of pathway id as key and its 
              associated reaction list as values
        
    """
    data = {}
    for p in pathways:
        try :
            reactions_list = getReactionIdsFromPathway(p)
            data[p] = reactions_list
        except :
            pass
    return data

def removeNotBiocycReactions(graph, notBR):
    """
    Deletes all reactions not found on BioCyc. Substrates are 
    also deleted if they are not used by found reactions.
    
    Parameters
    ----------
        graph (tlp.Graph): Ecoli K12 substrates - reactions graph 
        notBR (list): List of reaction IDs not found on BioCyc
    
    """
    biocycReactionNodes = hg.filterNodesWithIds(graph, notBR, excluded = True)
    nodeIdsToKeep = set(hg.getIdsFromNodes(graph, biocycReactionNodes))
    for n in biocycReactionNodes:
        neighborNodes = hg.getNeighborNodes(graph, n)
        nodeIdsToKeep.update(hg.getIdsFromNodes(graph, neighborNodes))
    hg.deleteNodes(graph, hg.filterNodesWithIds(graph, nodeIdsToKeep, excluded = True))

def getReactionIdsFromGraph(graph):
    """
    Returns all BioCyc ids corresponding to a reaction in the graph.
    
    Parameters
    ----------
        graph (tlp.Graph): Ecoli K12 substrates - reactions graph 
    
    ----------
    Return : reactions (list): List of associated reactions ids
    """
    ids=graph.getStringProperty('id')
    isReaction=graph.getBooleanProperty('reaction')
    reactions=[]
    for n in graph.getNodes():
        if isReaction[n]:
            reactions.append(ids[n])
    return reactions

def filterBiocycPathways(graph):
    """
    Removes the reactions not found on BioCyc and returns a dict of pathways as keys
    and their reactions as values.
    
    Parameters
    ----------
        graph (tlp.Graph): Ecoli K12 substrates - reactions graph 

    ----------    
    Return : data (list): dict of pathways (keys) and their associated reactions ids (values)
    """
    # measure query time
    global start
    start = time.time()
    reactions = getReactionIdsFromGraph(graph)
    pathways, reactions_to_del = getPathways(reactions)
    data = getPathwayIdsToReactions(pathways)
    removeNotBiocycReactions(graph, reactions_to_del)
    return data
