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
    ID : string
        the id of the requested object
    
    Returns
    -------
    doc : ElementTree.Element
        The BioCyc XML text of our requested object 
    """
    # monitor requests 
    URL = "https://websvc.biocyc.org/getxml?ECOLI:" + ID
    response = requests.get(URL, timeout = 5)
    print(f"""{round(time.time() - start, 2)}: {ID} ({response.status_code})""")
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
    pathways: list
        the list of pathways associated to our reactions 
        in BioCyc
    reaction: string
        the id of the reaction for which we search its
        pathways in BioCyc
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
    reactions: lists
        the list of reaction IDs retrieved from our Tulip
        graph
                          
    Returns
    -------
    pathways: list
        the list of pathways associated to our reactions
        in BioCyc
    to_delete: list
        the list of reactions to delete in our Tulip graph
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
    pathway : str
        the pathway id
        
    Returns
    -------
    reactions : list
        the list of associated reactions ids
    """
    reactions = []
    doc = requestBiocyc(pathway)
    for e in doc.findall(".//reaction-list/Reaction"):
        reactions.append(e.attrib['frameid'])
    return reactions

def getPathwayIdsToReactions(pathways):
    """
    Associate for each identified metabolic pathway the 
    reactions that are associated with it.
    
    Parameters
    ----------
    pathways: list
        the list of identified pathways
    
    Returns
    -------
    data: dict
        string of pathway id as key and its associated
        reaction list as values
    """
    data = {}
    for p in pathways:
        try :
            reactions_list = getReactionIdsFromPathway(p)
            data[p] = reactions_list
        except :
            pass
    return data

def getNodeIdsFromGraph(graph, targetReaction):
    """
    Returns all BioCyc ids corresponding to a substract
    or a reaction in the graph.

    Parameters
    ----------
    graph: tlp.Graph 
        the Ecoli K12 substrates - reactions graph 
    targetReaction: bool
        if True, returns the reactions, 
        if False, returns the substracts
    
    Returns
    -------
    substracts: list
        the list of associated ids
    """
    ids = graph.getStringProperty('id')
    isReaction = graph.getBooleanProperty('reaction')
    targetIds = []
    for n in graph.getNodes():
        if isReaction[n] == targetReaction:
            targetIds.append(ids[n])
    return targetIds

def removeNotBiocycReactions(graph, notBR):
    """
    Deletes all reactions not found on BioCyc. Substrates are 
    also deleted if they are not used by found reactions.
    
    Parameters
    ----------
    graph: tlp.Graph 
        the Ecoli K12 substrates - reactions graph 
    notBR: list 
        the list of reaction IDs not found on BioCyc
    """
    substracts = getNodeIdsFromGraph(graph,
                                     targetReaction = False)
    biocycReactionNodes = hg.filterNodesWithIds(graph,
                                                notBR + substracts,
                                                excluded = True)
    nodeIdsToKeep = set(hg.getIdsFromNodes(graph,
                                           biocycReactionNodes))
    for n in biocycReactionNodes:
        neighborNodes = hg.getNeighborNodes(graph, n)
        nodeIdsToKeep.update(hg.getIdsFromNodes(graph,
                                                neighborNodes))
    hg.deleteNodes(graph, hg.filterNodesWithIds(graph,
                                                nodeIdsToKeep,
                                                excluded = True))

def filterBiocycPathways(graph):
    """
    Removes the reactions not found on BioCyc and returns
    a dict of pathways as keys and their reactions as values.
    
    Parameters
    ----------
    graph: tlp.Graph
        the Ecoli K12 substrates - reactions graph 
        
    Returns
    -------
    data: dict
        dictionnary of pathways (keys) and their associated
        reactions ids (values)
    """
    # measure query time
    global start
    start = time.time()
    reactions = getNodeIdsFromGraph(graph,
                                    targetReaction = True)
    pathways, reactions_to_del = getPathways(reactions)
    data = getPathwayIdsToReactions(pathways)
    removeNotBiocycReactions(graph, reactions_to_del)
    print(f'''{len(reactions_to_del)} reactions not
          found on BioCyc and deleted.''')
    return data
