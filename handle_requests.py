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

def requestBiocyc(reactionId, list):
    URL = "https://websvc.biocyc.org/getxml?ECOLI:" + reactionId
    response = requests.get(URL, timeout = 2)
    if response.status_code not in list:
        list.append(response.status_code)
    if response.status_code == 200:
        doc = ET.fromstring(response.text)
    elif response.status_code == 429:
        print(response.headers)
        time.sleep(60)
        response = requests.get(URL, timeout = 2)
        doc = ET.fromstring(response.text)
    return doc, list

def getPathwayID(pathways, reaction, list):
    doc, liste = requestBiocyc(reaction, list)
    for e in doc.findall(".//in-pathway/Pathway"):
        ID = e.attrib['frameid']
        if ID not in pathways:
            pathways.append(ID)
   
def getPathways(reactions):
    pathways = []
    to_delete = []
    liste = []
    s = requests.Session()
    for r in reactions:
        try:
            getPathwayID(pathways, r, liste)      
        except :
            to_delete.append(r)
    return pathways, to_delete, liste

def getReactionsID(pathway):
    """
    Retrieves in BioCyc all reactions associated to a pathway.
    
    Parameters:
        pathway (string) : Pathway id
        
    Returns:
        reactions (list) : List of associated reactions
        
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
        pathways (list) : List of identified pathways
    
    Returns:
        data (dict) : String of pathway id as key and its 
                      associated reaction list as values.
        
    """
    data = {}
    for p in pathways:
        try :
            reactions_list = getReactionsID(data, p)
            data[p] = reactions_list
        except :
            pass
    return data

def subgraph(graph,pathway_name,dictionnaire):
    list_reactions = dictionnaire[pathway_name]
    id_list=graph.getStringProperty('id')
    name_list=graph.getStringProperty('name')
    pathwayNodes =[]
    for n in graph.getNodes():
        if id_list[n] in list_reactions:
            pathwayNodes.append(n)
            for n2 in graph.getInOutNodes(n):
                pathwayNodes.append(n2)
    pathwaySg = graph.inducedSubGraph(pathwayNodes)
    pathwaySg.setName(pathway_name)
    
def filterBiocycReactions(graph, notBR):
    """
    Supprime l'ensemble des reactions non trouvees sur BioCyc.
    Les substrats sont egalement supprimes s'ils ne sont pas 
    utilises par des reactions trouvees.
    
    Parameters:
        graph (tlp.Graph) : 
        notBR (list) : List of reaction IDs not found on BioCyc
    
    """
    biocycReactionNodes = hg.filterNodesWithIds(graph, notBR, excluded = True)
    nodeIdsToKeep = set(hg.getIdsFromNodes(biocycReactionNodes))
    for n in biocycReactionNodes:
        neighborNodes = hg.getNeighborNodes(graph, n)
        nodeIdsToKeep.update(hg.getIdsFromNodes(neighborNodes))
    hg.deleteNodes(graph, hg.filterNodesWithIds(graph, nodeIdsToKeep, excluded=True))
