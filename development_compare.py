"""
A script to compare Phallusia individuals developments.
ASLOUDJ Yanis
15/02/2022
"""

from xml.etree import ElementTree as ET

properties = [
    'all_cells',
    'cell_history',
    'cell_barycenter',
    'cell_contact_surface',
    'cell_labels_in_time',
    'problematic_cells',
    'cell_lineage',
    'cell_volume',
    'cell_name',
    'cell_fate',
    'tissuefate_guignard_2020',
    'tissuefate_lemaire_2009'
]

filename = 'Astec-pm1.xml'
doc = ET.parse(filename)

def printProperty(propertyName):
    """Prints the children of a given property on the .XML file.
    Temporary function."""
    nEntries = 0
    for child in doc.findall(f".//{propertyName}/cell"):
        nEntries += 1
        print(child.attrib, child.text)
    print(nEntries)

printProperty('cell_fate')