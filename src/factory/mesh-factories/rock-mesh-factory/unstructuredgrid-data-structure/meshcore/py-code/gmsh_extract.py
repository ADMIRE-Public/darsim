
import gmshparser
from load_parsers import parse

# from gmshparser.mesh_format_parser import MeshFormatParser
# from gmshparser.nodes_parser import NodesParser
# from gmshparser.elements_parser import ElementsParser
# from gmshparser.abstract_parser import AbstractParser
# from custom_parser import EntitiesParser

import numpy as np



class GmshCore:    
    def __init__(self, filepath):
        
        mesh = parse(filepath)
        self.mesh = mesh
        # import pdb; pdb.set_trace()
        # flags = EntitiesParser()
        # fop = flags.parse(mesh);

        self.coords, self.node_tag = self.get_coords(mesh)
        self.connectivities, self.el_dim, self.el_types, self.el_tag, self.el_id = self.get_elements(mesh)
        
        print("sucess")
    
    def get_elements(self,mesh):
        entities_types = []
        entities_dim = []
        entities_id = []
        entities_tag = []
        connectivity_matrix = []
        for index, entity in enumerate(mesh.get_element_entities()):
            local_entities_id = []
            local_connectivity_matrix = []
            eltype = entity.get_element_type()
            entities_types.append(eltype)
            entities_dim.append(entity.get_dimension())
            entities_tag.append(entity.get_tag())
            # print("Element type: %s" % eltype)
            # import pdb; pdb.set_trace()                    
            for element in entity.get_elements():
                elid = element.get_tag()
                # import pdb; pdb.set_trace()
                # entities_num_tag.append(len(elid))
                local_entities_id.append(elid)
                elcon = element.get_connectivity()
                local_connectivity_matrix.append(elcon)
                # print("Element id = %s, connectivity = %s" % (elid, elcon))    
            connectivity_matrix.append(np.array(local_connectivity_matrix))
            entities_id.append([local_entities_id])
        return connectivity_matrix, np.array(entities_dim), np.array(entities_types), np.array(entities_tag), entities_id
    
    
    def get_coords(self,mesh):
        coords = []
        node_tag = []
        for entity in mesh.get_node_entities():
            for node in entity.get_nodes():
                nid = node.get_tag()
                ncoords = node.get_coordinates()
                coords.append(ncoords)
                node_tag.append(nid)
                # print("Node id = %s, node coordinates = %s" % (nid, ncoords))
        return np.array(coords), np.array(node_tag)



