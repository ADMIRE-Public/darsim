from gmshparser.mesh import Mesh
from gmshparser.main_parser import MainParser
from gmshparser.mesh_format_parser import MeshFormatParser
from gmshparser.nodes_parser import NodesParser
from gmshparser.elements_parser import ElementsParser
from custom_parsers import EntitiesParser

PARSERS_LIST = [
    MeshFormatParser,
    NodesParser,
    ElementsParser,
    EntitiesParser,
]


def parse(filename: str) -> Mesh:
    """Parse Gmsh .msh file and return `Mesh` object."""
    mesh = Mesh()
    mesh.set_name(filename)
    parser = MainParser(PARSERS_LIST)
    with open(filename, "r") as io:
        parser.parse(mesh, io)
    return mesh


