from gmshparser.abstract_parser import AbstractParser

class EntitiesParser(AbstractParser):
    @staticmethod
    def get_section_name():
        return "Entities"

    @staticmethod
    def parse(mesh, io) -> None:
        import pdb; pdb.set_trace()
        s = io.readline().strip().split(" ")
        mesh.set_version(float(s[0]))
        mesh.set_ascii(int(s[1]) == 0)
        mesh.set_precision(int(s[2]))