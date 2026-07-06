import startup
import sys
from pathlib import Path



def convert_mesh(mesh_name):
    file_name = Path(mesh_name)
    if file_name.exists():
        import gmsh
        gmsh.initialize()
        gmsh.open(str(file_name))
        gmsh.write(str(file_name.with_suffix(".m")))
        print("Matlab Mesh Successfuly created.")
    else:
        print("Mesh file not found.")
    

def main():
    print("Gmsh to Matlab Mesh Converter")
    print(f"Converting Mesh: {sys.argv[1]}")    
    startup.check_and_install_requirements()
    convert_mesh(sys.argv[1])
    
if __name__ == "__main__":
    main()