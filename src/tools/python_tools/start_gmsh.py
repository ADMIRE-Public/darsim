# start_gmsh.py
"""Python module demonstrates passing MATLAB types to Python functions"""

import importlib

dependency = 'gmsh';
def check_and_install_requirements():
    try:
        importlib.import_module(dependency)
        print("All dependencies installed.")
    except ImportError:
        print("Installing dependencies...")
        try:
            import pip
            pip.main(['install', dependency])
        except ImportError:
            print("Failed to install dependencies. Please make sure you have pip installed.")


def update_pip():
    try:
        import pip
    except ImportError:
        print("Failed to install dependencies. Please make sure you have pip installed.")
        pip.main(['install', "--upgrade pip"])



# Calling the function to check and install NumPy

#check_and_install_requirements()
