![DARSim Logo](Logo/DARSimLogo2.png)

Delft Advanced Reservoir Simulator (TU Delft)

Reservoir simulation framework for CCS and UHS in geological porous media

---

## Overview

DARSim (Delft Advanced Reservoir Simulator) is a research-focused reservoir simulation framework developed at TU Delft. It is written in MATLAB and designed for modeling multiphase flow and transport in porous geological media, with applications in carbon capture and storage (CCS) and underground hydrogen storage (UHS).
DARSim is based on a fully implicit finite volume discretization framework and supports advanced numerical methods including:

- Algebraic Dynamic Multilevel (ADM) methods  
- Projection-based Embedded Discrete Fracture Models (pEDFM)

The simulator is designed for flexibility in grid representation and physical modeling, enabling research into complex subsurface flow systems.

![DARSim GIF](Logo/DARSimLogo.gif)

---

## Key Features

### Numerical Methods
- Fully implicit finite volume formulation
- ADM (Algebraic Dynamic Multilevel) methods
- pEDFM (Projection-based Embedded Discrete Fracture Modeling)

### Grid Support
- Cartesian grids
- Corner-point grids
- General unstructured meshes (externally generated)

### Fluid Models
- CO₂–brine systems (black-oil and compositional)
- H₂–brine systems (black-oil and compositional)
- H₂–CH₄–brine (fully compositional, 3-component systems)

### Physics Options
- Isothermal and non-isothermal simulations
- Diffusion and dispersion (Cartesian grids)

---

## Installation

To install DARSim:

1. Clone repository to your local machine.
2. Open MATLAB.
3. Add `/src` and all subdirectories to the MATLAB path.
4. Execute DARSim from within the `/src` directory.

---

## Running Simulations

DARSim can be executed either via a graphical user interface (GUI) or directly from the MATLAB command window.

### Command-line execution

```matlab
DARSim("<path-to-simulation-files>", "<simulation-file-name.txt>")
```

### Example

DARSim("examples/<example_name>", "Sim.txt")

Each simulation requires at minimum:

- A main simulation input file (e.g., Sim.txt as found in the included example problems)
- A SimulatorSettings.txt file defining global simulation parameters

### Input and Output
1. Inputs
- Simulation configuration files (TXT format)
- Grid and geometry definitions
- Fluid and rock property definitions
- Boundary and initial conditions, including wells
2. Outputs
- VTK files for visualization
- Results compatible with ParaView

### Example Problems

DARSim includes a suite of example cases demonstrating various scenarios, such as:

- Basic multiphase flow
- CO₂ injection scenarios using ADM
- 3-component hydrogen storage systems
- Fractured reservoir simulations

---

### Documentation

Please refer to the Documentation directory for full documentation. 

An overview of code organization and structure is provided in:

```matlab
/documentation/code-organization.html
```

---

### License and Usage Disclaimer

DARSim is provided for **research use only**.

The authors and Delft University of Technology assume **no liability** for the use of this software or its results.

### Principal Investigators

- Hadi Hajibeygi (TU Delft, Civil Engineering and Geosciences)
- Cornelis Vuik (TU Delft, Delft Institute of Applied Mathematics)

### Contributors List

Current and previous members of the DARSim development team include (but is not limited to):

- Mengjie Zhao
- Amin Masaghi Bonabi
- Ryan Haagenson
- Yuhang Wang
- Artur Castiel de Souza
- Mousa HosseiniMehr
- Janio Piguave
- Matei Ţene
- Matteo Cusini

### Citations

If you use DARSim in academic work, please cite the relevant publications associated with the project.

### Contact

For questions or collaboration inquiries, please contact the DARSim development team at TU Delft.

