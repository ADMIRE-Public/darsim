
# DARSim User Guide

![DARSim Logo](../Logo/DARSimLogo2.pdf)

## Table of Contents

- [Introduction](#1-introduction)
- [Quick Setup Guide](#2-quick-setup-guide)
- [Mathematical Model](#3-mathematical-model)
- [Discretization](#4-discretization)
- [Reservoir Properties](#5-reservoir-properties)
- [Fluid Models](#6-fluid-models)
- [Thermal Properties](#7-thermal-properties)
- [Initial Conditions](#8-initial-conditions)
- [Boundary Conditions](#9-boundary-conditions)
- [Wells](#10-wells)
- [Fractured Reservoirs](#11-fractures)
- [ADM](#12-algebraic-dynamic-multilevel-adm)
- [Output Controls](#13-output-controls)
- [Visualization](#14-visualization)

---

## 1. Introduction

DARSim (Delft Advanced Reservoir Simulator) is an open-source, research-focused reservoir simulation framework developed at TU Delft for modeling multiphase flow and transport in porous geological media. The simulator is written in MATLAB and is primarily intended for research in carbon capture and storage (CCS), underground hydrogen storage (UHS), and related subsurface flow applications.

DARSim employs a fully implicit finite-volume formulation and supports structured, corner-point, and general unstructured grids. Advanced numerical capabilities include Algebraic Dynamic Multilevel (ADM) methods and projection-based Embedded Discrete Fracture Models (pEDFM).

This guide describes the mathematical and numerical models implemented in DARSim, explains how to install and run the simulator, and provides an overview of the main physical models available to users. Internal implementation details and source code organization are described separately in the Developer Guide.

---

## 2. Quick Setup Guide

### 2.1 Installation

1. Clone or download the DARSim repository from GitHub.
2. Open MATLAB.
3. Add the `src` directory from the DARSim repository (with all subdirectories) to the MATLAB environment path. This can be done with the "Set Path" button in the MATLAB window.
4. Change the current working directory to `src` before running any simulations.

### 2.2 Running DARSim

DARSim can be executed directly from the MATLAB command window:

```matlab
DARSim("<path/to/simulation>","<InputFileName.txt>")
```

where the correct path to the simulation directory and the name of the main simulation input file are given as the two input arguments. All simulations require, at minimum:

- A main simulation input file (can have any name).
- A `SimulatorSettings.txt` file containing the global simulation settings.

More details on the contents of these files is provided below. Please refer to the `Examples` directory of the DARSim repository for examples of the simulation files. Note that the main simulation input file is named `Sim.txt` in the `Examples` directory.

### 2.2 Parallel Simulations

On Windows and Linux machines, DARSim can be run in parallel. This capability makes use of the native parallelized solver capabilities available through MATLAB. To run a simulation in parallel, pass the number of desried processors (`n`) to the DARSim function:

```matlab
DARSim("<path/to/simulation>","<InputFileName.txt>",n)
```

---

## 3. Mathematical Model

DARSim solves the conservation equations governing multiphase, multicomponent flow in porous media using a finite-volume formulation.

### 3.1 Mass Conservation

For each fluid component ($c$), DARSim solves a nonlinear mass conservation equation of the form

$$
\frac{\partial}{\partial t}\left(\phi S_\alpha x_{c,\alpha}\rho_\alpha\right) + \nabla \cdot \left(\sum_{\alpha} x_{c,\alpha}\rho_{\alpha}u_{\alpha}\right) + \nabla \cdot \left(\sum_{\alpha} S_{\alpha}\rho_{\alpha}J_{\alpha}\right) - \sum_{\alpha} x_{c,\alpha}\rho_\alpha(q_{w,\alpha}+q_{f,\alpha}) = 0
$$

where:

- $\phi$ is the porosity
- $S_\alpha$ is the saturation of phase $\alpha$
- $\rho_\alpha$ is the density of phase $\alpha$
- $x_{c,\alpha}$ is the component mass fraction of $c$ in phase $\alpha$
- $u_{\alpha}$ is the darcy flux of phase $\alpha$
- $J_\alpha$ is the diffusion/dispersion flux of phase $\alpha$
- $q_w$ is the well source/sink term
- $q_f$ is the fracture-matrix or fracture-fracture exchange terms

### 3.2 Darcy Flow

Phase velocities are computed using Darcy's law,

$$
u_\alpha = - \lambda_\alpha (\nabla p_\alpha - \rho_\alpha g \nabla z)
$$

where $\lambda_\alpha$ is the phase mobility. In DARSim, Darcy's law is approximated with an upwinded, two-point flux approximation in the discretized formulation. 

### 3.3 Phase Behaviour

Depending on the selected fluid model, DARSim performs flash calculations to determine phase compositions and updates fluid properties including density, viscosity, saturation, capillary pressure, and thermodynamic quantities required for compositional and thermal simulations.

---

## 4. Discretization

### 4.1 Spatial Discretization

DARSim uses a finite-volume discretization in which the governing conservation equations are integrated over control volumes. Fluxes between neighboring cells are evaluated with a two-point flux approximation, using transmissibilities together with an upwind treatment of phase mobilities.

The simulator supports:

- Cartesian grids;
- Corner-point grids;
- General unstructured meshes.

Regardless of grid type, the governing equations remain unchanged; only the geometric representation and connectivity differ. The following subsections describe how to use each grid type. 

#### 4.1.1 Cartesian Grids

To employ a cartesian grid, the main simluation input file should include:

```text
GRID_DISCRETIZATION
CartesianGrid

DIMENSION
[<DX>, <DY>, <DZ>]

SPECGRID
[<NX>, <NY>, <NZ>]
```

where DX, DY, and DZ are the spacial extent in the x, y and z directions and NX, NY, and NZ are the number of grid blocks to use in each direction. 

#### 4.1.2 Corner-Point Grids

There are two options for using a corner-point grid. 
 
#### Option A: Read an external CPG file (most common)

For this option, the main simulation input file should include:

```text
GRID_DISCRETIZATION
CornerPointGrid
my_model.grdecl

DIMENSION
[<DX>, <DY>, <DZ>]

SPECGRID
[<NX>, <NY>, <NZ>]
```

or a DARSim corner-point `.txt` file.

#### Option B: Create corner-point grid based on a Cartesian grid (for debugging)

In the main simulation input file, include:

```text
GRID_DISCRETIZATION
CornerPointGrid

DIMENSION
[<DX>, <DY>, <DZ>]

SPECGRID
[<NX>, <NY>, <NZ>]
```

In this approach, `build_cornerpoint()` interprets the setup as a `CartesianBox` and generates `CornerPointGridData` based on the CartesianGrid described. Often, the purpose of this approach is to debug any development work using the corner-point grid capability. 



#### Note: If you also want rock properties from the CPG file...

From the implementation of `build_cornerpoint.m`, the code also supports the following input block for `GRID_DISCRETIZATION`:

```text
GRID_DISCRETIZATION
CornerPointGrid
my_model.grdecl
INCLUDE_ROCK_PROPERTIES
```

This allows permeability and porosity contained in the CPG grid file to be included in the DARSim simulation.

#### Note: Cache behavior for CPG...

After the first successful simulation using a corner-point grid, the code will generate a file called:

```text
CornerPointGridData.mat
```

If you later make changes to the `grdecl` file, delete this `.mat` file before rerunning the simulation. Otherwise, the cached information about the corner-point grid will be used instead of the updated `grdecl` file. 

#### 4.1.3 Unstructured Grids

There are two options for using the unstructured grid capability in DARSim, similar to corner-point grids. 

#### Option A: Read an external Unstructured grid file (most common)

In this approach, the main simulation input file should include:
```text
GRID_DISCRETIZATION
UnstructuredGrid
unstructured_mesh.m

DIMENSION
[<DX>, <DY>, <DZ>]

SPECGRID
[<NX>, <NY>, <NZ>]
```

where `unstructured_mesh.m` is a helper script that must be included in the simulation directory. This script points to the extrenal unstructured mesh file (e.g. `uns_mesh.m`). For more information about running a simulation with an unstructure mesh, refer to the unstructured simulation examples provided in the `Examples` directory. 

Note that `DIMENSION` and `SPECGRID` are still required keywords in the main simulation input file, even though they are not used by the simulation when the grid is defined by an external file. In this case, use dummy values for the `DIMENSION` and `SPECGRID` keywords.

#### Option B: Unstructured grid based on a Cartesian grid (for debugging)

```text
GRID_DISCRETIZATION
UnstructuredGrid

DIMENSION
[<DX>, <DY>, <DZ>]

SPECGRID
[<NX>, <NY>, <NZ>]
```

In this approach, DARSim first creates an MRST `cartGrid` object, then wraps it through `MeshCore + translateMeshCore` into an unstructured grid based on the Cartesian grid desscribed in the main simulation input file. Often, the purpose of this approach is to debug any development work using the corner-point grid capability. 

### 4.2 Time Discretization

DARSim employs a fully implicit time discretization. At each time step:

1. Fluid properties are updated.
2. The nonlinear residual equations are assembled.
3. Newton iterations are performed until convergence.
4. The timestep size is adjusted based on convergence behaviour and reporting requirements.

The total simulation time is defined in the main simulation input file using the keyword `TOTAL_TIME`. For example:

```
TOTAL_TIME
	1000 DAYS 0 HOURS 0 MINUTES 0 SECONDS
```

This would set up a 1,000 day simulation. 

The time step size is dynamic in the simulation and based on the convergence behavior of the previous time step or the need for shorter time steps to avheive convergence of the Newton solver. The minimum and maximum time step sizes are included in the `SimulatorSettings.txt` file through the keyword `MINMAXDT`. For example, you could include the following block in the `SimulatorSettings.txt` file:

```
MINMAXDT
	1e-2 5.00e0
```

This defines the minimum and maximum time step size to be 1e-2 and 5 days, respectively. 

---

## 5. Reservoir Properties

Reservoir properties define the physical characteristics of the porous medium. Typical properties include:

- Porosity
- Permeability
- Relative permeability curves
- Capillary pressure curves

Each of these properties are included in the main simluation input file. 

### 5.1 Porosity

To define a single porosity value for the entire domain, simply include a single porosity value (between 0 and 1) with the `POROSITY` keyword, such as

```
POROSITY
    0.2
```

If you'd like to define a heterogeneous porosity field, this is done with an external porosity file. Instead of a single porosity value, the `POROSITY` keyword is given a file name. This file must be present in the simulation directory for DARSim to find it.

```
POROSITY
    het_porosity.txt
```

### 5.2 Permeability

For permeability, the units (milliDarcy or square meters) must be defined through the `PERMEABILITY_UNIT` keyword. Next, the scale (`Linear` or `Logarithmic`) must be defined with `PERMEABILITY SCALE`.

```
PERMEABILITY_UNIT
    m2

PERMEABILITY_SCALE
    Linear
```

Permeability values can be define as either homogeneous (single value) or heterogeneous (from and external file), similar to the porosity value. Permeablity, though, can also be anisotropic with different values in the x, y, and z direcitons. For example:

```
PERMEABILITY_X
    8.0e-14


PERMEABILITY_Y
    8.0e-14


PERMEABILITY_Z
    8.0e-14
```

To use external files for the permeability values (to account for heterogeneity), use the keywords as follows:

```
PERMEABILITY_X
	INCLUDE
    Perm_het_x.txt

PERMEABILITY_Y
	INCLUDE
    Perm_het_y.txt

PERMEABILITY_Z
	INCLUDE
    Perm_het_z.txt
```

### 5.3 Relative Permeability and Capillary Pressure Curves

Generally, the relative permeability and capillary pressure curves are supplied to DARSim with a lookup table. The RELPERM and CAPILARITY keywords then appear in the main simulation input file as follows:

```
RELPERM
	Table
	<LookupTableName.txt>

CAPILLARITY
	Table
	2
```

Here, the `2` under `CAPILARITY` tells DARSim which fluid phase is considered the wetting phase. Typically this is the brine phase (i.e. phase `2`). 

The lookup table file can have any name (given in place of `<LookupTableName.txt>`) and located in the simulation directory for DARSim to find it. The look up table will include five separate tables:

- Wetting phase relative permeability curve: --wetting_kr
- Non-wetting phase relative permeability curve (drainage): --nonwetting_drainage_kr
- Non-wetting phase relative permeability curve (imbibition): --nonwetting_imbibation_kr
- Capillary pressure curve (drainage): --Pc_drainage
- Capillary pressure curve (imbibition): --Pc_imbibation

Each table consists of two columns. The first column is the saturation value for a given data point and the second column is the associated relative permeability of capillary pressure value. For an example of the `LookupTable.txt` file, please refer to an example in the `Examples` directory of the DARSim repository. 

To use hysteresis, include the keyword `HYSTERESIS` in the main simulation input file:

```
HYSTERESIS
    ON
```

For cases without hysteresis, set this keyword to `OFF` and give identical curves for the drainage and imbitition curves inside `LookupTable.txt`. 

When including fractures in the simulation, the relative permeability model also needs the wetting-phase residual saturation value in the fracture cells. To include this, modify the `RELPERM` block to:

```
RELPREM
	Table
	<LookupTableName.txt>
	srw_frac
	0.2
```

### 5.4 Multiple Facies Types

It is also possible to include multiple lookup tables in DARSim, one for each rock/facies type in the simulation. The rock/facies types must first be given to DARSim with the keyword `FACIES_INDEX` in the main simulation input file, as follows:

```
FACIES_INDEX
    <facies_file.txt>
```

Here, the user should put the actual name of the facies file, which should be located in the simulation directory for DARSim to find. This file is a single column, where the row number corresponds to the global cell index in the rock mesh and the value given in a given row provides the rock/facies type to use for that rock matrix cell. There can be an arbitrary number of rock/facies types; however, each type requires a unique lookup table for the relative permeability and capillary pressure curves. If there are two rock/facies types, for example, the main simulation input file should include:

```
RELPERM
	Table
    <table_1.txt>
    <table_2.txt>

CAPILLARITY
	Table
	2
```

Here, `<table_1.txt>` and `<table_2.txt>` are standing in for the actual lookup table names to use for each rock/facies type. 

---

## 6. Fluid Models

DARSim supports several fluid formulations to address a wide range of subsurface flow problems related to CCS and UHS.

- CO₂–brine systems (black-oil and fully compositional)
- H₂–brine systems (black-oil and fully compositional)
- H₂–CH₄–brine systems (fully compositional)

The appropriate fluid model is selected through the main simulation input file with the keywords FLUID_MODEL and COMPONENT_LIST. For example, a compositional CO₂–brine system would be set with

```
FLUID_MODEL
	Compositional
	--Number of phases
	2
	--Number of components
	2

COMPONENT_LIST
	CO2-Brine
```

### 6.1 Fluid properties: Density, Viscosity, Salinity

The governing equations depend on two important fluid properties for each fluid phase: density and viscosity. 

For all possible fluid models (black-oil and compositional), fluid density is determined by the associated thermodynamic model employed by DARSim. Therefore, fluid density is not given to DARSim as an input parameter. (Note: DARSim assumes the density of the brine phase to be 1060 kg/m<sup>3</sup> at atmospheric conditions.)

Viscosity, however, is not generally determined by the thermodynamic model in DARSim and is taken as a constant value provided in the main simulation input file with the keyword FLUID_VISCOSITY:

```
FLUID_VISCOSITY (Pa.sec)
	6.5e-5
	5.4e-4
```

In this example, the viscosity of the gas phase (the first phase) is 6.5e-5 Pa.sec, while the viscosity of the brine phase (the second phase) is 5.4e-4 Pa.sec. The units for viscosity in the input file are assumed to be Pascal-seconds. 

For the three component system of H2-CH4-Brine, the viscosity of the gas phase can vary significantly depending on the composition. Therefore, DARSim employs a composition-dependent relationship to define the gas phase viscosity in the three component fluid model. In this case, the gas phase viscosity in the main simulation input file can be listed as "NaN". 

The salinity of the brine phase is controlled through the following keyword in the main simulation input file:

```
SALINITY (ppm by mass)
	100e3
```

---

## 7. Thermal Properties

DARSim allows for non-isothermal simulations for the CO₂–brine compositional fluid model. To use a non-isothermal simulation, define the fluid model as `Compositional_Energy` with the FLUID_MODEL keyword. 

A few additional fluid and reservoir properties are required for the non-isothermal simulations. These include:

- Heat conducvitiy for each fluid phase (FLUID_HEAT_CONDUCTIVITY)
- Specific heat for each fluid phase (FLUID_SPECIFIC_HEAT)
- Heat conductivity for the rock matrix (ROCK_HEAT_CONDUCTIVITY)
- Specific heat for the rock matrix (ROCK_SPECIFIC_HEAT)
- Mass density of the rock matrix (ROCK_DENSITY)

Here, the relevant keyword to include in the main simulation input file is shown in the parentheses. The units for thermal properties are SI: W/m/K for conductivity and J/Kg/K for specific heat. 

---

## 8. Initial Conditions

Initial conditions specify the reservoir state at the beginning of the simulation. Accurate initialization is essential for obtaining physically meaningful transient solutions. Depending on the selected physics, initial conditions may include:

- Pressure
- Phase saturations
- Component compositions
- Temperature

Initial pressure in the reservoir is typically treated as a constant value. For example, the follow keyword in the main simulation input file produces an initial pressure of 2.5e7 Pa. 

```
INITIAL_PRESSURE
	2.5e7
```

Similarly, the saturation of each phase can be set with the following:

```
INITIAL_SATURATION_1
	0.0
INITIAL_SATURATION_2
	1.0
```
while the overall mass fraction of each component (if using a compositional fluid model) is set using:

```
INITIAL_COMPOSITION_1
	0.0
INITIAL_COMPOSITION_2
	1.0
```

Alternatively, the initial conditions can be set using an external file. In this way, uniform values do not need to be applied across the entire reservoir. To use an external file for initial conditions, the main simulation input file should include:

```
INITIAL_STATE
    <Initial_State.txt>
```

Here, `Initial_State.txt` is the file name containing the initial conditions, which can have any name. The comma-separated text file should include one column for the cell index, one column for the brine phase pressure, and one column for the mass fractions of each component in the system. For example, a two-component CO2-Brine system would have an initial state file:

```
1, 1.0e7, 0.1, 0.9
2, 1.1e7, 0.2, 0.8
3, 1.5e7, 0.5, 0.5
...
```

There should be one row in the initial state file for every cell in the grid. Note that the mass fraction values must sum to unity. 

When the initial conditions are supplied from an external file, the simulation still requires a value for the `INITIAL_PRESSURE` keyword. This is a single pressure value that will be used to initialize the fluid properties of the injection fluid at the wells; however, fluid properties will rapidly come into equilibirum with the reservoir fluid pressure during the simulation. 

#### Note: Component and phase indexing...

For the black-oil fluid models, the first phase is the gaseous phase and the first component is CO2 or H2 (depending on the composition list) while the second phase and component is brine. For the three-component H2-CH4-Brine system, the first two components are H2 and CH4 respectively, while the last component is brine. 

---

## 9. Boundary Conditions

By default, DARSim applies a no-flux boundary condition on all external boundaries. To replicate a constant pressure boundary condition, it is suggested to use either a pore volume multiplier on the boundary cells or place constant pressure injection or production wells in the appropriate locations along the domain boundary. (More information about wells is given below.)

---

## 10. Wells

DARSim models injection and production wells through source and sink terms in the governing conservation equations and implemented with the well-known Peaceman correction term. Well models are coupled fully implicitly with the reservoir equations, allowing stable simulation of strongly coupled flow processes.

Default well controls include:

- Fixed bottom-hole pressure (BHP)
- Fixed injection or production rate

An example of a well definition in the main simulation input file is given below:

```
WELL_START
	TYPE ("INJ" or "PROD")
		INJ
	COORDINATE ("IJK_LIST" or "XYZ_LIST" or "CELL_INDEX_LIST")
		XYZ_LIST
		[ <x1> , <y1> , <z1> ]
		[ <x2> , <y2> , <z2> ]
	CONSTRAINT ("PRESSURE" or "RATE", Units: Pa or Pore Volumes/Day)
		RATE
		1e-1
	FORMULA (Well Index: "WI")
		WI
		9.2
	SCHEDULE ([Open, Close, Open, Close, ...] in days)
		0, 30
	COMPOSITION (Which component to inject -- only required for injectors)
		1
	BOUNDARY_CONDITION
		TEMPERATURE (K)
		300
WELL_END
```
Here, we have defined an injection (`INJ`) well that is perforated along the line segment with end points (x1, y1, z1) and (x2, y2, z2). Alternatively, the well could be defined as a production (`PROD`) well with end points being defined using the I, J, and K indices of the cells within a Cartesidan grid or with the individual cell indices of a corner-point or unstructured grid. The well is defined as rate (`RATE`) controlled, opposed to pressure (`PRESSURE`) controlled, with a rate (`RATE`) of 1e-1 pore volumes per day. (Note: "Pore volumes" reffers to the total pore volume in the computational domain.) The so-called "Well Index" (`WI`) of the Peaceman well model is given as 9.2. The well has an optional `SCHEDULE` that can be defined. Here, the well opens at day 0 and closes on day 30. If the `SCHEDULE` keyword is not included in the well definition, DARSim will assume the well remains open throughout the entire simulation. The keyword `COMPOSITION` defines which fluid component is being injected. Finally, the temperature of the injected fluid is defined with the `TEMPERATURE` keyword inside the `BOUNDARY_CONDITION` header. 

For more examples on how to define wells in DARSim, please refer to the `Examples` directory. 

---

## 11. Fractures

DARSim supports naturally fractured reservoirs through Embedded Discrete Fracture Models (EDFM) and Projection-based EDFM (pEDFM). 

Fractures are represented as lower-dimensional entities that exchange mass with the surrounding porous matrix through transmissibility-based coupling terms. This approach allows complex fracture networks to be modeled without requiring the computational grid to conform to individual fracture geometries.

To run a simulation without fractures, include the following the main simulation input file:

```
FRACTURED
	OFF
	NA
```

The `NA` is required here, which stands for "not applicable" with regard to the name of the fracture file name. 

When applying pEDFM to a Cartesian or corner-point grid, the user must first call DARSimFracGen.m located inside `src/darsim_legacy/FracturesGenerator`. To call this function, run the following command from within the `FracturesGenerator` directory:

```matlab
DARSimFracGen("<path/to/fracture_input>","<FractureInputFileName.txt>")
```

The file named `FractureInputFileName.txt` (which can have any name) contains the details of the fracture network and paramters defining the pEDFM options. For an example of this file, please refer to the `Examples` directory. 

After `DARSimFracGen` is called, an output file describing matrix and fracture network connectivity is generated, which is called `Fracture_Output.txt`. This file should be located in the same directory as the other simulation files. To include the fractures in the simulation, include the following block in the main simulation input file:

```
FRACTURED
	ON
	Fracture_Output.txt
```

This tells DARSim to apply the fracture connections and where to find the description of the fracture-matrix connectivity (i.e. directs DARSim to use the `Fracture_Output.txt` file). 

When applying pEDFM to unstructured grids, the workflow is somewhat streamlined. Instead of calling `DARSimFracGen` as a pre-processing step to the simulation run, the fracture network and connectivity is generated within the `DARSim` workflow. In this case, the `FRACTRUED` block inside the main simulation input file should read:

```
FRACTURED
	ON
	<FractureInputFileName.txt>
```

where `FractureInputFileName.txt` is the fracture input file that is generally provided to `DARSimFracGen` for Cartesian or corner-point grids. 

#### Note: Cache behavior for Fractures...

For the Cartesian and corner-point grids, `DARSimFracGen` generates the `Fracture_Output.txt` file, which contains all information regarding the connectivity of the fracture network. If changes are made to the fracture properties inside the fracture input file, the user should regeneerate the `Fracture_Output.txt` file. 

For unstructured grids, calling `DARsim` alone is sufficient for running a pEDFM simulation. The fracture network connectivity information is then cached inside a file saved to the simulation directory called `fracture_process.mat`. If changes are made to the fracture properties inside the fracture input file, the suer should first delete `fracture_process.mat` befor running `DARSim`. 

---

## 12. Algebraic Dynamic Multilevel (ADM)

ADM provides dynamic multilevel model reduction to improve computational efficiency while maintaining solution accuracy.

During a simulation, regions requiring high spatial resolution are solved on fine grids, whereas smoother regions are represented on coarser grids. The adaptive refinement strategy reduces computational cost while preserving important flow features near wells and displacement fronts.

**NEEDS MORE INFORMATION ABOUT SETTING UP AN ADM SIMULATION HERE**

---

## 13. Output Controls

Simulation output is controlled through the simulation settings file, which is always names `SimulatorSettings.txt`. Most commonly, the output is defined to be binary VTK format.

```
OUTPUT
	VTK
	BINARY
```

Output frequency may be specified independently of the nonlinear timestep size, allowing detailed post-processing while maintaining efficient time integration. For example, 100 equally spaced output reports is stipulated with:

```
REPORTS
	100
```

The time step reported with the VTK file is the report number. Therefore, to convert the report number to the simulation time, the user must multiply the report number by the output reporting period. The output reporting period is equal to the total simulation time divided by the number of stipulated output reports. 

---

## 14. Visualization

DARSim produces VTK output files that can be visualized using ParaView.

Typical visualization results include:

- Pressure
- Phase saturation
- Component concentration
- Relative permeability
- Capillary pressure
- Temperature
- Fracture variables

Visualization provides an effective means of interpreting reservoir behavior and evaluating simulation results. For more information on using Paraview, please refer to the software documentation.

---

# Additional Resources

Additional information, examples, and developer documentation are provided with the DARSim repository. Users interested in extending or modifying the simulator should refer to the accompanying **Developer Guide**, which describes the internal software architecture, execution workflow, and principal source files.
