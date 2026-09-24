
# DARSim Developer Guide

![DARSim Logo](../Logo/DARSimLogo2.pdf)

## Table of Contents

- [Introduction](#1-introduction)
- [DARSim Workflow](#2-darsim-workflow)
- [Key Source Files](#3-key-source-files)
- [Git Setup](#4-git-setup)
- [Best Practices for Development](#5-best-practices-for-development)
- [Example Development Workflow](#6-example-development-workflow)
- [Integration Tests](#7-integration-tests)

---

## 1. Introduction

This document describes a general overview of the DARSim simulation workflow. Additionally, the code architecture is summarized and key source code files are identified and described. Finally, the version control workflow (Git) used in the DARSim development environment hosted on GitLab is presented. It provides definitions of basic Git concepts and outlines the standard development procedure for contributors. 

---

## 2. DARSim Workflow

This section summarizes the principal execution workflow within DARSim. Unlike the User Guide, the emphasis here is on how the software is assembled and how the principal modules interact during a simulation.

### 2.1 Overall Execution

A DARSim simulation is initiated from the MATLAB `src` directory by calling:

```matlab
DARSim("<path/to/simulation>","<InputFileName.txt>")
```

The overall execution sequence is:

```text
DARSim()
 ├── startup_darsim
 ├── create_darsim_settings()
 │    ├── reader_darsim3()
 │    ├── blueprint_factory()
 │    └── darsim2_settings_factory()
 ├── create_simulator()
 │    ├── simulation_factory()
 │    └── BuildSimulation()
 ├── prepare_simulation()
 └── Run()
```

The entry point, `DARSim.m`, initializes the MATLAB environment, reads the simulation inputs, constructs the simulation settings, creates the simulator object, prepares the simulation, and starts the solver for the transient solution.

### 2.2 Simulation Assembly

The simulation is assembled primarily through `simulation_factory.m`, which constructs the major simulation objects, including:

- Discretization model
- Fluid model
- Production system
- Formulation
- Time driver
- Writer
- Summary
- Initializer

These components are combined to create a `SimulationDARSim` object. Before time stepping begins, the simulator initializes the computational grid, state variables, wells, mappings, and data structures required for the reservoir and fracture unknowns.

### 2.3 Phase Behaviour and Flash Calculations

For compositional, black-oil, geothermal, and thermal simulations, DARSim updates the thermodynamic state after each nonlinear state update.

Depending on the selected fluid model, this process may include:

- Flash calculations;
- Phase composition updates;
- Density and viscosity calculations;
- Saturation updates;
- Capillary pressure calculations;
- Thermal property calculations for non-isothermal simulations.

These property evaluations are fully coupled with the nonlinear reservoir solution.

### 2.4 Time Stepping and Nonlinear Solution

DARSim employs a fully implicit solution strategy.

For each timestep, the simulator:

1. Selects the current timestep size;
2. Updates fluid properties;
3. Assembles the nonlinear residual and Jacobian;
4. Performs Newton iterations;
5. Solves the linearized system;
6. Checks convergence;
7. Updates the timestep size if necessary;
8. Writes simulation output.

This workflow provides robust performance for strongly coupled multiphase and multicomponent flow problems.

---

## 3. Key Source Files

The following source files provide the primary entry points into the DARSim codebase. They are grouped according to their role within the simulator.

Additionally, there is a code organization map in HTML format provided with the source code. To view, please open `code-organization.html` found inside the `/src` directory. Here, you can navigate and search through the various source code files and get an overview of the code base structure. 

### 3.1 Entry Points and Simulation Assembly

These files control the overall execution of a DARSim simulation.

| File | Purpose |
|------|---------|
| `src/DARSim.m` | Main entry point for all simulations. |
| `src/reader/create_darsim_settings.m` | Creates the simulation settings from the user input files. |
| `src/reader/reader_darsim3.m` | Reads the primary simulation input files. |
| `src/reader/blueprint-factory/blueprint_factory.m` | Organizes user input into blueprint objects. |
| `src/reader/darsim2_settings_factory.m` | Translates blueprint objects into the legacy solver input structure. |
| `src/simulation/create_simulator.m` | Creates the simulator object. |
| `src/simulation/SimulationDARSim.m` | Main simulation class responsible for initialization and execution. |
| `src/factory/simulation_factory.m` | Constructs the principal simulation components. |

### 3.2 Grid Construction

These files construct the computational mesh and associated data structures.

| File | Purpose |
|------|---------|
| `src/factory/mesh-factories/rock-mesh-factory/build_mesh.m` | Main mesh construction routine. |
| `build_cornerpoint.m` | Corner-point grid construction. |
| `manage_extension.m` | Handles external mesh file formats. |
| `MeshCore.m` | Core unstructured mesh representation. |
| `translateMeshCore.m` | Converts `MeshCore` objects into DARSim grid structures. |

### 3.3 Mathematical Formulation and Solver

The nonlinear reservoir equations and time integration are implemented primarily within the following files.

| File | Purpose |
|------|---------|
| `Compositional_formulation.m` | Fully compositional formulation. |
| `Overall_Composition_formulation.m` | Overall-composition formulation. |
| `Immiscible_formulation.m` | Immiscible flow formulation. |
| `TimeLoop_Driver.m` | Fully implicit timestep driver and nonlinear solution loop. |

### 3.4 Fracture Implementation

The fracture workflow is implemented through dedicated reader, preprocessing, and mesh construction routines.

| File | Purpose |
|------|---------|
| `legacy_fracture_reader.m` | Reads legacy fracture descriptions. |
| `legacy_scanfracture_data.m` | Processes legacy fracture connectivity. |
| `scanfracture_data_uns.m` | Processes fractures for unstructured grids. |
| `FractureFactory.m` | Constructs fracture objects. |
| `FractureGeometryBuilder.m` | Generates fracture geometry. |

### 3.5 ADM Implementation

** WE NEED MORE DETAILED INFORMATION ON THE CODE ARCHITECTURE/WORKFLOW FOR ADM **

---

## 4. Git Setup

This section preserves the Git workflow adopted for DARSim development.

### 4.1 First-Time Setup

1. Register for a GitLab account.
2. Ask a DARSim code maintainer to add you as a team member to the GitLab project.
3. Install Git.
   - For Mac/Linux, Git is typically pre-installed and available through the Terminal.
   - For Windows, install Git and use Git Bash to execute Git commands.
4. Generate an SSH key for authentication (see the GitLab documentation).
5. Navigate to the location where you would like to clone the DARSim repository with the following command (*It is often a good idea to make a directory called /Codes somewhere safe.*):

```bash
cd /path/to/Codes
```

6. Clone the repository:

```bash
git clone git@gitlab.com:darsim2.0/darsim.git
```

7. Verify that the `main` branch is active with the following command (*You should see `* main` printed on the screen.*):

```bash
git branch
```

8. Open MATLAB.
   - Change the working directory to `src`.
   - Add `src` and all subfolders to the MATLAB path by selecting "Set Paths" from the toolbar.
9. Run an example simulation:

```matlab
DARSim("Examples/<example name>", "Sim.txt")
```

Congratulations! You now have a local Git-enabled repository of DARSim. 

### 4.2 Basic Git Concepts
Below is a list of basic Git principles that are important to understand for any DARSim contributor. This is not intended to be a comprehensive review of Git concepts. For a deeper grasp of the Git principles, please refer to one of the many online resources that are freely available. 

### Concept 1: Remote Repository

Central project containing code, documentation, and history hosted on GitLab.

### Concept 2: Local Repository

A clone of the remote repository hosted on your local computer, where edits can be made. The local repository can also be “synced” with the remote one, if desired. 

### Concept 3: Branches

Independent versions (or copies) of the code, used for development purposes. 

- The `main` branch is the default, stable version of DARSim and is protected (meaning it can only be edited through a “merge” after request and review by a code maintainer). 
- `feature` branches (with a user-given name) can be created to perform development work by developers/contributors. These branches will later be merged with the main branch after a testing and review. 

### Concept 4: Checkout

Command for switching between different branches on your local repository. The files in your repository will reflect only the active branch. You can check your active branch with 

```bash
git branch
```

and switch to a different branch with:

```bash
git checkout <branch_name>
```

### Concept 5: Commit

A saved snapshot of changes made to a feature branch in your local repository, always with a descriptive message. You can restore your local feature branch to any of your commits at a later time, so these are helpful when your code starts to break and you just want to reset to some previous point.  

Typical workflow:

```bash
git add .
git commit -m "Descriptive commit message"
```

### Concept 6: Merge

Combines two branches. This is usually performed via a "Merge Request" made on GitLab, but can also be done locally between two local branches. To merge with `main`, a "Merge Request" must be created on GitLab, after performing the required testing and a code review from a code maintainer. To merge a specified `feature` branch into your currently active branch, use the command:

```bash
git merge <feature_branch>
```

### Concept 7: Push

Sends your local commits made to a feature branch on your computer back to the remote repository on GitLab. This makes the feature branch of the remote repository match the feature branch in the local repository. If this is the first push, then the feature branch will be created on the remote repository automatically. Note: If conflicts exist between the local and remote branch, GitLab will ask to resolve these manually. To perform a push:

```bash
git push origin <feature_branch>
```

### Concept 8: Fetch

Retrieve updates of a branch from the remote repository without modifying any of your local branches. A **fetch** downloads the latest commits from the specified remote branch and updates Git's local reference to that remote branch (for example, origin/feature). Your local `feature` branch, if one exists, is not updated or merged automatically, and your currently activated branch and working files remain unchanged. This makes **fetch** useful for retrieving and inspecting changes from the remote repository before deciding whether to merge them into a local branch. A **fetch** is also useful when a branch exists on the remote repository but not yet in your local repository: you can first **fetch** the remote `feature` branch and then perform a checkout to create a corresponding local `feature` branch. To perform a fetch:

```bash
git fetch origin <branch_name>
```

### Concept 9: Pull

Retrieve updates of a branch from the remote repository and merges them with the current active branch on your local repository. Notice that this means whichever branch is pulled from the remote will be fully merged into whichever branch is currently activated locally. If you **pull** a `feature` branch while you have the `main` branch activated locally, the `feature` branch will be merged into you local `main` branch; in most cases this is an issue and not the desired outcome. Therefore, it is typically best to first **checkout** the branch that you are about to **pull**. Note: If conflicts exist between the local and remote branch, GitLab will ask to resolve these manually. To perform a pull:

```bash
git pull origin <branch_name>
```

### Concept 10: Roles

- **Maintainer:** Reviews and merges contributions with the `main` branch.
- **Developer/Contributor:** Develops and tests new functionality on feature branches.

---

## 5. Best Practices for Development

- Commit frequently using *descriptive messages*.
- Develop on feature branches rather than directly on `main`.
- Always pull the latest version of the `main` branch before creating a new `feature` branch.
- Avoid committing simulation output files or other large generated data.
- Ensure that all integration tests pass before requesting a merge with `main`.

---

## 6. Example Development Workflow

Here is an example workflow of DARSim development using Git. Suppose you need to fix a bug in the relative permeability model files. Here is a description of the workflow:

### Update the local `main` branch.

You first want to make sure that your local “main” branch matches the remote branch exactly before starting your development work, so you perform a pull. 

Navigate to your local repository: 

```bash
cd </your/local/repo/path>
```

Make sure you have the correct branch active:

```bash
git branch
```

If needed, switch the branch to `main`:

```bash
git checkout main
```

**Pull** from the remote repository to update the local version of `main`:

```bash
git pull origin main
```

### Create `feature` branch

Next, you want to create a new feature branch based on `main`. 

With the `main` branch active, create a new `feature` branch with **Checkout**:

```bash
git checkout -b RelPerm_bug_fix
```

### Edit source code

With the new `feature` branch active (called `RelPerm_bug_fix`), you can now start editing the source code as usual. You can run DARSim with these edits from MATLAB as well. 

### Save progress

Once you’ve reached a good stopping point (e.g. you fixed one problem, but now a new problem has come up), you should save your progress with a **Commit**.

Check what changes have been made since the last commit: 

```bash
git status
```

Add all the changes to the next commit:

```bash
git add .
```

Make your **Commit** with a *descriptive message*:

```bash
git commit -m "Fixed typo in RelPerm.m on line 106"
```

Repeat `git status` and confirm that you see "nothing to commit, working tree clean". This confirms that you **Commit** worked correctly. 

### Let a colleague contribute as well

Your colleague wants to help with the development, so you can share your up-to-date `feature` branch with them over GitLab using a **Push**.

Send your new feature branch to GitLab: 

```bash
git push origin RelPerm_bug_fix
```

Now, go to GitLab to make sure the new feature branch appears. You collegauge can now grab the feature branch with a **fetch**:

```bash
git fetch origin RelPerm_bug_fix
```

Then, you colleague can **checkout** this branch locally to activate it:

```bash
git checkout RelPerm_bug_fix
```

Your colleague can make edits on their computer and push them back to GitLab when they are done. 

### Grab your colleagues edits

Now you need to grab your colleagues edits before finishing your work. Do this with another **pull**, but this time for the `feature` branch (and not the `main` branch). 

```bash
git checkout RelPerm_bug_fix
git pull origin RelPerm_bug_fix
```

Reminder: If you don't first checkout the `feature` branch, the **pull** command will merge the remote `feature` branch code into whichever branch you have currently active. 

### Finish your edits

From here, you can finish make the required edits to the source code files to wrap up the development work. 

### Save your work

When you've finished your edits, save your work with a **commit**. 

```bash
git add .
git commit -m "Updated Brooks Corey model, now has hysterisis"
```

You can check that your commit worked properly with a status check

```bash
git status
```

You should see "nothing to commit, working tree clean" if your commit worked properly. 

### Run the integration tests

The directions below provide details on running the integration tests. These make sure that nothing else in the code was broken by the changes you and your colleague made. If all tests pass without issue, you can proceed to the merging step. 

### Send your new edits to the remote repository

Just as before, send your updated `feature` branch to the remote repository with a **push** command:

```bash
git push origin RelPerm_bug_fix
```

### Merge request with main

Finally, make a “Merge Request” on GitLab, telling the code maintainers that the new `feature` branch is ready to be included in the stable `main` branch. Follow these steps:

1.	From the remote repo on GitLab.com, select the feature branch from the drop-down menu. 
2.	Click “create merge request” and fill in the form with a description for the merge request, assign a code maintainer to the request, and request review from a team member if desired. Select the “delete branch after merge” option, which will automatically delete RelPerm_bug_fix branch after it is merged. 
3.	A code maintainer can finalize the merge of RelPerm_bug_fix with main branch. 

### After merging, update your local `main` branch to reflect the updates

Once the merge is complete, **pull** the new `main` branch to your local repo to reflect the changes. This must include the `--prune` option, which removes the branch tracking for the `feature` branch. 

```bash
git checkout main
git pull origin main --prune
```

### Delete your local `feature` branch

You no longer need the `feature` branch locally now that it is fully merged with the `main` branch. To delete it:

```bash
git branch -d RelPerm_bug_fix
```

If this does not work, it is likely because they merge with `main` was performed as a “squash” merge (meaning all commits in the feature branch were squashed into a single commit in the main branch). In this case, you need to force the delete with:

```bash
git branch -D RelPerm_bug_fix
```

---

## 7. Integration Tests

DARSim includes a suite of lightweight integration tests intended to verify that the simulator's core functionality remains operational after development changes. These test are summarized in Table 1 below. 

Note: These tests *do not* test all of DARSim's capabilities and *do not* gaurantee that a new feature branch of DARSim is bug-free. They simply provide a quick assessment of DARSim's basic functionality to provide at least *some* confidence in a new feature branch performance. 

To execute the integration tests:

1. Navigate to the `IntegrationTests` directory in MATLAB.
2. Run:

```matlab
RunIntegrationTests(<LogFileName.txt>)
```

where `<LogFileNAme.txt>` is an arbitrary file name for storing the integration test results. By default, this log file is named `IntegrationTestLog.txt`.

After execution, review the generated log file. The summary reports:

- Total number of tests;
- Number of skipped tests (with explanations);
- Number of failed tests (with explanations);
- Number of passed tests.

If a test fails, this is either due to a bug in the new feature branch you are testing, or because the benchmark output of the integration test should be updated (i.e. the new feature branch is fixing a bug in the main branch, and we expect the benchmark output to change.)

The integration tests provide a rapid assessment of software integrity but do not replace comprehensive validation for new physical models or numerical methods.

**Table 1. Review of the integration tests included with DARSim.**

| Test | Grid     | Fluid model       | ADM | pEDFM | Heterogeneity | Anisotropy |
|-----:|----------|-------------------|:---:|:-----:|:-------------:|:----------:|
| 1  | Cart     | CO2 black oil      | – | – | – | – |
| 2  | CPG      | CO2 black oil      | – | – | – | – |
| 3  | Unstruct | CO2 black oil      | – | – | – | – |
| 4  | Cart     | H2 black oil       | – | – | – | – |
| 5  | CPG      | H2 black oil       | – | – | – | – |
| 6  | Unstruct | H2 black oil       | – | – | – | – |
| 7  | Cart     | CO2 comp           | – | – | – | – |
| 8  | CPG      | CO2 comp           | – | – | – | – |
| 9  | Unstruct | CO2 comp           | – | – | – | – |
| 10 | Cart     | H2-brine comp      | – | – | – | – |
| 11 | CPG      | H2-brine comp      | – | – | – | – |
| 12 | Unstruct | H2-brine comp      | – | – | – | – |
| 13 | Cart     | H2-CH4-brine comp  | – | – | – | – |
| 14 | CPG      | H2-CH4-brine comp  | – | – | – | – |
| 15 | Unstruct | H2-CH4-brine comp  | – | – | – | – |
| 16 | Cart     | CO2 black oil      | ✓ | – | – | – |
| 17 | Cart     | H2 black oil       | ✓ | – | – | – |
| 18 | Cart     | CO2 comp           | ✓ | – | – | – |
| 19 | Cart     | H2-brine comp      | ✓ | – | – | – |
| 20 | Cart     | CO2 black oil      | – | ✓ | – | – |
| 21 | CPG      | CO2 black oil      | – | ✓ | – | – |
| 22 | Unstruct | CO2 black oil      | – | ✓ | – | – |
| 23 | Cart     | CO2 black oil      | – | – | ✓ | – |
| 24 | CPG      | CO2 black oil      | – | – | ✓ | – |
| 25 | Unstruct | CO2 black oil      | – | – | ✓ | – |
| 26 | Cart     | CO2 black oil      | ✓ | – | ✓ | – |
| 27 | Cart     | CO2 black oil      | – | – | – | ✓ |
| 28 | CPG      | CO2 black oil      | – | – | – | ✓ |
| 29 | Unstruct | CO2 black oil      | – | – | – | ✓ |

---
