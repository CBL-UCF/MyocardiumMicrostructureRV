# Sheetlet-Based Myocardial Microstructure Generation Toolbox

This MATLAB toolbox generates synthetic sheetlet-based myocardial microstructures for representative unit cell (RUC) models. The code creates sheetlet arrangements, generates cardiomyocyte-like fiber networks, adds intra- and inter-sheetlet branching, and prepares voxelized masks and geometry data for downstream simulation or analysis.

## Main Features

- Generate sheetlet arrangements using Voronoi/Thiessen tessellation
- Shrink sheetlet polygons to define cleavage-plane regions
- Generate cardiomyocyte base points inside sheetlets
- Extrude cardiomyocyte chains along the fiber direction
- Add branching between cardiomyocyte chains
- Build global and sheetlet-wise connectivity tables
- Convert geometric cell data into voxelized masks
- Visualize sheetlet boundaries and microstructural networks
- Save connectivity and geometry data for downstream mechanical simulations and mask generation

## Toolbox Structure

| File | Description |
|---|---|
| `Sheetlet_arrangement_Toolbox.m` | Generates sheetlet seed points, Voronoi polygons, and shrunk sheetlet regions |
| `Sheetlet_microstructure_Toolbox.m` | Generates intra-sheetlet cardiomyocyte networks and branching |
| `Coordinate_Toolbox.m` | Provides coordinate conversion, translation, rotation, and distance-calculation utilities |
| `Mask_generate_Toolbox.m` | Converts cell geometry and sheetlet regions into voxel-level mask information |
| `Mask_reshape_Toolbox.m` | Provides mask post-processing and slice perturbation tools |
| `Plotting_Toolbox.m` | Visualizes sheetlet boundaries and microstructural connectivity |
| Main script | Defines model parameters, generates connectivity, applies transformations, and prepares mask generation |

## Requirements

This code requires MATLAB with support for:

- `voronoin`
- `inpolygon`
- `graph`
- `parfor` for parallel collision detection

Recommended toolboxes:

- MATLAB Parallel Computing Toolbox
- MATLAB Statistics and Machine Learning Toolbox

## Basic Workflow

1. Define the representative unit cell dimensions.
2. Define sheetlet width, thickness, and shrink ratios.
3. Generate sheetlet arrangements using Voronoi tessellation.
4. Generate cardiomyocyte chains within each sheetlet.
5. Add branches to approach the target intra-cellular volume fraction.
6. Translate the generated structure to the desired coordinate origin.
7. Generate cell detail lists and voxel masks.
8. Save connectivity and geometry data for later use.

## Example Parameters

```matlab
ModelName = 'Microstructure_test1';

RUC_w = 200;
RUC_t = 80;
RUC_l = 500;

shtlt_w = 200;
shtlt_t = 45;

shtlt_shrink_ratios = [0.955, 0.825];

innitial_myocyte_diameter = 10;
myocyte_length_range = [100, 120];

shtlt_ICV = 0.95;
appx_target_ICV = 1.36;

branching_reg_thickness = 40;
randomization_breakpoint = 500000;
```

Running the Code

Run the main script after placing all toolbox files in the MATLAB path.

```matlab

run main_script_name.m
```

To generate a new microstructure, set:

```matlab
Construct_connectivity = 1;
```

To load a previously saved microstructure, set:
```matlab
Construct_connectivity = 0;
```

Saved Output Files

The workflow saves intermediate and final data files using the selected ModelName

Example outputs:

```matlab
Micro_test1_conData.mat
Micro_test1_mechInput.mat
```

The connectivity file stores generated microstructural connectivity, sheetlet regions, and branch information.

The mechanics input file stores translated coordinates and sheetlet information for downstream modeling.

Units

Most geometric inputs are defined in micrometers. Some cell detail lists convert coordinates, lengths, and radii to meters for volume calculations.

Notes
Cardiomyocytes are approximated as cylindrical segments.
Approximate branching volume fraction may exceed one because overlapping cylindrical nodes are not corrected at this stage.
The current implementation is intended for synthetic microstructure generation and require further validation before direct biological interpretation.

Citation


Author
