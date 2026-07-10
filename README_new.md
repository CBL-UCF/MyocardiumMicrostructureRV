# Myocardium Microstructure RV

A MATLAB toolbox for generating geometry-based myocardial microstructure models with realistic micro- and mesoscale organization. The framework generates sheetlet-based representative unit cell (RUC) models, constructs cardiomyocyte connectivity with branching, voxelizes the generated geometry, and performs mask reshaping to create realistic cellular architectures for downstream diffusion and mechanical simulations.

## Overview

The toolbox provides an end-to-end pipeline for:

- Generating myocardial sheetlet organization
- Creating cardiomyocyte networks inside sheetlets
- Adding intra- and inter-sheetlet branching
- Generating volumetric voxel-level intracellular/extracellular masks
- Improving voxelized geometry using mask reshaping
- Extracting regions of interest (ROIs) for diffusion simulations
- Exporting data for visualization and computational modeling

## Workflow

```text
Sheetlet Arrangement
          |
          v
Generate Intra-Sheetlet Cardiomyocytes
          |
          v
Add Inter-Sheetlet Branches
          |
          v
Generate Global Connectivity Network
          |
          v
Voxelize Microstructure
          |
          v
Generate ICV / ECV Mask
          |
          v
Mask Reshaping
          |
          v
Extract Simulation ROI
          |
          v
Diffusion / Mechanics Simulation Input
```

Toolbox Structure
|File |	Description |
|---|---|
|`RUC_microstructure_Toolbox.m`|	`Main toolbox for generating complete representative unit cell myocardial microstructures`|
|`Sheetlet_arrangement_Toolbox.m`|	`Generates sheetlet organization, Voronoi regions, and cleavage planes`|
|`Sheetlet_microstructure_Toolbox.m`|	`Generates cardiomyocyte chains, intra-sheetlet networks, and branching`|
|`Coordinate_Toolbox.m`|	`Coordinate transformations, geometry operations, and distance calculations`|
|`Mask_generation_Toolbox.m`|	`Converts geometric microstructure into voxelized masks`|
|`Mask_reshape_Toolbox.m`|	`Improves voxelized masks by enforcing cellular continuity and boundaries`|
|`Plotting_Toolbox.m`|	`Visualization functions for sheetlets and microstructure networks`|

1. Microstructure Generation

The first stage generates the geometrically defined myocardial RUC.

Define RUC parameters
```matlab
ModelName = 'Branched_myocardium_healthy';

RUC_w = 500;
RUC_t = 500;
RUC_l = 500;

shtlt_w = 200;
shtlt_t = 45;

shtlt_shrink_ratios = [0.955, 0.825];

innitial_myocyte_diameter = 10;
myocyte_length_range = [110,130];

shtlt_ICV = 0.95;
appx_target_ICV = 1.36;

branching_reg_thickness = 40;
Generate connectivity
[RUC_con_table_branched,...
 RUC_v_coord,...
 Final_branches,...
 RUC_con_table_sheetlet_wise] = ...
 RUC_microstructure_Toolbox.generate_ruc_myoconnectivity(...);
```
The generated connectivity represents:

Cardiomyocyte centerlines
Sheetlet-wise organization
Intra-sheetlet branches
Inter-sheetlet branches


2. Coordinate Processing

The generated structure is translated to a defined coordinate origin.
```matlab
min_point = [0,0,0];

coord_details = Coordinate_Toolbox.Modification_amount_and_dimension(...);

Modified_point_list = ...
Coordinate_Toolbox.Modify_RUC_coord_mask(...);
```
Outputs:

Modified node coordinates
Translated sheetlet regions
Geometry ready for voxelization


3. Voxel Mask Generation

Mask_generation_Toolbox converts geometric cardiomyocyte structures into voxel-level representations.

The generated mask contains:

Intracellular volume (ICV)
Extracellular volume (ECV)
Cell identity labels
Mask representation
Binary mask:

ICV voxel = 1
ECV voxel = 0


Cell-ID mask:

Cell voxel = cell ID
ECV voxel  = 0
Generate cell geometry
```matlab
shtlt_wise_cell_list = ...
Mask_generation_Toolbox.Sheetletwise_cell_details(...);

inter_shtlt_branch_list = ...
Mask_generation_Toolbox.gen_cell_detail_list(...);
```
Generate Intra-Sheetlet Mask

Each sheetlet is processed independently.
```matlab
[cell_mask,...
 cell_mask_id,...
 sheetlet_wise_voxel_cell_id,...
 n_cells] = ...
Mask_generation_Toolbox.update_inside_sheetlets(...);
```
This step:

Converts voxel indices into physical coordinates
Detects voxel-cylinder intersections
Assigns voxels to cardiomyocytes
Generates cell-ID mapping
Generate Inter-Sheetlet Branch Mask

Candidate voxels in cleavage planes are identified:
```matlab
cleavage_voxels = ...
Mask_generation_Toolbox.filter_unassigned_voxels(...);
```
Branches are added:
```matlab
[mask_before_reshape,...
 mask_id_before_reshape] = ...
Mask_generation_Toolbox.update_inter_sheetlet_branches(...);
```
4. Mask Reshaping

Voxelization can introduce artificial discontinuities.
Mask_reshape_Toolbox improves cellular geometry after mask generation.

Reshaping workflow
Initial Mask
     |
     v
Slice Perturbation
     |
     v
Recover Missing Cell Voxels
     |
     v
Improve Cell Continuity
     |
     v
Define Cell Boundaries
     |
     v
Final Reshaped Mask
Slice-Based Reshaping
```matlab
reshaped_mask_id = ...
Mask_reshape_Toolbox.reshape_mask_shake_slices(...);
```
This step:

Applies controlled slice perturbations
Improves slice-to-slice continuity
Reduces voxelization artifacts
Define Cell Boundaries
```matlab
final_mask_id = ...
Mask_reshape_Toolbox.reshape_mask_define_cellboundary(...);
```
This step:

Detects neighboring cell IDs
Checks true connectivity
Separates unrelated adjacent cells
Preserves connected branches

5. ROI Extraction for Diffusion Simulation

A smaller simulation domain can be extracted from the generated model domain:
```matlab
domain_min_point = [60,100,100];
domain_dim = [500,500,500];

mask_ROI = Diff_domain_mask(cell_mask,...
                            domain_min_point,...
                            domain_dim);
```
Volume fraction:
```matlab

vol_frac = nnz(mask_ROI)/numel(mask_ROI);
```
Visualization

Generated structures can be visualized using:
```matlab
Plotting_Toolbox.plot_microstructure(...)
```

or exported to external visualization tools.

Requirements

Tested with MATLAB.

Recommended:

Parallel Computing Toolbox
Statistics and Machine Learning Toolbox

Main MATLAB functions used:

voronoin
delaunay
graph
inpolygon
parfor
Output Files

Typical generated files:

*_conData.mat
    Connectivity information
    Branch information
    Sheetlet geometry


*_mechInput.mat
    Coordinates
    Sheetlet regions
    Mechanical simulation inputs


Units
Most geometric inputs are defined in micrometers. Some cell detail lists convert coordinates, lengths, and radii to meters for volume calculations.

Notes
Cardiomyocytes are represented using connected cylindrical segments.
Sheetlet separation represents cleavage-plane organization.
Branches provide intra- and inter-sheetlet connectivity.
Final voxel masks preserve realistic intracellular architecture.