# epframe
Elastic-Plastic analysis of 2D structural frames
=======

## Overview

EPFRAME performs incremental elastic-plastic analysis of 2D plane frames using the plastic hinge method. This implementation is translated from the original FORTRAN code by Hacksoo Lee (1986) into Python with modern enhancements.

The program tracks sequential formation of plastic hinges as loads increase, automatically adjusting member stiffnesses until a collapse mechanism forms.

## Features

- **Incremental Load Analysis**: Progressive loading until collapse mechanism forms
- **Plastic Hinge Tracking**: Sequential formation of plastic hinges with load factors
- **Automatic Stiffness Modification**: Member stiffnesses adjust as hinges form
- **Reaction Force Calculations**: Support reactions computed at each load step
- **Visualization Suite**: Automatic generation of deformed shapes, moment diagrams, shear diagrams, and axial force diagrams
- **CSV Data Export**: Compact numerical data exported to CSV for post-processing
- **Comment Support**: Input files can include `#` comments for documentation
- **Modern Python**: Uses numpy for efficient matrix operations

## Files

| File                  | Description                                 |
| --------------------- | ------------------------------------------- |
| `epframe.py`          | Main analysis program                       |
| `epframe_viz.py`      | Visualization and plotting tools            |
| `example_frame_5.dat` | Example input file (7-element portal frame) |

## Installation

### Requirements

```bash
pip install numpy matplotlib
```

### Quick Start

```bash
# Run analysis
python epframe.py example_frame_5.dat results.dat

# Generate visualizations
python epframe_viz.py results.dat
```

## Usage

### Analysis

```bash
python epframe.py input_file output_file
```

**Outputs:**

- `output_file` - Human-readable results with deformations, moments, and reactions
- `output_file.csv` - Compact numerical data for post-processing

### Visualization

```bash
python epframe_viz.py output_file
```

**Generated Plots:**

- `frame_N_geometry.png` - Original frame layout with supports
- `frame_N_deformed_hinge_X.png` - Deformed shape after each hinge
- `frame_N_moments_hinge_X.png` - Bending moment diagrams
- `frame_N_shear_final.png` - Shear force diagram (final state)
- `frame_N_axial_final.png` - Axial force diagram (final state)
- `frame_N_summary.png` - 4-panel progressive collapse summary

## Input File Format

```
Frame_Number
NCT  NE  E                           # Nodes, Elements, Young's modulus

# Node Data (NCT lines)
Node  X  Y  DFX  DFY  DFZ            # Coordinates and DOF flags

# Element Data (NE lines)  
Elem  N1  N2  I  A  Mp               # Connectivity, properties, plastic moment

LN                                   # Number of loaded nodes

# Load Data (LN lines)
Node  FX  FY  MZ                     # Applied forces and moment
```

### Coordinate System

- **X-axis**: Horizontal (→)
- **Y-axis**: Vertical (↑)
- **Z-axis**: Out of plane (rotation follows right-hand rule)

### Degree of Freedom Flags

- **0** = Fixed (restrained)
- **1** = Free (unrestrained)

### Units

Maintain consistent units throughout. Example using US customary:

- Length: inches
- Force: kips
- Moment: in-kips
- Stress: ksi
- Moment of Inertia (I): in⁴
- Area: in²

## Example Input File

```
# EPFRAME Example: 7-Element Portal Frame
# Units: inches, kips, ksi

5                              # Frame number

8   7   29000                  # 8 nodes, 7 elements, E=29000 ksi

# Node data: fixed supports at nodes 1 and 8
1     0     0   0   0   0      # Fixed support (left)
2     0   168   1   1   1      # Column top (left)
3   120   252   1   1   1      # Roof nodes
4   216   252   1   1   1
5   312   252   1   1   1
6   408   252   1   1   1
7   528   168   1   1   1      # Column top (right)
8   528     0   0   0   0      # Fixed support (right)

# Element data: W14x68 sections
1   1   2   954   19.7   4680  # Left column
2   2   3   954   19.7   4680  # Roof elements
3   3   4   954   19.7   4680
4   4   5   954   19.7   4680
5   5   6   954   19.7   4680
6   6   7   954   19.7   4680
7   7   8   954   19.7   4680  # Right column

5                              # 5 loaded nodes

# Loads: lateral + gravity
2   0.50   0.00   0            # Lateral at node 2
3   0.25  -1.00   0            # Combined at node 3
4   0.00  -1.00   0            # Gravity loads
5   0.00  -1.00   0
6   0.00  -1.00   0
```

## Output Format

### Main Output File

```
%
%     ELASTIC PLASTIC ANALYSIS OF FRAME NO   5
%     ---------------------------------------
%
%     * GENERAL DATA
%          NUMBER OF NODES              8
%          NUMBER OF ELEMENTS           7
%          MOD OF ELASTICITY       29000.0
%
%     * DATA FOR NODES
%           NODE   X-COORD   Y-COORD   DFX   DFY   DFZ
%             1        0.00      0.00     0     0     0
%             ...
%
%     * PLASTIC HINGE   1 FORMED IN ELEMENT   7 NEAR NODE   8 WHEN LOAD FACTOR IS       30.795
%
%          CUMULATIVE DEFORMATIONS
%                NODE    X-DISP       Y-DISP       ROTN
%                  1      0.00000      0.00000      0.00000
%                  2     -0.15074     -0.01714      0.00306
%                  ...
%
%          CUMULATIVE MOMENTS
%             ELEMENT       END MOMENTS             NODES     PLASTIC MOM
%                  1       1895.62   2904.68      1 AND 2       4680.00
%                  ...
%
%          REACTIONS AT SUPPORTS
%                NODE       RX           RY           MZ
%                  1        28.57        58.29     -1895.62
%                  8       -51.67        64.89      4680.00
```

### CSV Output File

The CSV file contains one header row followed by data rows for each load step:

```csv
NCYCL,EL,ND,CLG,CX1_X,CX1_Y,CX1_R,...,CM1_1,CM1_2,...,CT1,...
0,0,0,0.000000E+00,0.000000E+00,...
1,7,8,3.079452E+01,-1.507380E-01,...
```

Columns:

- `NCYCL` - Cycle number (0=initial, 1+=after each hinge)
- `EL` - Element where hinge formed
- `ND` - Node where hinge formed
- `CLG` - Cumulative load factor
- `CXn_X/Y/R` - Cumulative displacements at node n
- `CMn_1/2` - Cumulative moments at element n ends
- `CTn` - Cumulative axial force in element n

## Algorithm

1. **Initialization**
   
   - Build compatibility matrix K from geometry
   - Calculate initial member stiffnesses

2. **Load Increment Loop**
   
   - Form global stiffness matrix: K_global = K^T · S · K
   - Solve for displacements: K_global · δ = P
   - Calculate member forces: F = S · K · δ
   - Find load factor λ to first plastic hinge: λ = (Mp - M) / ΔM
   - Update cumulative values (displacements, moments, forces)
   - Modify stiffness at plastic hinge location
   - Repeat until mechanism forms

3. **Plastic Hinge Modification**
   
   - Near-end (hinge location): Stiffness → 0, moment locked at Mp
   - Far-end: Stiffness reduced to 75% (cantilever effect)

4. **Termination Conditions**
   
   - Collapse mechanism forms (singular stiffness matrix)
   - Deformations exceed limit (1000)
   - All elements have formed plastic hinges

## Example Results

For the 7-element portal frame example:

| Hinge | Element | Node | Load Factor |
| ----- | ------- | ---- | ----------- |
| 1     | 7       | 8    | 30.795      |
| 2     | 7       | 7    | 34.400      |
| 3     | 4       | 4    | 38.937      |
| 4     | 2       | 2    | 40.297      |

The frame forms a collapse mechanism after 4 plastic hinges at a load factor of 40.297.

## Visualization Examples

### Progressive Collapse

![Progressive Collapse Summary](frame_5_summary.png)

### Bending Moment Diagram

![Moment Diagram](frame_5_moments_hinge_4.png)

## Limitations

- 2D plane frames only (no 3D analysis)
- Concentrated loads at nodes only (no distributed loads)
- Small displacement theory (no geometric nonlinearity)
- Elastic-perfectly plastic material model
- No member instability (buckling) checks
- No P-Δ effects

## Troubleshooting

**"Singular matrix" or "DIVISION BY ZERO" error:**

- Check for unstable geometry (mechanism before loading)
- Verify support conditions provide adequate restraint
- Ensure all elements are properly connected

**"Deformations exceed 1000":**

- Collapse mechanism has formed (expected behavior)
- If unexpected, check plastic moment capacities

**Unexpected results:**

- Verify consistent units throughout input
- Check element connectivity (N1, N2 assignments)
- 
- Confirm DOF flags (0=fixed, 1=free)
- Review applied load directions and magnitudes

## References

1. Hacksoo Lee, "EPFRAME.F", University of Michigan, 1986
2. Neal, B.G., "The Plastic Methods of Structural Analysis", Chapman & Hall
3. Chen, W.F. and Sohal, I., "Plastic Design and Second-Order Analysis of Steel Frames", Springer

## License

This implementation is based on public domain FORTRAN code
and is now licensed from the MIT License.  

## Author

Original FORTRAN: Hacksoo Lee (1986)  
Python Translation: Duke University Civil & Environmental Engineering

