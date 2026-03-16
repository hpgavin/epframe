# epframe

Elastic-Plastic analysis of 2D structural frames
=======

## Overview

EPFRAME performs incremental elastic-plastic analysis of 2D plane frames using the plastic hinge method. This implementation is translated from the original FORTRAN code by Hacksoo Lee (1986) into Python. 

The program tracks sequential formation of plastic hinges as loads increase, automatically adjusting member stiffnesses until a collapse mechanism forms.

## Features

- **Incremental Load Analysis**: Progressive loading until collapse mechanism forms
- **Plastic Hinge Tracking**: Sequential formation of plastic hinges with load factors
- **Moment-Axial Interaction**: Axial loads reduce ultimate moments
- **Geometric Nonlinearity**: Tension forces increase stability, compression forces decreases stability
- **Unidirectional Reactions**: Reaction forces can be specified to act in only one direction
- **Automatic Stiffness Modification**: Member stiffnesses adjust as hinges form
- **Reaction Force Calculations**: Support reactions computed at each load step
- **Visualization Suite**: Automatic generation of deformed shapes, moment diagrams, shear diagrams, and axial force diagrams
- **Computation of Element Displacements Between Nodes**: Double integrating M(x)/EI between nodal displacements
- **CSV Data Export**: Compact numerical data exported to CSV for post-processing
- **Comment Support**: Input files can include `#` comments for documentation
- **Modern Python**: Uses numpy for efficient matrix operations, scipy for QP, and matplotlib for plotting

## Limitations

- 2D plane frames only (no 3D analysis)
- Concentrated loads at nodes only (no distributed loads)
- Elastic-perfectly plastic material model

## Files

| File                  | Description                                      |
| --------------------- | ------------------------------------------------ |
| `epframe.py`          | Main analysis program                            |
| `epframe_viz.py`      | Visualization and plotting tools                 |
| `square_tube.py`      | Section properties of square tube cross sections |
| `beam_4_example.dat`  | Example input file (4-node beam)                 |
| `beam_7_example.dat`  | Example input file (7-node beam)                 |
| `frame_7_example.dat` | Example input file (7-node gable frame)          |

## Installation

### Requirements

```bash
pip install numpy matplotlib
```

## Usage

### Analysis Step

```bash
# Run analysis
python epframe.py input_file output_file
```

**Outputs:**

- `output_file` - Human-readable results with deformations, moments, and reactions
- `output_file.csv` - Compact numerical data for post-processing

### Visualization Step

```bash
# Generate visualization figures from analysis results 
python epframe_viz.py output_file
```

**Generated Plots in sub-directory ./plots/ :**

- `output_file-geometry.pdf` - Original frame layout with supports
- `output_file-deformed_hinge_X.pdf` - Deformed shape after each hinge
- `output_file-moments_hinge_X.pdf` - Bending moment diagrams
- `output_file-shear_hinge_X.pdf` - Shear force diagram (final state only)
- `output_file-axial_hinge_X.pdf` - Axial force diagram (final state only)
- `output_file-load_displacement.pdf` - load factor vs max displacement 
- `output_file-summary.pdf` - 4-panel progressive collapse summary

## Input File Format

[Example Input File - a gable frame](examples/frame_7_example.txt)

## Output Format

### Main Output File

[Example Output File - a gable frame](examples/frame_7_example.out)

### CSV Output File

[Example Output File - a gable frame](examples/frame_7_example.out.csv)

The CSV file contains one header row followed by data rows for each load step:

Columns:

- `NCYCL` - Cycle number (0=initial, 1+=after each hinge)
- `EL` - Element where hinge formed
- `NH` - Node where hinge formed
- `CLF` - Cumulative load factor (lambda)
- `CDn_X/Y/R` - Cumulative displacements at node n
- `CMn_1/2` - Cumulative moments at element n ends
- `CTn` - Cumulative tension forces in element n

## Example Results

For the 7-element gable frame example:

| Hinge | Element | Node | Load Factor |
| ----- | ------- | ---- | ----------- |
| 1     | 7       | 8    | 18.001      |
| 2     | 6       | 7    | 20.007      |
| 3     | 4       | 4    | 22.334      |
| 4     | 1       | 2    | 22.708      |

## Visualization Examples

### Progressive Collapse

![Progressive Collapse Summary](examples/plots/frame_7_example-summary.png)

### Bending Moment Diagram

![Moment Diagram](examples/plots/frame_7_example-moments_hinge_04.png)

### Load Displacement Diagram

![Load Displacement](examples/plots/frame_7_example-load_displacement.png)

## Algorithm

**1. Initialization**

- Parse title, material properties (E, Fy), section properties (I, A, Z), and support conditions (bidirectional, unidirectional +/−, or free) from input file
- Derive plastic moment Mp = Z·Fy and axial yield force Py = A·Fy for each element
- Build compatibility matrix **K** from frame geometry and element connectivity
- Calculate initial elastic member flexibilities SF (bending) and SA (axial)
- Compute degree of static indeterminacy DI = 3·NE − ND and displacement limit DLMT = 0.1 × max frame dimension
- Write initial state to output and CSV files

**2. Load Increment Loop**

- Form elastic element stiffness matrix **S** from current SF and SA (modified at prior hinge locations)

- Form global elastic stiffness: **K**_e = **K** · **S** · **K**^T

- Assemble geometric stiffness **K**_g by transforming and scattering the 6×6 element geometric stiffness matrices (eq. 75, Bernoulli-Euler) weighted by current cumulative axial forces CT; add to form total stiffness: **K**_sat = **K**_e + **K**_g

- Check for geometric instability: if **K**_sat has a negative eigenvalue, report buckling and stop

- Solve for displacement rates δ using active-set iteration to enforce unidirectional reaction constraints: **K**_sat · δ = **P**, releasing any unidirectional support whose reaction force would act in the wrong direction

- Calculate member force rates: **F** = **S** · **K**^T · δ; separate into moment rates (SATX_e2) and axial force rates (SATX_ct)

- Check for compression yield: if |CT[el]| ≥ Py[el] for any element, report and stop

- Warn if hinge count already exceeds DI (kinematic mechanism masked by geometric stiffness)

- Find load factor α to next P-M hinge by solving the quadratic at each element end:
  
  (Mp·ṗ²/Py²)·α² + (s·ṁ + 2·Mp·P₀·ṗ/Py²)·α + (s·M₀ − Mp(1 − P₀²/Py²)) = 0
  
  where s = sign(M₀), ṁ = moment rate, ṗ = axial force rate, M₀ and P₀ are cumulative values; reduces exactly to the linear formula (Mp − |M₀|)/|ṁ| when ṗ ≈ 0

- Scale all rates by α; update cumulative displacements CD, moments CM, and axial forces CT

- Check cumulative displacements against DLMT; stop if exceeded

- Write hinge location, active support status, deformations, moments (with Mp, Mu = Mp(1−(P/Py)²), and P/Py), axial forces, and reactions to output file

**3. Plastic Hinge Modification**

- At the hinge end (near-end): set both SF entries to zero — bending stiffness goes to zero, moment is locked at current CM value
- At the far end of the same element: reduce the diagonal SF entry to 75% of its current value (4EI/L → 3EI/L for an intact element, reflecting the loss of far-end rotational restraint) and set the off-diagonal SF entry to zero

**4. Termination Conditions** (in order of priority)

- **Geometric instability (buckling):** minimum eigenvalue of **K**_sat < 0 before a hinge is found
- **Collapse mechanism:** **K**_sat is singular (solve returns None); distinguished from buckling by the sign of the minimum eigenvalue — zero eigenvalue is a mechanism, negative is buckling
- **Compression yield:** |P| ≥ Py in any element, signalling material instability in axial compression
- **Displacement limit:** any cumulative nodal displacement exceeds DLMT = 0.1 × max(frame span, frame height) — enforced on CD after each update, not on the unscaled rate vector
- **No further hinge possible:** all α values exceed 10⁹ (all element ends are either already at Mu or relieving)
- **Kinematic mechanism warning:** hinge count exceeds DI; analysis continues only if geometric stiffness keeps **K**_sat non-singular, and the warning is written to the output file at each such increment

## Troubleshooting

**"Singular matrix" or "DIVISION BY ZERO" error:**

- Check for unstable geometry (mechanism before loading)
- Verify support conditions provide adequate restraint
- Ensure all elements are properly connected

**"Deformations exceed Limit":**

- Collapse mechanism has formed (expected behavior)
- If unexpected, check plastic moment capacities

**Unexpected results:**

- Verify consistent units throughout input
- Check element connectivity (N1, N2 assignments)from
- Confirm DOF flags (0=fixed, 1=free)
- Review applied load directions and magnitudes

## References

1. Hacksoo Lee and Subhash Goel, "EPFRAME.F", University of Michigan, 1986
2. Neal, B.G., "The Plastic Methods of Structural Analysis", Chapman & Hall
3. Chen, W.F. and Sohal, I., "Plastic Design and Second-Order Analysis of Steel Frames", Springer

## License

This implementation is based on public domain FORTRAN code and is now licensed under the MIT License.  

## Author

Original FORTRAN: Hacksoo Lee and Subhash Goel, University of Michigan (1986) 
Python Translation: Duke University Civil & Environmental Engineering
