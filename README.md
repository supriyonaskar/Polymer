# Polymer
A Python-based interface for generating customizable polymer structures. This tool allows you to define monomer units, polymerization parameters, and chain architectures, facilitating seamless integration with molecular simulation workflows. Ideal for creating complex polymer systems for computational material science.

# Python Script for Optimizing Monomers and Building Polymers

This script performs the following tasks:
1. Optimizes a monomer's geometry using Psi4.
2. Calculates ESP and RESP charges.
3. Builds a polymer using Pysimm.
4. Please install Pysimm_by_SN from my repository and install it to run this code and add the below lines to your bashrc file.
   export LAMMPS_EXEC=path_to_your/lmp_mpi
   export PYTHONPATH=$PYTHONPATH:/path_to_your/pysimm
   export PATH=$PATH:/path_to_your/bin
   export ANTECHAMBER_EXEC=~/path_to_your/antechamber



## Prerequisites

### Required Python Modules
- `numpy`
- `math`
- `mdtraj`
- `pysimm`
- `psi4`
- `resp`
- `pandas`

### Required Files
- `input_parameters.dat`: This file should have 7 columns:
  1. Index of the monomer's head atom (index starts from 1).
  2. Index of the monomer's tail atom.
  3. Hydrogen attached to the head atom of the monomer.
  4. Hydrogen attached to the tail atom of the monomer.
  5. Name of the molecule (Ensure the initial XYZ file of the monomer has the name format `{molecule_name}_initial.xyz` and contains two remark lines at the beginning).
  6. Flag to keep (0) or remove (1) unnecessary files during execution.
  7. Optimization flag: 
      - 1: Optimize and prepare the polymer.
      - 0: Only build the polymer without optimizing.

### Input Structure
- `{molecule_name}_initial.xyz`: Initial monomer structure (from Avogadro or another tool).

## Key Steps

1. **Geometry Optimization:**
   - Uses Psi4 with the `scf/cc-pVDZ` method.
   - Produces optimized geometry in `{molecule_name}_optimized.xyz`.

2. **RESP Charge Calculation:**
   - Conducts ESP and RESP charge fitting.
   - Outputs charges in the `.mol2` file format for Pysimm.

3. **Polymer Building:**
   - Builds a polymer using Pysimm's random walk algorithm.
   - Caps polymers with hydrogen atoms if necessary.
   - Minimizes the polymer structure.
   - Outputs polymer configurations in `.lmps` format.

## Usage

### Execution Steps
1. Prepare the `input_parameters.dat` file and `{molecule_name}_initial.xyz` file.
2. Run the script:
   ```bash
   python script_name.py
   ```

### Example: Input Parameters
```
1  # Head atom index
2  # Tail atom index
3  # Hydrogen attached to head atom
4  # Hydrogen attached to tail atom
molecule_name  # Name of the molecule
1  # Remove unnecessary files (1: Yes, 0: No)
1  # Optimize and build polymer (1: Yes, 0: No)
```

## Notes
- The script modifies the `gaff.json` file for force-field parameters if required.
- Ensure sufficient memory and threads for Psi4:
  ```python
  psi4.set_memory('100 GB')
  psi4.set_num_threads(20)
  ```

### Troubleshooting
- If force-field parameters are missing, update `gaff.json` accordingly.
- Ensure the input XYZ file is properly formatted.

---

### Written by:
**Supriyo Naskar**
Institute: ICGM, CNRS, ENSCM, University of Montpellier
