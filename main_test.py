import os
from functions_test import * # Import the updated functions

workDir = os.getcwd()
#POTENTIALLY NOT WORKING? FIX ME
#check_gromacs_availability()
file = select_pdb(workDir)
if file:
    top_choice = select_topology(workDir)
    if top_choice:
        top_file = (
            top_choice[1]
            if top_choice[0] == 'Auto'
            else top_choice[1]
        )  # Get the topology file name
        ion_choice = select_ions()
        
        # Run the GROMACS pipeline
        protein_gro_file, top_file = pdb2gmx(file, top_file)
        solvated_gro_file, top_file = solvate(protein_gro_file, top_file)
        ions_tpr, top_file = pgenion(solvated_gro_file, top_file)
        ionized_gro_file, top_file = genion(ions_tpr, top_file, ion_choice)
        em_gro_file, top_file = energy_minimization(ionized_gro_file, top_file)

        # Get temperature and pressure for equilibration
        temp_input = inquirer.prompt([inquirer.Text('temp', message="Enter the temperature (K)", validate=validate_num)])['temp']
        pressure_input = inquirer.prompt([inquirer.Text('pressure', message="Enter the pressure (atm)", validate=validate_num)])['pressure']

        nvt_gro_file, top_file = nvt_equilibration(gro_file, top_file, temp_input)
        npt_gro_file, top_file = npt_equilibration(gro_file, top_file, temp_input, pressure_input)

        # Run MD simulation
        trajectory_file = md_simulation(nvt_gro_file, top_file)

        # Analysis and Graphing
        result_graphs(trajectory_file)  # Generate graphs

else:
    print('No .pdb file selected. Exiting.')
    exit()
