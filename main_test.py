import os
from functions_test import * # Import the updated functions

workDir = os.getcwd()
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
        gro_file, top_file = pdb2gmx(file, top_file)
        gro_file, top_file = solvate(gro_file, top_file)
        gro_file, top_file = genion(gro_file, top_file, ion_choice)
        gro_file, top_file = energy_minimization(gro_file, top_file)

        # Get temperature and pressure for equilibration
        temp_input = inquirer.prompt(
            [inquirer.Text('temp', message="Enter the temperature (K)", validate=validate_num)])['temp']
        pressure_input = inquirer.prompt(
            [inquirer.Text('pressure', message="Enter the pressure (atm)", validate=validate_num)])['pressure']

        gro_file, top_file = nvt_equilibration(gro_file, top_file, temp_input)
        gro_file, top_file = npt_equilibration(gro_file, top_file, temp_input, pressure_input)

        # Run MD simulation
        trajectory_file = md_simulation(gro_file, top_file)

        # Analysis and Graphing
        result_graphs(trajectory_file)
