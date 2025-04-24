import inquirer
import os
import re
import subprocess  # For running GROMACS commands

answers = None


def ex_con(temp_str):
    """Exit program if 'Exit' is in the input."""
    if 'Exit' in temp_str:
        exit()


def validate_num(temp_str):
    """Validates if the input is a positive number (optionally with K/k)."""

    if not isinstance(temp_str, str):
        return False
    match = re.fullmatch(r'(\d+)\s*[Kk]?', temp_str.strip())
    if not match:
        return False
    try:
        temp_num = int(match.group(1))
        return temp_num > 0
    except ValueError:
        return False


def select_pdb(workdir):
    """Selects a .pdb file from the Proteins folder."""

    global answers
    protein_dir = os.path.join(workdir, "Proteins")
    if not os.path.isdir(protein_dir):
        print("Proteins folder is missing or not a directory. Please create it.")
        exit()

    try:
        read_sel = [f for f in os.listdir(protein_dir) if f.endswith(".pdb")]
    except FileNotFoundError:
        print(f"Error: Directory not found: {protein_dir}")
        exit()
    except PermissionError:
        print(f"Error: Permission denied accessing: {protein_dir}")
        exit()
    if not read_sel:
        print("No .pdb files found in the Proteins folder.")
        exit()

    questions = [
        inquirer.List(
            'file',
            message="Select the .pdb file from the Proteins folder",
            choices=read_sel,
            validate=lambda _, x: ex_con(x),
        ),
    ]
    answers = inquirer.prompt(questions)
    if answers:
        file_loc = os.path.join(protein_dir, answers['file'])
        return file_loc
    else:
        exit()


def select_topology(workdir):
    """Selects or customizes topology settings."""

    global answers
    questions = [
        inquirer.List(
            'top_choice',
            message='Would you like the program to automatically configure the topology based off the project, or would you like to customize?',
            choices=['Auto', 'Customize', 'Exit'],
            default='Auto',
            validate=lambda _, x: ex_con(x),
        ),
    ]
    answers = inquirer.prompt(questions)

    if answers:
        top_choice = answers['top_choice']
        if top_choice == 'Auto':
            questions = [
                inquirer.List(
                    'project_type',
                    message="Select the MD simulation type",
                    choices=[
                        'Thermal simulations',
                        'Ligand binding simulations',
                        'Nucleic acids binding simulation',
                        'Exit',
                    ],
                    validate=lambda _, x: ex_con(x),
                ),
            ]
            answers = inquirer.prompt(questions)
            if answers:
                project_type = answers['project_type']
                topology_settings = {  # Define topologies here
                    'Thermal simulations': 'OPLS-AA/L',
                    'Ligand binding simulations': 'CHARMM27',
                    'Nucleic acids binding simulation': 'AMBER96',
                }
                if project_type in topology_settings:
                    return ('Auto', topology_settings[project_type])  # Return 'Auto' and the topology
                else:
                    exit()
        elif top_choice == 'Customize':
            top_dir = os.path.join(workdir, "Topology")
            try:
                top_files = [
                    f
                    for f in os.listdir(top_dir)
                    if f.endswith(
                        ".ff"
                    )
                ]  # or other topology file extensions
            except FileNotFoundError:
                print(
                    "Topology folder not found. Please create a 'Topology' folder in the working directory."
                )
                exit()
            except NotADirectoryError:
                print(
                    "Error: 'Topology' is not a directory. Please ensure it's a folder."
                )
                exit()
            if not top_files:
                print("No force field files found in the Topology folder.")
                exit()
            questions = [
                inquirer.List(
                    'custom_topology',
                    message="Choose the force field:",
                    choices=top_files + ['Exit'],
                    validate=lambda _, x: ex_con(x),
                ),
            ]
            answers = inquirer.prompt(questions)
            if answers and answers['custom_topology'] != 'Exit':
                return ('Custom', answers['custom_topology'])  # Return 'Custom' and the file
            else:
                exit()
        else:
            exit()
    else:
        exit()


def select_ions():
    """Selects ions to add to the system."""

    questions = [
        inquirer.Checkbox(
            'ions',
            message="Select which ions you would like to add to the system. The amount of each ion needed for a neutral system will be automatically configured.",
            choices=['Na+', 'Cl-', 'K+', 'Ca+2'],
        )
    ]
    answers = inquirer.prompt(questions)
    if answers:
        return answers['ions']
    else:
        return []


def run_gmx_command(command, input_files=None, output_files=None, input_data=None):
    """
    Executes a GROMACS command.

    Args:
        command (list): The GROMACS command as a list of strings.
        input_files (dict, optional): Dictionary of input files (e.g., {'-f': 'input.mdp'}).
        output_files (dict, optional): Dictionary of output files (e.g., {'-o': 'output.gro'}).
        input_data (str, optional): Data to pipe into the command's stdin.

    Returns:
        tuple: (stdout, stderr, returncode)
    """
    full_command = ['gmx', *command]

    if input_files:
        full_command.extend(
            [key, value] for key, value in input_files.items()
        )  # Add input file flags
    if output_files:
        full_command.extend(
            [key, value] for key, value in output_files.items()
        )  # Add output file flags

    process = subprocess.Popen(
        full_command,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    stdout, stderr = process.communicate(input=input_data)

    return stdout, stderr, process.returncode


def pdb2gmx(pdb_file, top_file):
    """Runs pdb2gmx to generate .gro and .top files."""

    output_files = {'-o': 'protein.gro', '-p': 'topol.top'}
    input_files = {'-f': pdb_file, '-i': 'pdb2gmx.itp'}  # Include itp file
    stdout, stderr, returncode = run_gmx_command(
        ['pdb2gmx', '-ff', top_file],
        input_files=input_files,
        output_files=output_files,
        input_data='4',  # Choose force field (e.g., OPLS-AA/L) - You might need to make this configurable
    )
    if returncode != 0:
        print(f"Error in pdb2gmx: {stderr}")
        exit()
    return 'protein.gro', 'topol.top'



def solvate(gro_file, top_file):
    """Runs solvate to add water to the system."""
    # Create a box
    run_gmx_command(
        ['editconf', '-f', gro_file, '-o', 'box.gro', '-c', '-d', '1.0', '-bt', 'cubic'],
    )  # Adjust box dimensions as needed

    # Solvate
    output_files = {'-o': 'solvated.gro', '-p': 'topol.top'}
    input_files = {'-cp': 'box.gro', '-cs': 'spc216.gro', '-p': top_file}
    stdout, stderr, returncode = run_gmx_command(
        ['solvate'], input_files=input_files, output_files=output_files
    )
    if returncode != 0:
        print(f"Error in solvate: {stderr}")
        exit()
    return 'solvated.gro', 'topol.top'


def genion(gro_file, top_file, ion_choice):
    """Runs genion to add ions to neutralize the system."""
    output_files = {'-o': 'ionized.gro', '-p': 'topol.top'}
    input_files = {'-p': top_file, '-s': 'solvated.gro'}
    command = ['genion', '-conc', '0.1']  # Example: 0.1 M ion concentration
    if 'Na+' in ion_choice:
        command.extend(['-pname', 'NA', '-nname', 'CL'])
    elif 'Cl-' in ion_choice:
        command.extend(['-pname', 'K', '-nname', 'CL'])
    elif 'K+' in ion_choice:
        command.extend(['-pname', 'K', '-nname', 'CL'])
    elif 'Ca+2' in ion_choice:
        command.extend(['-pname', 'CA', '-nname', 'CL'])
    else:
        print('No ions selected')
        return gro_file, top_file
    stdout, stderr, returncode = run_gmx_command(
        command,
        input_files=input_files,
        output_files=output_files,
        input_data='0',  # Neutralize the system
    )
    if returncode != 0:
        print(f"Error in genion: {stderr}")
        exit()
    return 'ionized.gro', 'topol.top'


def energy_minimization(gro_file, top_file):
    """Runs energy minimization."""

    # Generate .mdp file (you might want to provide a default or let the user customize)
    with open('em.mdp', 'w') as f:
        f.write(
            """;
; Run control
title                    = Energy Minimization
integrator               = steep
emtol                    = 10.0
emstep                   = 0.01
nsteps                   = 50000
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 500
nstenergy                = 500
; Neighborsearching
nstlist                  = 10
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.0
; Electrostatics
coulombtype              = PME
pme_order                = 4
fourierspacing           = 0.16
; Temperature coupling
tcoupl                   = No
; Pressure coupling
pcoupl                   = No
; Periodic boundary conditions
pbc                      = xyz
; Dispersion correction
dispcorr                 = EnerPres
; Velocity generation
gen_vel                  = no
"""
        )  # Basic EM parameters

    # Run grompp to create .tpr file
    run_gmx_command(
        ['grompp', '-f', 'em.mdp', '-c', gro_file, '-p', top_file, '-o', 'em.tpr']
    )

    # Run energy minimization
    output_files = {'-o': 'em.gro'}
    input_files = {'-s': 'em.tpr'}
    stdout, stderr, returncode = run_gmx_command(
        ['mdrun', '-v'], input_files=input_files, output_files=output_files
    )  # -v for verbose output
    if returncode != 0:
        print(f"Error in energy minimization: {stderr}")
        exit()
    return 'em.gro', top_file



def nvt_equilibration(gro_file, top_file, temp):
    """Runs NVT equilibration."""
    # Create .mdp file for NVT
    with open('nvt.mdp', 'w') as f:
        f.write(
            f""";
; Run control
title                    = NVT equilibration
integrator               = md
dt                       = 0.002
nsteps                   = 50000 ; 100 ps
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 500
nstenergy                = 500
; Neighborsearching
nstlist                  = 10
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.0
; Electrostatics
coulombtype              = PME
pme_order                = 4
fourierspacing           = 0.16
; Temperature coupling
tcoupl                   = V-rescale
tc-grps                  = System
tau_t                    = 0.1
ref_t                    = {temp}
; Pressure coupling
pcoupl                   = No
; Periodic boundary conditions
pbc                      = xyz
; Dispersion correction
dispcorr                 = EnerPres
; Velocity generation
gen_vel                  = yes
gen_temp                 = {temp}
"""
        )

    # Run grompp
    run_gmx_command(
        ['grompp', '-f', 'nvt.mdp', '-c', gro_file, '-p', top_file, '-o', 'nvt.tpr']
    )

    # Run mdrun
    output_files = {'-o': 'nvt.gro'}
    input_files = {'-s': 'nvt.tpr'}
    stdout, stderr, returncode = run_gmx_command(
        ['mdrun', '-v'], input_files=input_files, output_files=output_files
    )
    if returncode != 0:
        print(f"Error in NVT equilibration: {stderr}")
        exit()
    return 'nvt.gro', top_file



def npt_equilibration(gro_file, top_file, temp, pressure):
    """Runs NPT equilibration."""

    # Create .mdp file for NPT
    with open('npt.mdp', 'w') as f:
        f.write(
            f""";
; Run control
title                    = NPT equilibration
integrator               = md
dt                       = 0.002
nsteps                   = 50000 ; 100 ps
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 500
nstenergy                = 500
; Neighborsearching
nstlist                  = 10
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.0
; Electrostatics
coulombtype              = PME
pme_order                = 4
fourierspacing           = 0.16
; Temperature coupling
tcoupl                   = V-rescale
tc-grps                  = System
tau_t                    = 0.1
ref_t                    = {temp}
; Pressure coupling
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau_p                    = 2.0
ref_p                    = {pressure}
compressibility          = 4.5e-5
; Periodic boundary conditions
pbc                      = xyz
; Dispersion correction
dispcorr                 = EnerPres
; Velocity generation
gen_vel                  = no
gen_temp                 = {temp}
"""
        )

    # Run grompp
    run_gmx_command(
        ['grompp', '-f', 'npt.mdp', '-c', gro_file, '-p', top_file, '-o', 'npt.tpr']
    )

    # Run mdrun
    output_files = {'-o': 'npt.gro'}
    input_files = {'-s': 'npt.tpr'}
    stdout, stderr, returncode = run_gmx_command(
        ['mdrun', '-v'], input_files=input_files, output_files=output_files
    )
    if returncode != 0:
        print(f"Error in NPT equilibration: {stderr}")
        exit()
    return 'npt.gro', top_file



def md_simulation(gro_file, top_file):
    """Runs the molecular dynamics simulation."""

    # Create .mdp file for MD
    with open('md.mdp', 'w') as f:
        f.write(
            """; Run control
title                    = Molecular dynamics
integrator               = md
dt                       = 0.002
nsteps                   = 500000 ; 1 ns (adjust as needed)
; Output control
nstxout                  = 500
nstvout                  = 500
nstfout                  = 500
nstenergy                = 500
nstlog                   = 500
; Neighborsearching
nstlist                  = 20
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.0
; Electrostatics
coulombtype              = PME
pme_order                = 4
fourierspacing           = 0.16
; Temperature coupling
tcoupl                   = V-rescale
tc-grps                  = System
tau_t                    = 0.1
ref_t                    = 300
; Pressure coupling
pcoupl                   = Parrinello-Rahman
pcoupltype               = isotropic
tau_p                    = 2.0
ref_p                    = 1.0
compressibility          = 4.5e-5
; Periodic boundary conditions
pbc                      = xyz
; Dispersion correction
dispcorr                 = EnerPres
; Velocity generation
gen_vel                  = no
"""
        )

    # Run grompp
    run_gmx_command(
        ['grompp', '-f', 'md.mdp', '-c', gro_file, '-p', top_file, '-o', 'md.tpr']
    )

    # Run mdrun
    output_files = {'-o': 'md.gro', '-x': 'md.xtc', '-g': 'md.log'}  # Important outputs
    input_files = {'-s': 'md.tpr'}
    stdout, stderr, returncode = run_gmx_command(
        ['mdrun', '-v'], input_files=input_files, output_files=output_files
    )
    if returncode != 0:
        print(f"Error in MD simulation: {stderr}")
        exit()
    return 'md.xtc'  # Return the trajectory file



def graph_choice():
    """Choose which graphs, if any, you would like to generate for data analysis."""
    questions = [
        inquirer.Checkbox(
            'graphs',
            message="Select which graphs you would like to generate (GROMACS tools will be used):",
            choices=[
                'RMSD',
                'RMSF',
                'Radius of Gyration',
                'Hydrogen Bonds',
                'Exit'
            ],
        )
    ]
    answers = inquirer.prompt(questions)
    if answers:
        if 'Exit' in answers['graphs']:
            exit()
        return answers['graphs']
    else:
        return []


def result_graphs(trajectory_file):
    """Generate graphs using GROMACS analysis tools."""
    # Example: RMSD
    graphs = graph_choice()
    if not graphs:
        return
    if 'RMSD' in graphs:
        run_gmx_command(
            ['gmx', 'rms', '-f', trajectory_file, '-s', 'md.tpr', '-o', 'rmsd.xvg'],
            input_data='0',
        )  # 0 for protein
        print("RMSD graph generated (rmsd.xvg)")
    if 'RMSF' in graphs:
        run_gmx_command(
            ['gmx', 'rmsf', '-f', trajectory_file, '-s', 'md.tpr', '-o', 'rmsf.xvg'],
            input_data='0',
        )
        print("RMSF graph generated (rmsf.xvg)")
    if 'Radius of Gyration' in graphs:
        run_gmx_command(
            ['gmx', 'gyrate', '-f', trajectory_file, '-s', 'md.tpr', '-o', 'gyrate.xvg'],

        )
        print("Radius of Gyration graph generated (gyrate.xvg)")

    if 'Hydrogen Bonds' in graphs:
        run_gmx_command(
            ['gmx', 'hbond', '-f', trajectory_file, '-s', 'md.tpr', '-o', 'hbond.xvg'],

        )
        print("Hydrogen Bond graph generated (hbond.xvg)")
    # Add other analysis tools as needed (e.g., gmx rmsf, gmx gyrate, gmx anaeig)
    #  Remember to handle user input for selecting groups for analysis.
    return