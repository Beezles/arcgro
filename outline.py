#pdb file select
import gmxapi as gmx
import inquirer
import os

#Select pdb from menu
def select_pdb(workdir):
    protein_dir = os.path.join(workdir, "Proteins")
    try:
        access_check = os.access(protein_dir, os.F_OK)  #check if accessible
    except:
        print("Proteins folder is missing or not accessible")
    print("your path is: \n" + protein_dir + "\n")
    read_sel=os.listdir(protein_dir) #read files in proteins folder
    questions = [
        inquirer.List('file',
                        message="Select the .pdb file from the Proteins folder",
                        choices=read_sel,
                        ),
    ]
    answers = inquirer.prompt(questions)
    if answers:
        print(f"You selected: {answers['file']}")
        file_loc = os.path.join(protein_dir, answers['file'])
        read_check = os.access(file_loc, os.R_OK) #check if file is readable
        if read_check is False:
            print(f"ERROR: The file {answers['file']} is NOT readable")
    return file_loc if answers else None

simulation_input = gmx.read_tpr(input)
md = gmx.mdrun(simulation_input)

def pdb(file.pdb):
    input: test_clean.pdb
    return test_pro.gro, topol.top
#topology file preparation and auto selection, give option of custom topology options
def top():
    #print("do you want auto topology selection, or would you like to customize your topology selection?")
    top_choice = inquirer.prompt([
        {
            type: 'checkbox'
            message: 'do you want auto topology selection, or would you like to customize your topology selection?'
            name: 'top_sel'
            choices: [
                new.inquirer
            ]
        }
    ])
#solvation system parameters
def solvate(pro.gro):
    input: test_pro.gro 
    output: test_solvate.gro & topol.top (updated)
#ionization
def ionz(solvate.gro):
    input: test_solvate.gro
    output: test_ions.gro (ionization of system) & topol.top (updated)
def energy_minimization(ions.gro):
    input: test_ions.gro
    output: em.gro (energy minimized system) & topol.top (updated)
def NPT_NVT(em.gro):
#consider parallelization
    input: em.gro
    output: nvt.gro, npt.gro, topol.top (updated)
def MD_sim(npt.gro):
    input: npt.gro
    output: md_0_1.xtc (trajectory, all data which is analyzed to get results) & topol.top (updated)