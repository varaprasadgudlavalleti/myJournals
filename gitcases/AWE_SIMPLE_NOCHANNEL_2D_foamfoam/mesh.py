import os
import subprocess
import numpy as np
import CalGrading  # CalGrading.py
import math
from pathlib import Path
from colorama import Fore, Style

# Function to run a command and capture the output
def run_command(command, print_output=False, cwd=None):
    result = subprocess.run(command, shell=True, cwd=cwd,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output = result.stdout + result.stderr
    if result.returncode != 0:
        print(f"Error running command: {command}\n{result.stderr}")
    if print_output:
        print(f"Running command: {command}")
        print(output)
    return output

def clean_dir(directory):
    print(f"Cleaning directory: {directory}")
    os.chdir(directory)
    run_command("./Allclean")
    os.chdir("..")

def sed(find, replaceWith, filePath, newFilePath=None):
    if newFilePath is not None:
        return f"sed -e 's|{find}|{replaceWith}|g' {filePath} > {newFilePath}"
    else:
        return f"sed -i 's|{find}|{replaceWith}|g' {filePath}"

# Get the current script's directory dynamically
case = Path(__file__).resolve().parent

# Instantiate the all_funcs object from CalGrading module.
all_funcs = CalGrading.AllFunctions()

# -------------------------- Mesh Calculations --------------------------
# Definitions
H   = 30.0
W   = 30.0
tp  = 1.5
tdi = 0.5
tn  = 1.5
z0  = 0
z1  = 30.0

# Cumulative thickness calculations
t1 = tp + tdi      # positive electrode + diaphragm
t2 = t1 + tn       # positive + diaphragm + negative electrode

# Mesh settings
NX = 80
NY = 140
NZ = 1

# Mesh calculations in X-direction
NPX  = float(round(tp / t2 * NX))
NDIX = float(round((tdi * 2.0) / t2 * NX))
NNX  = float(round(tn / t2 * NX))

# Mesh calculations in Y-direction (here simply the total NY)
NPY  = float(round(NY))
NDIY = float(round(NY))
NNY  = float(round(NY))

# Define dictionaries for X-direction blocks
meshPX = {
    "NCELLS": NPX,
    "WSTARTCELL": 0.05,
    "THICKNESS": tp  # positive electrode thickness
}
meshDIX = {
    "NCELLS": NDIX,
    "WSTARTCELL": 0.05,
    "THICKNESS": tdi  # diaphragm thickness
}
meshNX = {
    "NCELLS": NNX,
    "WSTARTCELL": 0.05,
    "THICKNESS": tn  # negative electrode thickness
}

# Define dictionaries for Y-direction blocks
meshPY = {
    "NCELLS": NPY,
    "WSTARTCELL": 0.2,
    "LENGTH": H  # positive electrode length
}
meshDIY = {
    "NCELLS": NDIY,
    "WSTARTCELL": 0.2,
    "LENGTH": H  # diaphragm length
}
meshNY = {
    "NCELLS": NNY,
    "WSTARTCELL": 0.2,
    "LENGTH": H  # negative electrode length
}

# Calculate additional mesh parameters using CalGrading functions
# For X-direction blocks
meshPX["CELLTOCELLER"] = all_funcs.calculate('calculate_CELLTOCELLER_from_NCELLS_WSTARTCELL_LENGTH',
                                               meshPX['NCELLS'], meshPX['WSTARTCELL'], meshPX['THICKNESS'])
meshPX['TER'] = all_funcs.calculate('calculate_TER_from_NCELLS_CELLTOCELLER',
                                      meshPX['NCELLS'], meshPX['CELLTOCELLER'])
meshDIX["CELLTOCELLER"] = all_funcs.calculate('calculate_CELLTOCELLER_from_NCELLS_WSTARTCELL_LENGTH',
                                                meshDIX['NCELLS'], meshDIX['WSTARTCELL'], meshDIX['THICKNESS'])
meshDIX['TER'] = all_funcs.calculate('calculate_TER_from_NCELLS_CELLTOCELLER',
                                       meshDIX['NCELLS'], meshDIX['CELLTOCELLER'])
meshNX["CELLTOCELLER"] = all_funcs.calculate('calculate_CELLTOCELLER_from_NCELLS_WSTARTCELL_LENGTH',
                                               meshNX['NCELLS'], meshNX['WSTARTCELL'], meshNX['THICKNESS'])
meshNX['TER'] = all_funcs.calculate('calculate_TER_from_NCELLS_CELLTOCELLER',
                                      meshNX['NCELLS'], meshNX['CELLTOCELLER'])

# For Y-direction blocks
meshPY["CELLTOCELLER"] = all_funcs.calculate('calculate_CELLTOCELLER_from_NCELLS_WSTARTCELL_LENGTH',
                                               meshPY['NCELLS'], meshPY['WSTARTCELL'], meshPY['LENGTH'])
meshPY['TER'] = all_funcs.calculate('calculate_TER_from_NCELLS_CELLTOCELLER',
                                      meshPY['NCELLS'], meshPY['CELLTOCELLER'])
meshDIY["CELLTOCELLER"] = all_funcs.calculate('calculate_CELLTOCELLER_from_NCELLS_WSTARTCELL_LENGTH',
                                                meshDIY['NCELLS'], meshDIY['WSTARTCELL'], meshDIY['LENGTH'])
meshDIY['TER'] = all_funcs.calculate('calculate_TER_from_NCELLS_CELLTOCELLER',
                                       meshDIY['NCELLS'], meshDIY['CELLTOCELLER'])
meshNY["CELLTOCELLER"] = all_funcs.calculate('calculate_CELLTOCELLER_from_NCELLS_WSTARTCELL_LENGTH',
                                               meshNY['NCELLS'], meshNY['WSTARTCELL'], meshNY['LENGTH'])
meshNY['TER'] = all_funcs.calculate('calculate_TER_from_NCELLS_CELLTOCELLER',
                                      meshNY['NCELLS'], meshNY['CELLTOCELLER'])

# -------------------------- Table Printing Functions --------------------------
def print_table(title, rows, headers):
    # Dynamically calculate column widths based on headers and data
    num_cols = len(headers)
    col_widths = []
    for col in range(num_cols):
        items = [str(headers[col])] + [str(row[col]) for row in rows]
        col_widths.append(max(len(item) for item in items))

    total_width = sum(col_widths) + 3 * (num_cols - 1) + 4  # borders and spaces

    # Build title box (without an extra table top border)
    title_border = "─" * (total_width - 2)
    title_box = (
        "┌" + title_border + "┐\n" +
        "│" + title.center(total_width - 2) + "│\n" +
        "└" + title_border + "┘"
    )
    print(title_box)

    # Immediately print header row (attached to the title box)
    header_row = "│ " + " │ ".join(h.center(w) for h, w in zip(headers, col_widths)) + " │"
    header_sep = "├" + "┬".join("─" * (w + 2) for w in col_widths) + "┤"
    print(header_row)
    print(header_sep)

    # Print each data row
    for row in rows:
        row_str = "│ " + " │ ".join(str(cell).center(w) for cell, w in zip(row, col_widths)) + " │"
        print(row_str)

    table_bottom = "└" + "┴".join("─" * (w + 2) for w in col_widths) + "┘"
    print(table_bottom)
    print()  # extra newline after table for clarity

# -------------------------- Prepare Data for Tables --------------------------
# For X-direction blocks:
x_rows = [
    ["POSITIVE", meshPX["NCELLS"], meshPX["WSTARTCELL"], meshPX["THICKNESS"], meshPX["CELLTOCELLER"], meshPX["TER"]],
    ["DIAPHRAGM", meshDIX["NCELLS"], meshDIX["WSTARTCELL"], meshDIX["THICKNESS"], meshDIX["CELLTOCELLER"], meshDIX["TER"]],
    ["NEGATIVE", meshNX["NCELLS"], meshNX["WSTARTCELL"], meshNX["THICKNESS"], meshNX["CELLTOCELLER"], meshNX["TER"]],
]

# For Y-direction blocks:
y_rows = [
    ["POSITIVE", meshPY["NCELLS"], meshPY["WSTARTCELL"], meshPY["LENGTH"], meshPY["CELLTOCELLER"], meshPY["TER"]],
    ["DIAPHRAGM", meshDIY["NCELLS"], meshDIY["WSTARTCELL"], meshDIY["LENGTH"], meshDIY["CELLTOCELLER"], meshDIY["TER"]],
    ["NEGATIVE", meshNY["NCELLS"], meshNY["WSTARTCELL"], meshNY["LENGTH"], meshNY["CELLTOCELLER"], meshNY["TER"]],
]

headers = ["Block", "NCELLS", "WSTARTCELL", "Length", "CELLTOCELLER", "TER"]

# -------------------------- Print Mesh Tables --------------------------
print_table("MESH CELLS in X-DIRECTION", x_rows, headers)
print_table("MESH CELLS in Y-DIRECTION", y_rows, headers)

# -------------------------- File Substitutions and Execution --------------------------
# Construct blockMeshDict paths
blockMeshDictOrig = case / "system" / "blockMeshDict.orig"
blockMeshDict = case / "system" / "blockMeshDict"

# Apply substitutions for blockMeshDict using sed commands
run_command(sed('N_CELLS_POSITIVEELECTRODE_X', meshPX['NCELLS'], blockMeshDictOrig, blockMeshDict))
run_command(sed('EXPANSION_RATIO_POSITIVEELECTRODE_X', meshPX["TER"], blockMeshDict))
run_command(sed('N_CELLS_DIAPHRAGM_X', meshDIX['NCELLS'], blockMeshDict))
run_command(sed('EXPANSION_RATIO_DIAPHRAGM_X', meshDIX["TER"], blockMeshDict))
run_command(sed('N_CELLS_NEGATIVEELECTRODE_X', meshNX['NCELLS'], blockMeshDict))
run_command(sed('EXPANSION_RATIO_NEGATIVEELECTRODE_X', meshNX["TER"], blockMeshDict))

run_command(sed('N_CELLS_POSITIVEELECTRODE_Y', meshPY['NCELLS'], blockMeshDict))
run_command(sed('EXPANSION_RATIO_POSITIVEELECTRODE_Y', meshPY["TER"], blockMeshDict))
run_command(sed('N_CELLS_DIAPHRAGM_Y', meshDIY['NCELLS'], blockMeshDict))
run_command(sed('EXPANSION_RATIO_DIAPHRAGM_Y', meshDIY["TER"], blockMeshDict))
run_command(sed('N_CELLS_NEGATIVEELECTRODE_Y', meshNY['NCELLS'], blockMeshDict))
run_command(sed('EXPANSION_RATIO_NEGATIVEELECTRODE_Y', meshNY["TER"], blockMeshDict))

# Run the Allclean and Allrun commands
#run_command("./Allclean", cwd=case)
#run_command("./Allrun", cwd=case)
