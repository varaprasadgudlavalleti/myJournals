#!/usr/bin/env python3
import re
from pathlib import Path
from colorama import Fore, Style

def generate_substitutions():
    """
    Build a dictionary of substitution keys and their corresponding values.
    This includes:
      - Ionic Charges (z_*): Electrical charges for species.
      - Ion Mobilities (u_*): Mobilities in m²/(V·s).
      - Diffusion-Related Constants: Activation energies (Ea), reference diffusion coefficients (D_ref),
        and reference temperatures (T_ref).
      - Molecular Weights (MW): Standard molecular weights for the species.
    """
    subs = {
        # --- Ionic Charges (z_*) ---
        "z_OH": "-1",      # Hydroxide: -1
        "z_waterH2O": "0",      # Water: 0
        "z_H2":  "0",      # Hydrogen: 0
        "z_O2":  "0",      # Oxygen: 0
        "z_K":   "1",      # Potassium: +1

        # --- Ion Mobilities (u_*) in m²/(V·s) ---
        "u_OH": "2.05e-3",
        "u_K":  "7.619e-4",

        # --- Diffusion-Related Constants ---
        # Activation Energies (J/mol)
        "Ea_waterH2O": "16762.54",
        "Ea_OH":   "16514.4",
        "Ea_H2":   "19344.19",
        "Ea_O2":   "20171.89",
        # Reference Diffusion Coefficients (m²/s)
        "D_waterH2O_ref": "2.272e-9",
        "D_OH_ref":  "2.54e-9",
        "D_H2_ref":  "5.8e-9",
        "D_O2_ref":  "1.9e-9",
        # Reference Temperatures (K)
        "T_ref_H2":  "298",
        "T_ref_O2":  "298",
        "T_ref_waterH2O": "353",
        "T_ref_OH":  "353",

        # --- Molecular Weights (kg/mol) ---
        "MW_H2":   "2.01588e-3",
        "MW_O2":   "31.9988e-3",
        "MW_OH":   "17.007e-3",
        "MW_waterH2O":  "18.01528e-3",
        "MW_K":    "39.10e-3",
    }
    return subs

def print_table(title, headers, rows):
    """
    Print a formatted table with a title, given header names and row data.
    """
    # Determine column widths based on header and row content.
    col_widths = [max(len(str(item)) for item in col) for col in zip(headers, *rows)]
    total_width = sum(col_widths) + 3 * (len(headers) - 1) + 4

    # Print title box.
    print("┌" + "─" * (total_width - 2) + "┐")
    print("│" + title.center(total_width - 2) + "│")
    print("├" + "─" * (total_width - 2) + "┤")

    # Print header row.
    header_line = "│ " + " │ ".join(h.ljust(w) for h, w in zip(headers, col_widths)) + " │"
    print(header_line)
    print("├" + "─" * (total_width - 2) + "┤")

    # Print each row of data.
    for row in rows:
        row_line = "│ " + " │ ".join(str(item).ljust(w) for item, w in zip(row, col_widths)) + " │"
        print(row_line)
    print("└" + "─" * (total_width - 2) + "┘\n")

def print_unified_table(subs):
    """
    Build a unified table in which each row corresponds to a species (OH, H2O, H2, O2, K)
    and the columns list the associated properties:
      - Charge (z_*)
      - Mobility (u_*)
      - Activation Energy (Ea)
      - Reference Diffusion Coefficient (D_ref)
      - Reference Temperature (T_ref)
      - Molecular Weight (MW)

    For species where a property is not defined (e.g. mobility for H2O, H2, and O2), a dash ("–") is used.
    """
    species_data = {
        "OH": {
            "Charge": subs["z_OH"],
            "Mobility": subs["u_OH"],
            "Activation Energy": subs["Ea_OH"],
            "D_ref": subs["D_OH_ref"],
            "T_ref": subs["T_ref_OH"],
            "Molecular Weight": subs["MW_OH"]
        },
        "H2O": {
            "Charge": subs["z_waterH2O"],
            "Mobility": "–",
            "Activation Energy": subs["Ea_waterH2O"],
            "D_ref": subs["D_waterH2O_ref"],
            "T_ref": subs["T_ref_waterH2O"],
            "Molecular Weight": subs["MW_waterH2O"]
        },
        "H2": {
            "Charge": subs["z_H2"],
            "Mobility": "–",
            "Activation Energy": subs["Ea_H2"],
            "D_ref": subs["D_H2_ref"],
            "T_ref": subs["T_ref_H2"],
            "Molecular Weight": subs["MW_H2"]
        },
        "O2": {
            "Charge": subs["z_O2"],
            "Mobility": "–",
            "Activation Energy": subs["Ea_O2"],
            "D_ref": subs["D_O2_ref"],
            "T_ref": subs["T_ref_O2"],
            "Molecular Weight": subs["MW_O2"]
        },
        "K": {
            "Charge": subs["z_K"],
            "Mobility": subs["u_K"],
            "Activation Energy": "–",
            "D_ref": "–",
            "T_ref": "–",
            "Molecular Weight": subs["MW_K"]
        }
    }

    headers = [
        "Species",
        "Charge",
        "Mobility (m²/(V·s))",
        "Activation Energy (J/mol)",
        "D_ref (m²/s)",
        "T_ref (K)",
        "Molecular Weight (kg/mol)"
    ]

    rows = []
    for species in ["OH", "H2O", "H2", "O2", "K"]:
        row = [
            species,
            species_data[species]["Charge"],
            species_data[species]["Mobility"],
            species_data[species]["Activation Energy"],
            species_data[species]["D_ref"],
            species_data[species]["T_ref"],
            species_data[species]["Molecular Weight"]
        ]
        rows.append(row)

    print_table("SETTING species PROPERTIES: constant/speciesProperties.orig", headers, rows)

def perform_substitution(orig_file, target_file, subs):
    """
    Reads an original file, performs case-insensitive substitution of placeholders
    that appear in the form "set_<key>" with the corresponding values in subs, and writes
    the result to the target file.
    """
    try:
        with open(orig_file, "r") as f:
            content = f.read()
        # Loop over each substitution key and replace only occurrences that start with "set_"
        for key, replacement in subs.items():
            pattern = re.compile(r"set_" + re.escape(key), re.IGNORECASE)
            content = pattern.sub(replacement, content)
        with open(target_file, "w") as f:
            f.write(content)
        print(f"Substitution done for {target_file}")
    except Exception as e:
        print(f"Error processing {orig_file}: {e}")

def run_substitutions():
    """
    Perform substitutions on a single file named speciesProperties.orig
    and write the result to speciesProperties.
    """
    case = Path(__file__).resolve().parent

    orig_file = case / "constant" / "speciesProperties.orig"
    target_file = case / "constant" / "speciesProperties"

    subs = generate_substitutions()
    perform_substitution(orig_file, target_file, subs)

if __name__ == "__main__":
    # Generate the substitutions dictionary with all constants.
    subs = generate_substitutions()
    # Print the unified table.
    print_unified_table(subs)
    # Run file substitution for speciesProperties.
    run_substitutions()
