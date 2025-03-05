#!/usr/bin/env python3
import os
import subprocess
import math
import re
from pathlib import Path
from colorama import Fore, Style

# ANSI escape codes for colors (no background)
RED   = "\033[31m"
BLUE  = "\033[34m"
RESET = "\033[0m"

def get_table_width(widths=None):
    """
    Compute and return the total width of a table based on column widths.
    """
    if widths is None:
        widths = [10, 14, 12, 10, 24, 14, 30, 30, 14]
    top_border = "┌" + "┬".join("─" * (w + 2) for w in widths) + "┐"
    return len(top_border)

def calculate_properties(ppi, dp, tau, mu, C=5.0):
    """
    Calculate effective porosity, specific surface area, permeability, and Darcy source term.

    Using the modified Kozeny–Carman equation with tortuosity:
        K = (ε³ · dp²) / (16 · C · (1-ε) · τ²)
    and effective porosity computed as:
        ε = (ppi/0.0254)² · π · (dp/2)²
    Specific surface area (ssa) is calculated as:
        ssa = (ppi/0.0254)² · π · dp · τ
    """
    NL = ppi / 0.0254              # Convert ppi to pores/meter.
    NA = NL ** 2                   # Areal pore density (pores per m²).
    Ap = math.pi * (dp / 2)**2      # Pore cross-sectional area.
    epsilon = NA * Ap              # Effective porosity.
    if epsilon >= 1:
        epsilon = 0.9999           # Enforce physical limit.
    ssa = NA * math.pi * dp * tau   # Specific surface area.
    K = (epsilon**3 * dp**2) / (16 * C * (1 - epsilon) * tau**2)
    darcy = mu / K                 # Darcy source term.
    return ppi, dp, tau, ssa, mu, epsilon, K, darcy

def calculate_wire_diameter(ppi, dp):
    """
    Calculate the effective wire (strut) diameter using the difference method.

    In simple terms:
      1. Compute the cell spacing (L) from the pore density:
            L = (25.4 mm) / ppi   (converted to meters)
         For example, for 110 ppi, L ≈ 25.4e-3 / 110 ≈ 231e-6 m.
      2. The pore (open space) is given as dp.
         The remaining space in the cell is (L - dp), split equally on both sides.
         Thus, the effective wire diameter is: d = (L - dp) / 2
    """
    cell_spacing = 0.0254 / ppi  # cell spacing in meters
    wire_diameter = (cell_spacing - dp) / 2
    return wire_diameter, cell_spacing

def calculate_sherwood_values(epsilon, tau, Re):
    """
    Calculate two Sherwood numbers using the following correlations:

    1. Baseline (Diffusion‐Limited):
         Sh_diff_coeff = ε / τ

    2. Porous (with convection):
         Sh = (ε/τ)[1.1 + 0.43 * Re^(0.33)]
    """
    Sh_diff_coeff = epsilon / tau
    Sh = (epsilon / tau) * (1.1 + 0.43 * (Re ** 0.33))
    return Sh_diff_coeff, Sh

def calculate_Re_seepage(K, mu, dp, deltaP=1.0, L_e=1.5e-3, rho=1000):
    """
    Calculate the pore-scale (seepage) Reynolds number.

    Darcy (seepage) velocity: u = (K/μ) * (ΔP / L_e)
    Kinematic viscosity: ν = μ / ρ
    Then: Re_seepage = u * dp / ν

    Parameters:
      K     : Permeability (m²)
      mu    : Dynamic viscosity (Pa·s)
      dp    : Pore diameter (m)
      deltaP: Pressure drop (Pa) [default 1.0]
      L_e   : Electrode thickness (m) [default 1.5e-3]
      rho   : Fluid density (kg/m³) [default 1000]
    """
    u = (K / mu) * (deltaP / L_e)  # Darcy velocity (m/s)
    nu = mu / rho                  # Kinematic viscosity (m²/s)
    Re_seepage = (u * dp) / nu
    return Re_seepage

# Example input parameters for each component.
# We compute Re_seepage dynamically based on the available properties.
components = {
    "PE": {"ppi": 110, "dp": 200e-6, "tau": 1.2, "mu": 2.0e-3, "C": 5.0},
    "DI": {"ppi": 300,  "dp":50e-6, "tau": 1.5, "mu": 2.0e-3, "C": 5.0},
    "NE": {"ppi": 110, "dp": 200e-6, "tau": 1.2, "mu": 2.0e-3, "C": 5.0}
}

# Compute porous properties, wire diameters, and Sherwood numbers.
results = {}
wire_results = {}
sherwood_table = {}
for comp, params in components.items():
    results[comp] = calculate_properties(
        params["ppi"],
        params["dp"],
        params["tau"],
        params["mu"],
        params.get("C", 5.0)
    )
    wire_d, cell_spacing = calculate_wire_diameter(params["ppi"], params["dp"])
    wire_results[comp] = (wire_d, cell_spacing)
    # Calculate Re_seepage using computed permeability (K), viscosity, and pore size.
    Re_seepage = calculate_Re_seepage(results[comp][6], results[comp][4], params["dp"])
    Sh_diff_coeff, Sh = calculate_sherwood_values(results[comp][5], results[comp][2], Re=Re_seepage)
    # Store the computed Re_seepage along with Sh_diff_coeff and Sh.
    sherwood_table[comp] = (Re_seepage, Sh_diff_coeff, Sh)

def get_dynamic_widths_sherwood(sherwood_table):
    """
    Compute dynamic column widths for the sherwood table based on headers and data.
    """
    header_top = ["Comp", "Re_seepage", "Sh_diff_coeff", "Sh"]
    header_bot = ["", "mRe", "ε/τ", "(ε/τ)[1.1+0.43Re^0.33]"]
    data_rows = []
    for comp, vals in sherwood_table.items():
        data_rows.append([comp, f"{vals[0]:.4e}", f"{vals[1]:.4e}", f"{vals[2]:.4e}"])
    num_columns = len(header_top)
    widths = []
    for col in range(num_columns):
        max_width = max(len(header_top[col]), len(header_bot[col]))
        for row in data_rows:
            if len(row[col]) > max_width:
                max_width = len(row[col])
        widths.append(max_width)
    return widths

def get_dynamic_widths_wire(wire_results):
    """
    Compute dynamic column widths for the wire table based on headers and data.
    """
    header_top = ["Comp", "Wire Diameter (d, m)", "Cell Spacing (S, m)"]
    header_bot = ["Comp", "d=(L - dp)/2", "0.0254 / PPI"]
    data_rows = []
    for comp, vals in wire_results.items():
        data_rows.append([comp, f"{vals[0]:.4e}", f"{vals[1]:.4e}"])
    num_columns = len(header_top)
    widths = []
    for col in range(num_columns):
        max_width = max(len(header_top[col]), len(header_bot[col]))
        for row in data_rows:
            if len(row[col]) > max_width:
                max_width = len(row[col])
        widths.append(max_width)
    return widths

def print_table(results):
    """
    Print a table with the computed porous properties for each component.
    """
    widths = [10, 14, 12, 10, 24, 14, 30, 30, 14]

    def border_line(left, mid, right):
        return left + mid.join("─" * (w + 2) for w in widths) + right

    top_border    = border_line("┌", "┬", "┐")
    header_sep    = border_line("├", "┼", "┤")
    bottom_border = border_line("└", "┴", "┘")

    total_width = len(top_border)
    overall_title = ("SETTING POROUS PROPERTIES: eps.orig, tau.orig, porousProperties.orig, fvOptions.orig, setFieldsDict.orig")
    title_border = "─" * (total_width - 2)
    title_box = (
        "┌" + title_border + "┐\n" +
        "│" + overall_title.center(total_width - 2) + "│\n" +
        top_border
    )

    header_top = ["Comp", "Pore Density", "Pore Size", "Tortuosity", "Specific Surface Area", "Viscosity", "Porosity", "Permeability", "Darcy"]
    header_bot = ["", "(PPI)", "(dp, m)", "(τ)", "(ssa, m⁻¹)", "(μ, Pa·s)", "(ε)", "(K, m²)", "(μ/K)"]

    header_line1 = "│ " + " │ ".join(text.center(width) for text, width in zip(header_top, widths)) + " │"
    header_line2 = "│ " + " │ ".join(text.center(width) for text, width in zip(header_bot, widths)) + " │"
    header_rows = header_line1 + "\n" + RED + header_line2 + RESET

    formula_cells = [
        "Comp",
        "input",
        "input",
        "input",
        "ssa=(ppi/0.0254)²·π·dp·τ",
        "input",
        "ε=(ppi/0.0254)²·π*(dp/2)²",
        "K=ε³·dp²/(16C(1-ε)τ²)",
        "μ/K"
    ]
    formula_row = "│ " + " │ ".join(cell.center(width) for cell, width in zip(formula_cells, widths)) + " │"

    data_rows = ""
    for res in results:
        comp, ppi_val, dp_val, tau_val, ssa_val, mu_val, epsilon, K, darcy = res
        row_cells = [
            comp,
            f"{ppi_val:.0f}",
            f"{dp_val:.4e}",
            f"{tau_val:.2f}",
            f"{ssa_val:.4e}",
            f"{mu_val:.4e}",
            f"{epsilon:.4f}",
            f"{K:.4e}",
            f"{darcy:.4e}"
        ]
        data_row = "│ " + " │ ".join(cell.center(width) for cell, width in zip(row_cells, widths)) + " │"
        data_rows += data_row + "\n"

    print(title_box)
    print(RED + header_rows + RESET)
    print(header_sep)
    print(BLUE + formula_row + RESET)
    print(header_sep)
    print(data_rows.rstrip())
    print(bottom_border)

def print_wire_table(wire_results):
    """
    Print a table with the computed wire (strut) diameter and cell spacing for each component.
    This function dynamically adjusts the column widths so that the title row and all other rows are automatically adjusted.
    """
    header_top = ["Comp", "Wire Diameter (d, m)", "Cell Spacing (S, m)"]
    header_bot = ["Comp", "d=(L - dp)/2", "0.0254 / PPI"]
    data = []
    for comp, vals in wire_results.items():
        data.append([comp, f"{vals[0]:.4e}", f"{vals[1]:.4e}"])
    num_columns = len(header_top)
    widths = []
    for col in range(num_columns):
        max_width = max(len(header_top[col]), len(header_bot[col]))
        for row in data:
            if len(row[col]) > max_width:
                max_width = len(row[col])
        widths.append(max_width)

    def border_line(left, mid, right):
        return left + mid.join("─" * (w + 2) for w in widths) + right

    top_border    = border_line("┌", "┬", "┐")
    header_sep    = border_line("├", "┼", "┤")
    bottom_border = border_line("└", "┴", "┘")

    total_width = len(top_border)
    overall_title = "WIRE (STRUT) DIAMETER CALCULATION (Difference Method): transportProperties.orig - d_Ne, d_Pe"
    # Adjust column widths if overall_title is wider than the current table width minus 2
    desired_total_width = len(overall_title) + 2
    if total_width < desired_total_width:
        extra_width = desired_total_width - total_width
        widths[0] += extra_width
        top_border    = border_line("┌", "┬", "┐")
        header_sep    = border_line("├", "┼", "┤")
        bottom_border = border_line("└", "┴", "┘")
        total_width = len(top_border)

    title_border = "─" * (total_width - 2)
    title_box = (
        "┌" + title_border + "┐\n" +
        "│" + overall_title.center(total_width - 2) + "│\n" +
        top_border
    )

    header_line1 = "│ " + " │ ".join(text.center(width) for text, width in zip(header_top, widths)) + " │"
    header_line2 = "│ " + " │ ".join(text.center(width) for text, width in zip(header_bot, widths)) + " │"

    print(title_box)
    print(RED + header_line1 + RESET)
    print(header_sep)
    print(BLUE + header_line2 + RESET)
    print(header_sep)

    data_rows = ""
    for row in data:
        data_row = "│ " + " │ ".join(cell.center(width) for cell, width in zip(row, widths)) + " │"
        data_rows += data_row + "\n"
    print(data_rows.rstrip())
    print(bottom_border)

def print_sherwood_table(sherwood_table):
    """
    Print a separate table with the computed Sherwood numbers for each component.
    The table includes:
      - a header row (in red) with column names,
      - a separate second row (in blue) showing calculation formulas,
      - and the data rows (in white).
    Columns:
      - Comp
      - Re_seepage (computed from available properties)
      - Sh_diff_coeff = ε/τ
      - Sh = (ε/τ)[1.1+0.43Re^0.33]
    This function dynamically adjusts column widths.
    """
    header_top = ["Comp", "Re_seepage", "Sh_diff_coeff", "Sh"]
    header_bot = ["", "mRe", "ε/τ", "(ε/τ)[1.1+0.43Re^0.33]"]
    data = []
    for comp, vals in sherwood_table.items():
        data.append([comp, f"{vals[0]:.4e}", f"{vals[1]:.4e}", f"{vals[2]:.4e}"])
    num_columns = len(header_top)
    widths = []
    for col in range(num_columns):
        max_width = max(len(header_top[col]), len(header_bot[col]))
        for row in data:
            if len(row[col]) > max_width:
                max_width = len(row[col])
        widths.append(max_width)

    def border_line(left, mid, right):
        return left + mid.join("─" * (w + 2) for w in widths) + right

    top_border    = border_line("┌", "┬", "┐")
    header_sep    = border_line("├", "┼", "┤")
    bottom_border = border_line("└", "┴", "┘")

    total_width = len(top_border)
    overall_title = "SHERWOOD NUMBER CALCULATIONS"
    title_border = "─" * (total_width - 2)
    title_box = (
        "┌" + title_border + "┐\n" +
        "│" + overall_title.center(total_width - 2) + "│\n" +
        top_border
    )

    header_line1 = "│ " + " │ ".join(text.center(width) for text, width in zip(header_top, widths)) + " │"
    header_line2 = "│ " + " │ ".join(text.center(width) for text, width in zip(header_bot, widths)) + " │"

    print(title_box)
    print(RED + header_line1 + RESET)
    print(header_sep)
    print(BLUE + header_line2 + RESET)
    print(header_sep)

    data_rows = ""
    for row in data:
        data_row = "│ " + " │ ".join(cell.center(width) for cell, width in zip(row, widths)) + " │"
        data_rows += data_row + "\n"
    print(data_rows.rstrip())
    print(bottom_border)

def generate_substitutions():
    subs = {
        "PE_set_ppi":     results["PE"][0],
        "PE_set_dp":      results["PE"][1],
        "PE_set_tau":     results["PE"][2],
        "PE_set_ssa":     results["PE"][3],
        "PE_set_mu":      results["PE"][4],
        "PE_set_epsilon": results["PE"][5],
        "PE_set_K":       results["PE"][6],
        "PE_set_darcy":   results["PE"][7],
        "PE_set_wire_d":  wire_results["PE"][0],
        "PE_set_cellSpacing": wire_results["PE"][1],
        "PE_set_sherwood": sherwood_table["PE"][1],  # using Sh_diff_coeff for example
        "DI_set_ppi":     results["DI"][0],
        "DI_set_dp":      results["DI"][1],
        "DI_set_tau":     results["DI"][2],
        "DI_set_ssa":     results["DI"][3],
        "DI_set_mu":      results["DI"][4],
        "DI_set_epsilon": results["DI"][5],
        "DI_set_K":       results["DI"][6],
        "DI_set_darcy":   results["DI"][7],
        "DI_set_wire_d":  wire_results["DI"][0],
        "DI_set_cellSpacing": wire_results["DI"][1],
        "DI_set_sherwood": sherwood_table["DI"][1],
        "NE_set_ppi":     results["NE"][0],
        "NE_set_dp":      results["NE"][1],
        "NE_set_tau":     results["NE"][2],
        "NE_set_ssa":     results["NE"][3],
        "NE_set_mu":      results["NE"][4],
        "NE_set_epsilon": results["NE"][5],
        "NE_set_K":       results["NE"][6],
        "NE_set_darcy":   results["NE"][7],
        "NE_set_wire_d":  wire_results["NE"][0],
        "NE_set_cellSpacing": wire_results["NE"][1],
        "NE_set_sherwood": sherwood_table["NE"][1],
    }
    # New substitutions for wire diameter (d) using simplified syntax:
    subs["PE_set_d"] = wire_results["PE"][0]
    subs["DI_set_d"] = wire_results["DI"][0]
    subs["NE_set_d"] = wire_results["NE"][0]
    # New substitutions for Sherwood number (Sh) using simplified syntax:
    subs["PE_set_Sh"] = sherwood_table["PE"][2]
    subs["DI_set_Sh"] = sherwood_table["DI"][2]
    subs["NE_set_Sh"] = sherwood_table["NE"][2]

    # Format values as strings.
    subs["PE_set_ppi"]     = f"{subs['PE_set_ppi']:.0f}"
    subs["PE_set_dp"]      = f"{subs['PE_set_dp']:.4e}"
    subs["PE_set_tau"]     = f"{subs['PE_set_tau']:.2f}"
    subs["PE_set_ssa"]     = f"{subs['PE_set_ssa']:.4e}"
    subs["PE_set_mu"]      = f"{subs['PE_set_mu']:.4e}"
    subs["PE_set_epsilon"] = f"{subs['PE_set_epsilon']:.4f}"
    subs["PE_set_K"]       = f"{subs['PE_set_K']:.4e}"
    subs["PE_set_darcy"]   = f"{subs['PE_set_darcy']:.4e}"
    subs["PE_set_wire_d"]  = f"{subs['PE_set_wire_d']:.4e}"
    subs["PE_set_cellSpacing"] = f"{subs['PE_set_cellSpacing']:.4e}"
    subs["PE_set_sherwood"] = f"{subs['PE_set_sherwood']:.4e}"
    subs["PE_set_d"] = f"{subs['PE_set_d']:.4e}"
    subs["PE_set_Sh"] = f"{subs['PE_set_Sh']:.4e}"

    subs["DI_set_ppi"]     = f"{subs['DI_set_ppi']:.0f}"
    subs["DI_set_dp"]      = f"{subs['DI_set_dp']:.4e}"
    subs["DI_set_tau"]     = f"{subs['DI_set_tau']:.2f}"
    subs["DI_set_ssa"]     = f"{subs['DI_set_ssa']:.4e}"
    subs["DI_set_mu"]      = f"{subs['DI_set_mu']:.4e}"
    subs["DI_set_epsilon"] = f"{subs['DI_set_epsilon']:.4f}"
    subs["DI_set_K"]       = f"{subs['DI_set_K']:.4e}"
    subs["DI_set_darcy"]   = f"{subs['DI_set_darcy']:.4e}"
    subs["DI_set_wire_d"]  = f"{subs['DI_set_wire_d']:.4e}"
    subs["DI_set_cellSpacing"] = f"{subs['DI_set_cellSpacing']:.4e}"
    subs["DI_set_sherwood"] = f"{subs['DI_set_sherwood']:.4e}"
    subs["DI_set_d"] = f"{subs['DI_set_d']:.4e}"
    subs["DI_set_Sh"] = f"{subs['DI_set_Sh']:.4e}"

    subs["NE_set_ppi"]     = f"{subs['NE_set_ppi']:.0f}"
    subs["NE_set_dp"]      = f"{subs['NE_set_dp']:.4e}"
    subs["NE_set_tau"]     = f"{subs['NE_set_tau']:.2f}"
    subs["NE_set_ssa"]     = f"{subs['NE_set_ssa']:.4e}"
    subs["NE_set_mu"]      = f"{subs['NE_set_mu']:.4e}"
    subs["NE_set_epsilon"] = f"{subs['NE_set_epsilon']:.4f}"
    subs["NE_set_K"]       = f"{subs['NE_set_K']:.4e}"
    subs["NE_set_darcy"]   = f"{subs['NE_set_darcy']:.4e}"
    subs["NE_set_wire_d"]  = f"{subs['NE_set_wire_d']:.4e}"
    subs["NE_set_cellSpacing"] = f"{subs['NE_set_cellSpacing']:.4e}"
    subs["NE_set_sherwood"] = f"{subs['NE_set_sherwood']:.4e}"
    subs["NE_set_d"] = f"{subs['NE_set_d']:.4e}"
    subs["NE_set_Sh"] = f"{subs['NE_set_Sh']:.4e}"
    return subs

def perform_substitution(orig_file, target_file, subs):
    """
    Reads an original file, performs case-insensitive substitution of placeholders
    with the values in subs, and writes the result to the target file.
    """
    try:
        with open(orig_file, "r") as f:
            content = f.read()
        for placeholder, replacement in subs.items():
            pattern = re.compile(re.escape(placeholder), re.IGNORECASE)
            content = pattern.sub(replacement, content)
        with open(target_file, "w") as f:
            f.write(content)
        print(f"Substitution done for {target_file}")
    except Exception as e:
        print(f"Error processing {orig_file}: {e}")

# --- Substitution Commands Section ---
case = Path(__file__).resolve().parent

def run_substitutions():
    subs = generate_substitutions()
    eps_orig                = case / "0.orig" / "eps.orig"
    eps_target              = case / "0.orig" / "eps"
    tau_orig                = case / "0.orig" / "tau.orig"
    tau_target              = case / "0.orig" / "tau"
    fvOptions_orig          = case / "constant" / "fvOptions.orig"
    fvOptions_target        = case / "constant" / "fvOptions"
    porousProperties_orig   = case / "constant" / "porousProperties.orig"
    porousProperties_target = case / "constant" / "porousProperties"
    setFieldsDict_orig      = case / "system" / "setFieldsDict.orig"
    setFieldsDict_target    = case / "system" / "setFieldsDict"
    transportProperties_orig = case / "constant" / "transportProperties.orig"
    transportProperties_target = case / "constant" / "transportProperties"

    perform_substitution(eps_orig, eps_target, subs)
    perform_substitution(tau_orig, tau_target, subs)
    perform_substitution(fvOptions_orig, fvOptions_target, subs)
    perform_substitution(porousProperties_orig, porousProperties_target, subs)
    perform_substitution(setFieldsDict_orig, setFieldsDict_target, subs)
    perform_substitution(transportProperties_orig, transportProperties_target, subs)

if __name__ == "__main__":
    table_results = []
    for comp in components.keys():
        vals = calculate_properties(
            components[comp]["ppi"],
            components[comp]["dp"],
            components[comp]["tau"],
            components[comp]["mu"],
            components[comp].get("C", 5.0)
        )
        table_results.append((comp,) + vals)
    print_table(table_results)
    print_wire_table(wire_results)
    print_sherwood_table(sherwood_table)
    run_substitutions()
