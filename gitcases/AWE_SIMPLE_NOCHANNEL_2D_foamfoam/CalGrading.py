import math
import CalGrading #CalGrading.py

class AllFunctions:
    def __init__(self, precision=10, relTol=1e-5, rMax=1e5):
        self.precision = precision
        self.relTol = relTol
        self.rMax = rMax

    def sum_cells(self, x0, r, n):
        if n is None or n < 1:
            return 0
        else:
            LENGTH = 0
            x = x0
            for i in range(1, n + 1):
                LENGTH += x
                x *= r
            return LENGTH

    def root_by_bisection(self, f, x1, x2):
        max_steps = 500
        f1 = f(x1)
        f2 = f(x2)
        if f1 * f2 >= 0:
            print(f"{x1} and {x2} do not bracket the root ({f1}, {f2})")
            return None, "Internal error: Wrong start values for root finding"

        steps = 0
        while steps < max_steps:
            steps += 1
            x_mid = 0.5 * (x1 + x2)
            if (x2 - x1) < self.relTol:
                return x_mid
            f_mid = f(x_mid)
            if f1 * f_mid < 0:
                x2 = x_mid
                f2 = f_mid
            else:
                x1 = x_mid
                f1 = f_mid

        return None, "Internal error: root finding did not converge"

    def round_val(self, val):
        return [round(val), float(format(val, f".{self.precision}g"))]

    def calculate_TER_from_NCELLS_CELLTOCELLER(self, NCELLS, CELLTOCELLER):
        if NCELLS > 1:
            return math.pow(CELLTOCELLER, NCELLS - 1)
        else:
            return [CELLTOCELLER, "Number of cells must be >1"]

    def calculate_CELLTOCELLER_from_NCELLS_TER(self, NCELLS, TER):
        if NCELLS > 1:
            return math.pow(TER, 1 / (NCELLS - 1))
        else:
            return [TER, "Number of cells must be >1"]

    def calculate_TER_from_WSTARTCELL_WENDCELL(self, WSTARTCELL, WENDCELL):
        return WENDCELL / WSTARTCELL

    def calculate_WENDCELL_from_WSTARTCELL_TER(self, WSTARTCELL, TER):
        return WSTARTCELL * TER

    def calculate_WSTARTCELL_from_WENDCELL_TER(self, WENDCELL, TER):
        return WENDCELL / TER

    def calculate_LENGTH_from_WSTARTCELL_CELLTOCELLER_NCELLS(self, WSTARTCELL, CELLTOCELLER, NCELLS):
        return self.sum_cells(WSTARTCELL, CELLTOCELLER, NCELLS)

    def calculate_LENGTH_from_WENDCELL_CELLTOCELLER_NCELLS(self, WENDCELL, CELLTOCELLER, NCELLS):
        return self.sum_cells(WENDCELL, 1 / CELLTOCELLER, NCELLS)

    def calculate_NCELLS_from_TER_CELLTOCELLER(self, RATIO, CELLTOCELLER):
        if abs(CELLTOCELLER - 1) > 1e-5:
            return self.round_val(math.log(TER) / math.log(CELLTOCELLER) + 1)
        else:
            return [None, "Cell to cell ratio must not be 1"]

    def calculate_WSTARTCELL_from_NCELLS_CELLTOCELLER_LENGTH(self, NCELLS, CELLTOCELLER, LENGTH):
        if abs(CELLTOCELLER - 1) > 1e-5:
            return LENGTH * (1 - CELLTOCELLER) / (1 - math.pow(CELLTOCELLER, NCELLS))
        else:
            return LENGTH / NCELLS

    def calculate_NCELLS_from_CELLTOCELLER_WSTARTCELL_LENGTH(self, CELLTOCELLER, WSTARTCELL, LENGTH):
        if abs(CELLTOCELLER - 1) > self.relTol:
            return self.round_val(math.log(1 - LENGTH / WSTARTCELL * (1 - CELLTOCELLER)) / math.log(CELLTOCELLER))
        else:
            return self.round_val(LENGTH / WSTARTCELL)

    def calculate_NCELLS_from_CELLTOCELLER_WENDCELL_LENGTH(self, CELLTOCELLER, WENDCELL, LENGTH):
        if abs(CELLTOCELLER - 1) > self.relTol:
            return self.round_val(math.log(1 / (1 + LENGTH / WENDCELL * (1 - CELLTOCELLER) / CELLTOCELLER)) / math.log(CELLTOCELLER))
        else:
            return self.round_val(LENGTH / WENDCELL)

    def calculate_CELLTOCELLER_from_NCELLS_WENDCELL_LENGTH(self, NCELLS, WENDCELL, LENGTH):
        if abs(NCELLS * WENDCELL - LENGTH) / LENGTH < self.relTol:
            return 1
        else:
            if NCELLS * WENDCELL > LENGTH:
                cMax = math.pow(self.rMax, 1 / (NCELLS - 1))
                cMin = math.pow(1 + self.relTol, 1 / (NCELLS - 1))
            else:
                cMax = math.pow(1 - self.relTol, 1 / (NCELLS - 1))
                cMin = math.pow(1 / self.rMax, 1 / (NCELLS - 1))

            return self.root_by_bisection(
                lambda c: (1 / math.pow(c, NCELLS - 1)) * (1 - math.pow(c, NCELLS)) / (1 - c) - LENGTH / WENDCELL,
                cMin,
                cMax
            )

    def calculate_CELLTOCELLER_from_NCELLS_WSTARTCELL_LENGTH(self, NCELLS, WSTARTCELL, LENGTH):
        if abs(NCELLS * WSTARTCELL - LENGTH) / LENGTH < self.relTol:
            return 1
        else:
            if NCELLS * WSTARTCELL < LENGTH:
                cMax = math.pow(self.rMax, 1 / (NCELLS - 1))
                cMin = math.pow(1 + self.relTol, 1 / (NCELLS - 1))
            else:
                cMax = math.pow(1 - self.relTol, 1 / (NCELLS - 1))
                cMin = math.pow(1 / self.rMax, 1 / (NCELLS - 1))

            return self.root_by_bisection(
                lambda c: (1 - math.pow(c, NCELLS)) / (1 - c) - LENGTH / WSTARTCELL,
                cMin, cMax
            )

    def calculate_NCELLS_from_TER_WSTARTCELL_LENGTH(self, TER, WSTARTCELL, LENGTH):
        if abs(TER - 1) < self.relTol:
            return self.round_val(LENGTH / WSTARTCELL)
        else:
            return self.round_val(
                self.root_by_bisection(
                    lambda n: (1 - math.pow(TER, n / (n - 1))) / (1 - math.pow(TER, 1 / (n - 1))) - LENGTH / WSTARTCELL,
                    0, LENGTH / WSTARTCELL
                )
            )

    def calculate(self, function_name, *args):
        if hasattr(self, function_name):
            return getattr(self, function_name)(*args)
        else:
            raise ValueError(f"Function {function_name} not found")
############################################################################################################################################# SPLIT LENGTH CELLS #############################################################################################################################

    def calculate_S_CELLTOCELLER(self, width_start_cell, width_end_cell, S_NCELLS):
        return (width_end_cell / width_start_cell) ** (1 / (S_NCELLS - 1)) if S_NCELLS > 1 else 1

    def calculate_S_TER(self, width_start_cell, width_end_cell):
        return width_end_cell / width_start_cell

    def allocate_cells_by_percentage(self, total_length, sections, cell_distribution, initial_width_start_cell):
        # Ensure the sum of percentages is 100%
        total_percentage = sum(sections)
        if total_percentage != 100:
            raise ValueError(f"Sum of sections percentages should be 100%, but got {total_percentage}%.")

        allocated_cells = []
        remaining_length = total_length

        for i, section in enumerate(sections):
            S_LENGTH = (section / 100) * total_length  # Calculate the length of each section
            S_NCELLS = cell_distribution[i]  # Get the number of cells allocated to the section

            # Calculate the total expansion ratio
            S_TER = self.calculate_S_TER(initial_width_start_cell, initial_width_start_cell * S_LENGTH / initial_width_start_cell)

            # Calculate the cell-to-cell expansion ratio
            S_CELLTOCELLER = self.calculate_S_CELLTOCELLER(initial_width_start_cell, initial_width_start_cell * S_TER, S_NCELLS)

            # Append the calculated ratios for the section
            allocated_cells.append({
                'section': i + 1,
                'S_LENGTH': S_LENGTH,
                'S_NCELLS': S_NCELLS,
                'S_TER': S_TER,
                'S_CELLTOCELLER': S_CELLTOCELLER
            })

            # Update the remaining length for the next sections
            remaining_length -= S_LENGTH

        return allocated_cells

# Example usage
total_length = 1.0  # Total length of the block
sections = [30, 50, 20]  # Sections represent percentages (e.g., 30% of total length)
cell_distribution = [5, 10, 15]  # Allocate 5, 10, and 15 cells to the corresponding sections
initial_width_start_cell = 0.01  # The starting width of the first cell

# Create an instance of AllFunctions
all_funcs = AllFunctions()

# Get the allocated cells for each section with the new condition
allocated_cells = all_funcs.allocate_cells_by_percentage(total_length, sections, cell_distribution, initial_width_start_cell)
#print("Allocated cells per section:", allocated_cells)



############################################################################################################################################# SPLIT LENGTH CELLS #############################################################################################################################




# Example usage
#link to original https://openfoamwiki.net/index.php/Scripts/blockMesh_grading_calculation
if __name__ == "__main__":
    #---- overview of names ----

    #NCELLS -> number of cells
    #CELLTOCELLER -> cell-to-cell ratio // growth between cell
    #TER -> Total expansion ratio  // total growth relative between start and end // width of end cell/width of start cell
    #LENGTH -> distance for growth
    #WSTARTCELL -> width for start cell
    #WENDCELL -> width for end cell
    all_funcs = AllFunctions()
    known_values = {
        "NCELLS": 5,
        "WSTARTCELL": 2e-05,
        "CELLTOCELLER": 2
    }
    known_values["TER"] = all_funcs.calculate('calculate_TER_from_NCELLS_CELLTOCELLER', known_values['NCELLS'], known_values['CELLTOCELLER'])
    known_values["LENGTH"] = all_funcs.calculate('calculate_LENGTH_from_WSTARTCELL_CELLTOCELLER_NCELLS', known_values['WSTARTCELL'], known_values['CELLTOCELLER'], known_values['NCELLS'])
    known_values["WENDCELL"] = all_funcs.calculate('calculate_WENDCELL_from_WSTARTCELL_TER', known_values['WSTARTCELL'], known_values['TER'])

    print("known_values:", known_values)

    initial_values = {
        "NCELLS": 11,
        "WSTARTCELL": 0.01,
        "WENDCELL": 1
    }

    results = all_funcs.calculate('calculate_CELLTOCELLER_from_NCELLS_WSTARTCELL_LENGTH', initial_values['NCELLS'], initial_values['WSTARTCELL'], initial_values['WENDCELL'])

    print("Calculated results:", results)

