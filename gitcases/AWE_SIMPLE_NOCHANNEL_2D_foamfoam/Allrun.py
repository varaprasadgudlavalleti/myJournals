#!/usr/bin/env python3
import subprocess
from pathlib import Path
import importlib.util
from colorama import Fore, Style
from porous import get_table_width  # Import the helper function

def print_separator():
    # Automatically adjust the cyan line length to match the porous table width.
    width = get_table_width()
    print(Fore.CYAN + "─" * width + Style.RESET_ALL)

def print_comment(comment):
    print(Fore.RED + comment + Style.RESET_ALL)

def print_file_content(filepath):
    with open(filepath, "r") as f:
        print(f.read())

def run_shell_command(command, cwd=None):
    subprocess.run(command, shell=True, cwd=cwd, check=True)

def run_python_file(filepath):
    spec = importlib.util.spec_from_file_location("module", filepath)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

def run_python_file_capture(filepath):
    result = subprocess.run(["python3", filepath],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            text=True,
                            cwd=Path(filepath).parent)
    print(result.stdout)
    if result.stderr:
        print(result.stderr)

# Determine the case directory (where Allrun.py resides)
case = Path(__file__).resolve().parent

print_separator()
print_comment("Executing Allclean: Resetting files")
run_shell_command("./Allclean", cwd=case)

print_separator()
print_comment("Executing porous.py: Applying substitutions and printing table")
run_python_file_capture(str(case / "porous.py"))
# Reapply substitutions to ensure they persist
#from porous import run_substitutions
#run_substitutions()

print_separator()
print_comment("Executing species.py: Applying substitutions and printing table")
run_python_file_capture(str(case / "species.py"))

print_separator()
print_comment("Executing mesh.py: Building mesh")
run_python_file(str(case / "mesh.py"))

print_separator()
print_comment("Executing Allrun shell script: Launching simulation")
run_shell_command("./Allrun", cwd=case)
print_separator()
