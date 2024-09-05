import re
import os
import yaml
import pandas as pd
from tabulate import tabulate


# find all shell scripts recursively starting from the given directory
def find_matching_files(directory, extensions):
    shell_scripts = []
    for root, _, files in os.walk(directory):
        for file in files:
            if any(file.endswith(ext) for ext in extensions):
                shell_scripts.append(os.path.join(root, file))
    if len(shell_scripts) == 0:
        print("No shell scripts found using parent directory: " + directory)
    return shell_scripts


# find all shell function names from a shell script
def extract_function_names(script):
    function_pattern = re.compile(r"^\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*\(\)\s*{")
    function_names = []

    for line in script.split("\n"):
        match = function_pattern.match(line)
        if match:
            function_names.append(match.group(1))

    return function_names


# copy function contents from a shell script into string
def extract_shell_function(script, function_name):
    start_pattern = f"^\\s*{function_name}\\s*\\(\\)\\s*{{"
    end_pattern = "^}"
    function_started = False
    function_lines = []

    for line in script.split("\n"):
        if re.match(start_pattern, line):
            function_started = True
        if function_started:
            function_lines.append(line)
        if re.match(end_pattern, line) and function_started:
            break

    return "\n".join(function_lines)


# extract all shell variables from a shell script
def extract_shell_variables(script):
    pattern = r"\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?"
    variables = re.findall(pattern, script)
    defined_vars = list_defined_variables(script_str)
    return list(set(variables + defined_vars))


# recursively extract all defined variables from a yaml file
def extract_yaml_variables(yaml_file):
    with open(yaml_file, "r") as file:
        yaml_data = yaml.safe_load(file)

    variables = set()

    def extract_variables(data):
        if isinstance(data, dict):
            for key, value in data.items():
                variables.add(key)
                extract_variables(value)
        elif isinstance(data, list):
            for item in data:
                extract_variables(item)

    extract_variables(yaml_data)

    return list(variables)


# extract all defined variables from a shell script
def list_defined_variables(script):
    defined_vars = set()
    lines = script.split("\n")

    for line in lines:
        # Remove comments
        line = re.sub(r"#.*", "", line)

        # Extract variable assignments
        assignment_match = re.match(r"\s*([a-zA-Z_][a-zA-Z0-9_]*)\s*=", line)
        if assignment_match:
            var_name = assignment_match.group(1)
            defined_vars.add(var_name)

    return list(defined_vars)


# append the script name to the SourceScripts column for each variable
# define the meta dict if this is the first time encountering a variable
def append_source_script(variable_dict, script_suffix, script_str):
    if script_suffix.endswith(".sh"):
        all_variables = extract_shell_variables(script_str)
    else:
        all_variables = extract_yaml_variables(script_suffix)
    for variable in all_variables:
        if variable not in variable_dict.keys():
            meta_dict = {}
            meta_dict["SourceScripts"] = []
            meta_dict["DefinedIn"] = []
            meta_dict["DefinedInFunction"] = []
            meta_dict["UsedInFunction"] = []
            variable_dict[variable] = meta_dict
        else:
            meta_dict = variable_dict[variable]
        meta_dict["SourceScripts"].append(script_suffix)
    return variable_dict


# append the script name to the DefinedIn column for each defined variable
# all variables defined in yaml files are considered defined variables
def append_defined_script(variable_dict, script_suffix, script_str):
    if script_suffix.endswith(".sh"):
        defined_vars = list_defined_variables(script_str)
    else:
        defined_vars = extract_yaml_variables(script_suffix)
    for defined_var in defined_vars:
        meta_dict = variable_dict[defined_var]
        meta_dict["DefinedIn"].append(script_suffix)
    return variable_dict


# append the function name to the DefinedInFunction column for each defined variable
def append_function_info(variable_dict, script_str):
    defined_functions = extract_function_names(script_str)
    for func_name in defined_functions:
        shell_func = extract_shell_function(script_str, func_name)
        func_vars = extract_shell_variables(shell_func)
        func_defined_vars = list_defined_variables(shell_func)
        for func_var in func_vars:
            meta_dict = variable_dict[func_var]
            meta_dict["UsedInFunction"].append(func_name)
            if func_var in func_defined_vars:
                meta_dict["DefinedInFunction"].append(func_name)
    return variable_dict


if __name__ == "__main__":
    """
    Description:
        This script generates a Markdown table of all shell variables used in the IMI codebase.
        It must be run from the top level directory of the IMI codebase. It is preconfigured to
        ignore certain files and directories based on the specified blacklist. The outputted
        markdown file is saved to src/components/shell_variable_manifest.md
    Usage:
        python src/utilities/generate_var_manifest.py
    """
    # blacklist of files and directories not considered for the manifest
    blacklist = [
        "geoschem_config.yml",
        "resources/containers/",
        ".github/",
        "conda_env.yml",
        ".readthedocs.yaml",
    ]
    cwd = os.getcwd()
    output_file = f"{cwd}/src/components/shell_variable_manifest.md"

    # check that script is being run from top level directory of IMI codebase
    if cwd.split("/")[-1] != "integrated_methane_inversion":
        error_str = (
            "Please run this script from the top level directory "
            + "of the integrated_methane_inversion directory"
        )
        raise Exception(error_str)

    # find all shell scripts
    shell_scripts_found = find_matching_files(cwd, [".sh", ".yml", "yaml"])

    # remove any blacklisted files
    shell_scripts_found = [
        item
        for item in shell_scripts_found
        if not any(black in item for black in blacklist)
    ]

    # scrape all shell scripts for variables and define the columns of the table
    variable_dict = {}
    for script in shell_scripts_found:
        script_str = open(script).read()
        script_suffix = script.split("integrated_methane_inversion/")[-1]

        # Add SourceScripts column
        variable_dict = append_source_script(variable_dict, script_suffix, script_str)

        # Add DefinedIn column
        variable_dict = append_defined_script(variable_dict, script_suffix, script_str)

        # Add UsedInFunction and DefinedInFunction columns for .sh files
        if script_suffix.endswith(".sh"):
            variable_dict = append_function_info(variable_dict, script_str)

    # Convert to DataFrame and sort the rows
    df = pd.DataFrame.from_dict(variable_dict, orient="index")
    df["LowercaseName"] = df.index.str.lower()
    df = df.sort_values(by="LowercaseName").drop(columns="LowercaseName")

    # Convert to Markdown table
    header = ["ShellVariable"] + df.columns.tolist()
    markdown_table = tabulate(df, headers=header, tablefmt="pipe")

    # Print the Markdown table
    with open(output_file, "w") as file:
        file.write("# Shell Variable Manifest File\n")
        file.write(
            "This file is intended to automatically list shell script variables "
            + "for greater traceability of inherited variables. The file is "
            + "regenerated by running `$ python src/utilities/generate_var_manifest.py`"
            + " from the top level IMI directory. This happens automatically on new pull"
            + " requests via Github Actions\n\n"
        )
        file.write(
            "Note: This does not include variables defined in python scripts.\n\n"
        )
        file.write(markdown_table)
