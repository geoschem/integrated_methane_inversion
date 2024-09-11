#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import yaml
import sys

# Description: Parse yaml file and convert to shell variables
# Example Usage as a script:
#   $ python parse_yaml.py config.yml
# To actually export returned shell variables:
#   $ eval $(python parse_yaml.py config.yml)


def parse_yaml(file_path, prefix=""):
    """
    Parse yaml file and convert to shell variable format
    Arguments
        file_path [str] : path to yaml file
        prefix    [str] : prefix for shell variables
    """
    with open(file_path, "r") as stream:
        data = yaml.safe_load(stream)
        return yaml_to_shell_variables(data, prefix)


def yaml_to_shell_variables(data, prefix=""):
    """
    Traverse dictionary data and convert to shell variables
    Arguments
        data   [dict] : yaml file data dictionary
        prefix [str]  : prefix for shell variables
    """
    shell_vars = []

    def recurse(data, current_prefix):
        if isinstance(data, dict):
            for key, value in data.items():
                recurse(value, f"{current_prefix}{key}_")
        elif isinstance(data, list):
            array_elements = []
            for idx, value in enumerate(data):
                if isinstance(value, dict) or isinstance(value, list):
                    recurse(value, f"{current_prefix}{idx}_")
                else:
                    array_elements.append(f'"{value}"')
            if array_elements:
                array_declaration = (
                    f"{current_prefix[:-1]}=({' '.join(array_elements)})"
                )
                shell_vars.append(array_declaration)
        elif isinstance(data, bool):
            shell_vars.append(f'{current_prefix[:-1]}="{str(data).lower()}"')
        else:
            shell_vars.append(f'{current_prefix[:-1]}="{data}"')

    recurse(data, prefix)
    return shell_vars


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 parse_yaml.py <path_to_yaml>")
        sys.exit(1)

    yaml_file = sys.argv[1]
    shell_vars = parse_yaml(yaml_file)

    for var in shell_vars:
        print(var)
