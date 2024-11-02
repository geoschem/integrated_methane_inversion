import os
import sys
import yaml

"""
A simple utility script that replaces imi config file variables with the 
corresponding environment variables if they exist. The resultant config 
file is written to the desired path. environment variables have the format
IMI_<VARIABLE_NAME>.
the 
Arguments
    config_path        [String] : path to yaml config file
    output_config_path [String] : path to write output config file
"""

if __name__ == "__main__":
    config_path = sys.argv[1]
    config_ouput_path = sys.argv[2]

    # load config file
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    print(f"Overriding config variables in {config_path} with environment variables.")
    override_count = 0
    for key in config.keys():
        env_var = f"IMI_{key}"
        # check if environment variable existss
        if env_var in os.environ:
            print(f"Updating '{key}: {config[key]}' with '{os.environ[env_var]}'")
            updated_value = os.environ[env_var].strip("'")
            updated_value = os.environ[env_var].strip('"')
            config[key] = updated_value
            override_count += 1

    # write overwritten output config file to desired path
    if override_count > 0:
        print(f"Writing output config file to {config_ouput_path}")
        with open(config_ouput_path, "w") as f:
            yaml.safe_dump(config, f)
    else:
        print(f"No config variables overridden.")
