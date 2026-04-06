import os
import sys
import yaml
import re

def change_dates(config_fpath, new_start_date, new_end_date, run_duration, directory_to_edit):
    """
    Replace geoschem_config.yml start and end dates in selected directory,
    including all GEOS-Chem run directories if a parent directory

    Arguments
        config_fpath      [str] : path for config.yml
        new_start_date    [str] : Desired start date (yyyymmdd)
        new_end_date      [str] : Desired end date (yyyymmdd)
        run_duration      [str] : Run duration (yyyymmdd)
        directory_to_edit [str] : Directory containing GEOS-Chem run directories
    """

    config = yaml.load(open(config_fpath), Loader=yaml.FullLoader)
    # Define desired geoschem_config.yml date strings
    new_start_line = f"  start_date: [{new_start_date}, 000000]\n"
    new_end_line = f"  end_date: [{new_end_date}, 000000]\n"

    # Determine run directories to edit
    contents = os.listdir(directory_to_edit)
    if "geoschem_config.yml" in contents:
        rundirs = [directory_to_edit]
    else:
        rundirs = [
            os.path.join(directory_to_edit, r)
            for r in contents
            if os.path.isdir(os.path.join(directory_to_edit, r))
        ]

    # Process them
    for d in rundirs:
        if config['UseGCHP']:
            cap_restart_path = os.path.join(d, "cap_restart")
            with open(cap_restart_path, "w") as f:
                f.write(f"{new_start_date} 000000\n")
            
            set_common_path = os.path.join(d, "setCommonRunSettings.sh")
            with open(set_common_path, "r") as f:
                content = f.read()

            # Replace Run_Duration="YYYYMMDD 000000"
            new_duration_line = f'Run_Duration="{run_duration} 000000"'
            content_new = re.sub(r'Run_Duration="\d{8} 000000"', new_duration_line, content)

            with open(set_common_path, "w") as f:
                f.write(content_new)
        else:
            pth = os.path.join(d, "geoschem_config.yml")
            with open(pth, "r") as file:
                filedata = file.read()
            with open(pth, "r") as file:
                for line in file:
                    if line.startswith("  start_date:"):
                        old_start_line = line
                    if line.startswith("  end_date:"):
                        old_end_line = line
            # Replace the target string
            filedata = filedata.replace(old_start_line, new_start_line)
            filedata = filedata.replace(old_end_line, new_end_line)
            # Write the file out again
            with open(pth, "w") as file:
                file.write(filedata)


if __name__ == "__main__":

    # Inputs
    config_fpath = sys.argv[1]
    new_start_date = sys.argv[2]
    new_end_date = sys.argv[3]
    run_duration = sys.argv[4]
    directory_to_edit = sys.argv[5]

    # Run the script
    change_dates(config_fpath, new_start_date, new_end_date, run_duration, directory_to_edit)
