import os
import sys
import yaml


def change_dates(new_start_date, new_end_date, directory_to_edit):
    """
    Replace geoschem_config.yml start and end dates in selected directory,
    including all GEOS-Chem run directories if a parent directory

    Arguments
        new_start_date    [str] : Desired start date (yyyymmdd)
        new_end_date      [str] : Desired end date (yyyymmdd)
        directory_to_edit [str] : Directory containing GEOS-Chem run directories
    """

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
    new_start_date = sys.argv[1]
    new_end_date = sys.argv[2]
    directory_to_edit = sys.argv[3]

    # Run the script
    change_dates(new_start_date, new_end_date, directory_to_edit)
