#!/usr/bin/env python3

"""
Description:
------------
This Python script (assumes Python3) reads a GEOS-Chem or
HEMCO-standalone log file containing dry-run output and does
the following:

    (1) Creates a list of unique files that are required for the
        GEOS-Chem or HEMCO-standalone simulation;

    (2) Creates a bash script to download missing files from the AWS
        s3://gcgrid bucket or from a specified server;

    (3) Executes the bash script to download the necessary data;

    (4) Removes the bash script upon successful download.


Remarks:
--------
    (1) This script only requires the "os", "sys", "subprocess", and
        PyYaml packages.

    (2) Jiawei Zhuang found that it is much faster to issue aws s3 cp
        commands from a bash script than a Python script.  Therefore,
        in this routine we create a bash script with all of the
        download commands that will be executed by the main routine.
"""

# Imports
import os
import sys
import subprocess
import yaml

# Exit with error if we are not using Python3
assert sys.version_info.major >= 3, \
"ERROR: Python 3 is required to run download_data.py!"

# Define global variables
DATA_DOWNLOAD_SCRIPT = "./auto_generated_download_script.sh"
CONFIG_FILE = "./HEMCO_sa_Config.rc"
TIME_FILE = "./HEMCO_sa_Time.rc"

def read_config_file(
        config_file,
        to_str=False
):
    """
    Reads configuration information from a YAML file.

    Args:
    -----
    config_file : str
        The configuration file in YAML format
    to_str : bool
        Set this to True if you wish to return the data in the YAML
        file as strings, or False otherwise.

    Returns:
    --------
    config : dict
        Dictionary with the contents of the YAML file
    """
    try:
        with open(config_file, encoding="UTF-8") as stream:
            if to_str:
                return yaml.load(stream, Loader=yaml.loader.BaseLoader)
            return yaml.load(stream, Loader=yaml.loader.SafeLoader)
    except FileNotFoundError as err:
        msg = f"Error reading configuration in {config_file}: {err}"
        raise FileNotFoundError(msg) from err


def extract_pathnames_from_log(
        args
):
    """
    Returns a list of pathnames from a GEOS-Chem log file.

    Args:
    -----
    args : dict
        Contains output from function parse_args.

    Returns:
    --------
    paths : dict
        paths["comments"]: Dry-run comment lines.
        paths["found"] : List of file paths found on disk.
        paths["missing"]: List of file paths that are missing.
        paths["local_prefix"]: Local data directory root.

    Author:
    -------
    Jiawei Zhuang (jiaweizhuang@g.harvard.edu)
    Modified by Bob Yantosca (yantosca@seas.harvard.edu)
    """

    # Initialization
    comments = ["!"*79,
                "!!! LIST OF (UNIQUE) FILES REQUIRED FOR THE SIMULATION"]
    data_found = set()
    data_missing = set()
    dryrun_log = args["dryrun_log"]

    # Open file (or die with error)
    with open(dryrun_log, "r", encoding="UTF-8") as ifile:

        # Read data from the file line by line.
        # Add file paths to the data_list set.
        line = ifile.readline()

        while line:

            # Convert line to uppercase for string match
            upcaseline = line.upper()

            # Search for data paths that have been found
            if (": OPENING" in upcaseline) or (": READING" in upcaseline):
                data_found.add(line.split()[-1])

            # Search for data paths that are missing
            elif "FILE NOT FOUND" in upcaseline:
                data_missing.add(line.split()[-1])

            # Search for certain dry-run comment strings
            # (and make sure to prevent duplicates)
            elif ("!!! STA" in upcaseline) or ("!!! END" in upcaseline) or \
                 ("!!! SIM" in upcaseline) or ("!!! MET" in upcaseline) or \
                 ("!!! GRI" in upcaseline):
                if line.rstrip() not in comments:
                    comments.append(line.rstrip())

            else:
                pass

            # Read next line
            line = ifile.readline()

        # Add another line to the comment list
        comments.append("!"*79)

        # Convert sets to lists and sort in alphabetical order
        found = sorted(list(data_found))
        missing = sorted(list(data_missing))

        # Find the local data directory prefix (path to ExtData)
        local_prefix = ""
        for path in found + missing:
            if "ExtData" in path:
                index = path.find("ExtData")
                local_prefix = path[:index]
                if "ExtData" not in local_prefix:
                    local_prefix = os.path.join(local_prefix, "ExtData")
                break

        # Exit if the local path does not contain ExtData
        if len(local_prefix) == 0:
            msg = \
                "Could not locate the ExtData folder in your local disk space!"
            raise ValueError(msg)

        # Close file and return
        # The "sorted" command will return unique values
        ifile.close()

        paths = {
            "comments": comments,
            "found": found,
            "missing": missing,
            "local_prefix": local_prefix
        }
        return paths


def get_run_info():
    """
    Searches through the geoschem_config.yml file for GEOS-Chem
    simulation parameters.

    Returns:
    -------
    run_info : dict
        Contains the GEOS-Chem run parameters: start_date,
        start_time, end_date, end_time, met, grid, and sim.
    """

    # Create dictionary with GEOS-Chem simulation parameters
    # NOTE: Numbers are returned as strings, and need to be converted
    run_info = {}

    # Read start & end datetimes
    config = read_config_file(TIME_FILE, to_str=True)
    temp = config["START"].replace("-","").replace(":","").split()
    run_info["start_date"] = int(temp[0])
    run_info["start_time"] = int(temp[1])
    temp = config["END"].replace("-","").replace(":","").split()
    run_info["end_date"] = int(temp[0])
    run_info["end_time"] = int(temp[1])

    # Get root data dir, met field, and resolution
    with open(CONFIG_FILE, "r", encoding="UTF-8") as ifile:
        for line in ifile:
            line = line.strip()
            if "END SECTION SETTINGS" in line:
                break
            if "ROOT:" in line:
                run_info["root_data_dir"] = line.split(":")[-1].strip(" ")
            if "MET:" in line:
                run_info["met_field"] = line.split(":")[-1].strip(" ")
            if "RES:" in line:
                run_info["resolution"] = line.split(":")[-1].strip(" ")
                run_info["grid"] = run_info["resolution"]

    return run_info


def replace_entry_in_list(
        the_list,
        old_entry,
        new_entry
):
    """
    Replaces a string entry in a list with a new entry.

    Args:
    -----
    the_list : list of str
       The list
    old_entry : (str
        Entry to replace
    new_entry : str
        Replacement text

    Returns:
    --------
    the_list : list of str
        The modified list
    """
    return list(map(lambda x: x.replace(old_entry, new_entry), the_list))


def write_unique_paths(
        paths,
        unique_log
):
    """
    Writes unique data paths from dry-run output to a file.

    Args:
    -----
        paths : dict
            Contains output from function extract_pathnames_from_log.

        unique_log : str
            Log file that will hold unique data paths.
    """
    combined_paths = paths["found"] + paths["missing"]
    combined_paths.sort()

    try:
        with open(unique_log, "w", encoding="UTF-8") as ofile:
            for comment in paths["comments"]:
                print(comment, file=ofile)
            for path in combined_paths:
                print(path, file=ofile)
            for comment in paths["comments"]:
                print(comment, file=ofile)
        ofile.close()
        print(f"Log with unique file paths written to: {unique_log}")
    except RuntimeError as exc:
        raise RuntimeError(f"Could not write {unique_log}") from exc


def create_download_script(
        paths,
        args
):
    """
    Creates a data download script to obtain missing files
    from the s3://gcgrid bucket on the AWS cloud or from a
    specified server.

    Args:
    -----
    paths : dict
        Contains output from function extract_pathnames_from_log.
    args : dict
        Contains output from function parse_args.
    """

    # Extract portal parameters
    portal_name = args["portal"]
    portal = args["config"]["portals"][portal_name]
    is_s3_bucket = portal["s3_bucket"]
    remote_root = portal["remote"]
    quote = portal["quote"]
    cmd_prefix = portal["command"]
    if "@PATH@" in cmd_prefix:
        cmd_prefix = cmd_prefix.replace("@PATH@", paths["local_prefix"])

    # Create the data download script
    with open(DATA_DOWNLOAD_SCRIPT, "w", encoding="UTF-8") as ofile:

        # Write shebang line to script
        print("#!/bin/bash\n", file=ofile)
        print("# This script was generated by download_data.py\n", file=ofile)

        # Write download commands for only the missing data files
        for path in paths["missing"]:

            # We do not download HEMCO restart files
            if "HEMCO_restart" in path:
                continue

            # Write path
            index = path.find("ExtData") + 7
            local_dir = os.path.dirname(path)
            remote_path = remote_root + path[index:]
            cmd = cmd_prefix + quote + remote_path + quote
            if is_s3_bucket:
                cmd += " " + local_dir + "/"
            print(cmd, file=ofile)
            print(file=ofile)

        # Close file and make it executable
        ofile.close()
        os.chmod(DATA_DOWNLOAD_SCRIPT, 0o755)

def download_the_data(
        args
):
    """
    Downloads GEOS-Chem data files from the AWS s3://gcgrid bucket
    or from a specified server.

    Args:
    -----
    args : dict
        Output of runction parse_args.
    """

    # Get information about the run
    run_info = get_run_info()

    # Get a unique list of data paths, both found and missing:
    # Expand the data paths to include links to restart files
    paths = extract_pathnames_from_log(args)

    # Write a list of unique file paths
    write_unique_paths(paths, args["dryrun_log"] + ".unique")

    # Exit without downloading if skip-download flag was specified
    if args["skip_download"]:
        return

    # Print a message
    if len(args["portal"]) > 0:
        print(f"Downloading data from {args['portal']}")

    # Create script to download missing files from AWS S3
    create_download_script(paths, args)

    #### DEBUG: Uncomment this if you wish to see the download script
    #if args["skip_download"]:
    #    return

    # Run the data download script and return the status
    # Remove the file afterwards
    status = subprocess.call(DATA_DOWNLOAD_SCRIPT)
    os.remove(DATA_DOWNLOAD_SCRIPT)

    # Raise an exception if the data was not successfully downloaded
    if status != 0:
        msg = f"Error downloading data from {args['portal']}"
        raise RuntimeError(msg)


def parse_args():
    """
    Reads global settings from the download_data.yml configuration file.
    Also parses command-line arguments and returns a dictionary
    containing all of these settings.

    Returns:
    --------
    args : dict
        args["config"] : Dict with global settings from download_data.yml
        args["dryrun_log"] Name of the GEOS-Chem dry-run log file
        args["portal"]: Name of the remote portal for download
        args["skip_download"]: Are we skipping the download? (T/F)
    """
    dryrun_log = None
    dryrun_found = False
    portal_found = False
    portal_remote = None
    skip_download = False
    skip_found = False

    # Read the YAML configuration file
    config = read_config_file("download_data.yml")

    # Get a list of portal names + short names
    portal_list = list(config["portals"].keys())
    short_name_list = []
    for mir in portal_list:
        short_name_list.append(config["portals"][mir]["short_name"])

    # Parse command-line arguments (argument 0 is the program name)
    for i in range(1, len(sys.argv)):
        arg = sys.argv[i].lower()
        arg = arg.lstrip('-')

        if not dryrun_found:
            dryrun_log = arg
            dryrun_found = True
            continue

        if not portal_found:
            for mir in portal_list:
                portal = mir.lower()
                short_name = config["portals"][mir]["short_name"].lower()
                if arg in portal or arg in short_name:
                    portal_remote = portal
                    portal_found = True
                    continue

        if not skip_found:
            if "skip" in arg:
                skip_download = True
                skip_found = True
                continue


    if dryrun_log is None:
        msg = "The dryrun log file was not supplied!  Exiting ..."
        raise ValueError(msg)

    if portal_remote is None and not skip_download:
        msg = "Portal name missing or invalid!  Exiting ..."
        raise ValueError(msg)

    args = {
        "config": config,
        "dryrun_log": dryrun_log,
        "portal": portal_remote,
        "skip_download": skip_download
    }
    return args


def main():
    """
    Main program.  Gets command-line arguments and calls function
    download_the_data to initiate a data-downloading process.

    Calling sequence:
    -----------------
        ./download_data.py log PORTAL-NAME
        ./download_data.py log -skip-download  # Print unique log & exit
    """

    # Download the data files from the remote server
    download_the_data(parse_args())


if __name__ == "__main__":
    main()
