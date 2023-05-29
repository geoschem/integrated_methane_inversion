import re
import glob
import sys

def check_for_duplicate_orbit_numbers(Sat_datadir):

    """
    Function that checks for duplicate TROPOMI filenames in your directory (either on your cluster or after download to AWS)
    Takes advatnage of the fact that there should be the same number of unique orbit numbers as there are files
    """
    
    files = sorted(glob.glob(Sat_datadir + "/*.nc"))

    all_orbit_numbers = []
    for filename in files:
        matches = re.findall(r'_(\d{5})_', filename) # find orbit number
        assert len(matches) == 1, f"Please check the TROPOMI filenames, as they are not formatted in the way that the IMI expects. (e.g., {filename})"
        all_orbit_numbers.append(matches[0])
        
    number_of_unique_orbit_numbers = len(set(all_orbit_numbers)) # forming a set drops the duplicates
    number_of_files_in_Sat_datadir = len(files)

    assert number_of_unique_orbit_numbers == number_of_files_in_Sat_datadir, "Duplicate orbit numbers found in TROPOMI datafiles. Remove duplicate files from system before continuing."

if __name__ == "__main__":
    Sat_datadir = sys.argv[1]
    check_for_duplicate_orbit_numbers(Sat_datadir)