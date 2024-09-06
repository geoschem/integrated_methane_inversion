import numpy as np
import datetime
import os


def cull_restarts(jacobian_runs_dir, posterior_dir, start_day, end_day):
    """
    GEOS-Chem bugs mean we need to output daily restart files instead of just one
    at the end of each week. This function culls the extra restarts we don't need.

    Arguments
        jacobian_runs_dir [str] : directory containing Jacobian run directories
        posterior_dir     [str] : posterior run directory
        start_day         [str] : first day of current 1-week inversion ("yyyymmdd")
        end_day           [str] : last day of current 1-week inversion ("yyyymmdd")
    """

    # List Jacobian run directories
    contents = os.listdir(jacobian_runs_dir)
    jacobian_dirs = [r for r in contents if "CH4_" in r]
    run_dirs = [os.path.join(jacobian_runs_dir, r) for r in jacobian_dirs]

    # Add posterior run directory to the list
    run_dirs.append(posterior_dir)

    # Datetimes for the start and end days
    start = datetime.datetime.strptime(start_day, "%Y%m%d")
    end = datetime.datetime.strptime(end_day, "%Y%m%d")

    # Look in each dir and remove extra restarts
    for direc in run_dirs:
        contents_i = os.listdir(direc)
        restarts = [r for r in contents_i if "GEOSChem.Restart" in r]
        for restart in restarts:
            date = restart[17:25]
            date_dt = datetime.datetime.strptime(date, "%Y%m%d")
            if date_dt > start and date_dt < end:
                delete_path = os.path.join(direc, restart)
                os.system(f"rm {delete_path}")

    print("Culled extra daily restart files")


if __name__ == "__main__":
    import sys

    jacobian_runs_dir = sys.argv[1]
    posterior_dir = sys.argv[2]
    start_day = sys.argv[3]
    end_day = sys.argv[4]

    cull_restarts(jacobian_runs_dir, posterior_dir, start_day, end_day)
