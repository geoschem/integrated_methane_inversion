import xarray as xr
import datetime
from joblib import Parallel, delayed
from utils import zero_pad_num_hour


def setup_gc_cache(startday, endday, gc_source_path, gc_destination_path):
    """
    This script sets up a directory containing hourly GEOS-Chem output diagnostics
    files. The hourly files are convenient for computing virtual TROPOMI columns
    from the GEOS-Chem simulated atmosphere (to compare with the real TROPOMI columns).

    Arguments
        startday            [str] : First day of inversion period; formatted YYYYMMDD
        endday              [str] : Last day of inversion period; formatted YYYYMMDD
        gc_source_path      [str] : GEOS-Chem output directory
        gc_destination_path [str] : Target GEOS-Chem data directory in inversion workspace

    """

    # Make date range
    days = []
    dt = datetime.datetime.strptime(startday, "%Y%m%d")
    dt_max = datetime.datetime.strptime(endday, "%Y%m%d")
    while dt < dt_max:
        dt_str = str(dt)[0:10].replace("-", "")
        days.append(dt_str)
        delta = datetime.timedelta(days=1)
        dt += delta

    hours = range(24)

    # For each day:
    def process(d):
        # Load the SpeciesConc and LevelEdgeDiags data
        SpeciesConc_data = xr.load_dataset(
            f"{gc_source_path}/GEOSChem.SpeciesConc.{d}_0000z.nc4"
        )
        LevelEdgeDiags_data = xr.load_dataset(
            f"{gc_source_path}/GEOSChem.LevelEdgeDiags.{d}_0000z.nc4"
        )

        # For each hour:
        for h in hours:

            # Select data for that hour
            SpeciesConc_for_hour = SpeciesConc_data.isel(time=slice(h, h + 1, 1))
            LevelEdgeDiags_for_hour = LevelEdgeDiags_data.isel(time=slice(h, h + 1, 1))

            # Save to new .nc4 file at destination
            SpeciesConc_save_pth = f"{gc_destination_path}/GEOSChem.SpeciesConc.{d}_{zero_pad_num_hour(h)}00z.nc4"
            LevelEdgeDiags_save_pth = f"{gc_destination_path}/GEOSChem.LevelEdgeDiags.{d}_{zero_pad_num_hour(h)}00z.nc4"
            SpeciesConc_for_hour.to_netcdf(
                SpeciesConc_save_pth,
                encoding={v: {"zlib": True, "complevel": 1} for v in SpeciesConc_for_hour.data_vars},
            )
            LevelEdgeDiags_for_hour.to_netcdf(
                LevelEdgeDiags_save_pth,
                encoding={v: {"zlib": True, "complevel": 1} for v in LevelEdgeDiags_for_hour.data_vars},
            )

    results = Parallel(n_jobs=-1)(delayed(process)(day) for day in days)
    print(f"Set up hourly data files in {gc_destination_path}")


if __name__ == "__main__":
    import sys

    startday = sys.argv[1]
    endday = sys.argv[2]
    gc_source_path = sys.argv[3]
    gc_destination_path = sys.argv[4]

    setup_gc_cache(startday, endday, gc_source_path, gc_destination_path)
