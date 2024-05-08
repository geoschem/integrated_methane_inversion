import xarray as xr
import os


def multiply_posteriors(period_number, base_directory, lognormal):
    """
    Before running the posterior simulation to update initial conditions for the next period, need to
    apply the latest posterior scale factors to the dynamic ScaleFactors.nc file. The posterior
    simulation will use those updated scale factors.

    Arguments
        period_number   [int] : What period are we on? For the first period, period_number = 1
        base_directory  [str] : The base directory for the inversion, where e.g., "preview_sim/" resides

    """

    # Useful paths
    sf_path = os.path.join(base_directory, "ScaleFactors.nc")
    posterior_dir = os.path.join(base_directory, f"kf_inversions/period{period_number}")
    gridded_posterior_filename = (
        "gridded_posterior_ln.nc" if lognormal else "gridded_posterior.nc"
    )
    gridded_posterior_path = os.path.join(posterior_dir, gridded_posterior_filename)

    # Load newest gridded posterior and previous/final gridded posterior
    sf = xr.load_dataset(sf_path)
    latest_posterior = xr.load_dataset(gridded_posterior_path)

    # Multiply new and old posteriors
    sf["ScaleFactor"] = sf["ScaleFactor"] * latest_posterior["ScaleFactor"]

    # Ensure good netcdf attributes for HEMCO
    sf.lat.attrs["units"] = "degrees_north"
    sf.lat.attrs["long_name"] = "Latitude"
    sf.lon.attrs["units"] = "degrees_east"
    sf.lon.attrs["long_name"] = "Longitude"
    sf.ScaleFactor.attrs["units"] = "1"

    # Save final posterior
    save_path = os.path.join(base_directory, "ScaleFactors.nc")
    sf.to_netcdf(
        save_path,
        encoding={v: {"zlib": True, "complevel": 9} for v in sf.data_vars},
    )

    # Archive scale factors
    archive_path = os.path.join(base_directory, "archive_sf")
    sf.to_netcdf(
        os.path.join(archive_path, f"posterior_sf_period{period_number}.nc"),
        encoding={v: {"zlib": True, "complevel": 9} for v in sf.data_vars},
    )


if __name__ == "__main__":
    import sys

    period_number = sys.argv[1]
    base_directory = sys.argv[2]
    lognormal = sys.argv[3] == "true"

    multiply_posteriors(period_number, base_directory, lognormal)
