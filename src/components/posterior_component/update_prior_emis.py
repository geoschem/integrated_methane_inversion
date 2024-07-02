import sys
import xarray as xr

def get_posterior_emissions(prior, scale):
    """
    Function to calculate the posterior emissions from the prior 
    and the scale factors. Properly accounting for no optimization 
    of the soil sink.
    Args:
        prior  : xarray dataset
            prior emissions
        scales : xarray dataset scale factors
    Returns:
        posterior : xarray dataset
            posterior emissions
    """
    # keep attributes of data even when arithmetic operations applied
    xr.set_options(keep_attrs=True)
    
    # we do not optimize soil absorbtion in the inversion. This 
    # means that we need to keep the soil sink constant and properly 
    # account for it in the posterior emissions calculation.
    # To do this, we:
    
    # make a copy of the original soil sink
    prior_soil_sink = prior["EmisCH4_SoilAbsorb"].copy()
    
    # remove the soil sink from the prior total before applying scale factors
    prior["EmisCH4_Total"] = prior["EmisCH4_Total"] - prior_soil_sink
    
    # scale the prior emissions for all sectors using the scale factors
    posterior = prior * scale["ScaleFactor"]
    
    # But reset the soil sink to the original value
    posterior["EmisCH4_SoilAbsorb"] = prior_soil_sink
    
    # Add the original soil sink back to the total emissions
    posterior["EmisCH4_Total"] = posterior["EmisCH4_Total"] + prior_soil_sink
    return posterior

if __name__ == "__main__":
    prior = xr.load_dataset(sys.argv[1])
    scale = xr.load_dataset(sys.argv[2])
    target_file = sys.argv[3]
    updated_emis = get_posterior_emissions(prior, scale)
    updated_emis.to_netcdf(target_file)