Constructing an inversion ensemble
==================================

You can use the IMI to create a low-cost ensemble of sensitivity inversions with different 
inversion parameters. This is because the Jacobian matrix computed in the first inversion 
can easily be reused for additional  inversions. This can be useful for understanding the 
sensitivity of the inversion results to the choice of prior error, observation error, 
or other inversion parameters.

The IMI has multiple options for creating an inversion ensemble.

1. **Automated ensemble generation**:

The simplest way to generate an ensemble is to run the IMI a single time with
a configuration file that specifies vectors of the desired range of hyperparameters (eg.`PriorError: [0.5, 0.75]`). 
The IMI will then run multiple inversions with the various cominations of hyperparameter values. Each ensemble member is
saved to the inversion_results.nc and gridded_posterior.nc. This method is useful for quickly generating an ensemble without
having to manually run the IMI multiple times with new run directories and configuration files. However, vectors can only be
applied for the following hyperparameters: `PriorError`, `ObsError`, `Gamma`, `PriorErrorBCs`, `PriorErrorBufferElements`, 
`PriorErrorOH`.

In the result files (inversion_results.nc and gridded_posterior.nc), the data variables from the different
ensemble members are distinguished by the data variable suffix(eg. `xhat_1_ensemble_member`, `xhat_2_ensemble_member`). 
The hyperparameters used for each ensemble member are saved in the attributes of each data variable 
(eg. `ds[xhat_1_ensemble_member].attrs`). The data variables without a suffix are from the ensemble member that most 
closely matches the expected output of the chi-square distribution (:math:`J_a / n \approx 1`) (Lu et al., 2021), 
where :math:`n` is the number of state vector elements.

2. **Manual ensemble generation**:
In this scenario, you have already run your base inversion. You can then use the IMI to create an ensemble of inversions
by specifying a new run directory that references the precomputed Jacobian matrix from the base run directory.
This method is useful for generating an ensemble with inputs that cannot be specified as vectors in the configuration
file. For example, you may want to run an ensemble of sensitivity inversions that swap out prior emission inventories, 
or use Lognormal instead of Gaussian prior error distributions. Note: When applying state vector clustering you cannot 
use the precomputed Jacobian if you swap out prior emissions inventories. This is because the underlying distribution of 
emissions within individual state vector clusters may change, requiring a new Jacobian matrix to be computed.

See the `Common configurations page <../other/common-configurations.html#running-a-sensitivity-inversion>`__ 
for instructions on how to re-configure the IMI to use a pre-computed Jacobian in a new run directory. Then modify
the values of ``PriorError``, ``ObsError``, and/or ``Gamma`` in the configuration file and re-run the inversion.
