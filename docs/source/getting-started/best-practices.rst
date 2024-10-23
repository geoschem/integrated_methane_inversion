IMI Best Practices
===================

Choosing an inversion time period
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The IMI can be applied to any period of interest beginning 1 May 2018, when the TROPOMI methane record begins.
* Common choices for the length of the inversion period are one year, one season (~3-6 months), one month, or one week.
* We recommend choosing time periods of one week or more to ensure there are enough satellite observations for a successful inversion.
* The `IMI Preview feature <../getting-started/imi-preview.html>`_ can be used to refine the choice of inversion period.

Defining a region of interest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The IMI can be applied to any region of interest, from the global scale down to small focus areas such as cities, oil and gas basins, and agricultural areas.
* The region of interest can be specified in several ways:
* Setting latitude/longitude bounds for a rectangular domain.
* Using a shapefile.
* Interactively in the Integral Earth web user interface.
* We recommend users select regions of interest larger than about 10,000 km\ :sup:`2` (100x100 km\ :sup:`2`) to ensure there are enough satellite observations for a successful inversion.
* Larger regions of interest require more computational resources. This can be mitigated by optimally reducing the effective resolution of the inversion via `smart state vector clustering <../advanced/using-clustering-options.html>`_.

Configuring the inversion domain
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Regional inversions focus on a region of interest within a larger rectilinear inversion domain.
* The inversion domain includes both the region of interest and an external buffer region.
* The buffer region is broken into a collection of buffer emission elements representing emissions outside the region of interest.
* We recommend using ≥ 8 buffer elements to pad the region of interest by ≥ 2°. The default number is 8.

Reducing the dimension of the state vector for large regions of interest
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Inversions for large regions of interest at the IMI native 0.25°x0.3125° grid resolution can be computationally expensive.
* This can be mitigated by reducing the dimension of the state vector using the state vector clustering options.
* Smart state vector clustering combines 0.25°x0.3125° into coarser grid elements where the prior emission estimates are low and/or where TROPOMI provides few observations

Interpreting the IMI Preview
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Examine the expected information content for the region and period of interest. This includes the map of expected averaging kernel sensitivities and the expected degrees of freedom for signal (DOFS).
    * The averaging kernel sensitivities should be higher where the prior emission estimates are higher and where more observations are available.
    * DOFS > 0.5 is a bare minimum to achieve any solid information about emissions.
    * DOFS < 2 is marginal for most applications.
* If the expected information content is low, consider:
    * Increasing the inversion period to incorporate more observations.
    * Increasing the prior error estimate.

Choosing the TROPOMI data product for the inversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The IMI supports inversions with two versions of the TROPOMI methane record:
    * The operational TROPOMI retrieval product developed by the SRON Netherlands Institute for Space Research.
    * The Blended TROPOMI+GOSAT retrieval product developed by `Balasus et al. (2023) <https://amt.copernicus.org/articles/16/3787/2023/>`_ to mitigate retrieval artifacts in the operational product.
* Choosing a product depends on the application. The operational product is updated every few days. The blended product is updated intermittently and is currently available through 2023.
* We recommend using the blended product when available (currently until 2024-01-01) to mitigate retrieval artifacts.
