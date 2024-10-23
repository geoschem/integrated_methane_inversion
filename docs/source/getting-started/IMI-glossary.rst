IMI Glossary
============
Glossary of commonly used terms in the Integrated Methane Inversion (IMI) workflow.


.. list-table::
   :widths: 30, 70
   :class: tight-table

   * - ``Period of interest``
     - The period for which the IMI will optimize mean emissions based on TROPOMI observations. Specified by start and end dates. Only observations made during the period are considered. The period of interest can be split into shorter sub-periods (to be optimized sequentially) through the IMI Kalman filter feature.
   * - ``Region of interest``
     - The region over which the IMI will optimize mean emissions at up to 0.25°×0.3125° (≈25-km) resolution. Specified by rectilinear latitude/longitude bounds, shapefile, or interactively through the Integral Earth user interface. The region of interest can be rectilinear or irregular in shape.
   * - ``Inversion domain``
     - The region of interest and a surrounding buffer region. The inversion domain is always rectilinear in shape.
   * - ``Buffer emission elements``
     - The 2D emission elements that make the buffer region. Default number is 8.
   * - ``Emission state vector``
     - The collection of 2D emission elements (up to 0.25°×0.3125° resolution) to be optimized in the inversion. Includes elements within the region of interest and buffer elements.
   * - ``Prior emission estimate``
     - Best estimate of emissions before performing the inversion, based on a bottom-up inventory.
   * - ``Posterior emission estimate``
     - Best estimate of emissions after performing the inversion.
   * - ``Averaging kernel sensitivity``
     - Estimates how sensitive the posterior solution for a given state vector (emission) element is to observations as opposed to the prior estimate. An emission element with averaging kernel sensitivity 0 is not quantified by the observations at all, and the inversion returns the prior value for that element. An emission element with averaging kernel sensitivity 1 is fully quantified by the observations, and the inversion results for that element are independent of the prior estimate. 
   * - ``Degrees of freedom for signal (DOFS)``
     - The sum of the averaging kernel sensitivities for the region of interest. Measures the information content of the observations towards optimizing the state vector; represents the number of independent pieces of information on the state vector that the observations can quantify.
