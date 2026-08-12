.. _imi_glossary:

Terminology
==================================

.. glossary::
   :sorted:

   Period of interest
      The period for which the IMI will optimize mean emissions based on TROPOMI observations. Specified by start and end dates. Only observations made during the period are considered. The period of interest can be split into shorter sub-periods (to be optimized sequentially) through the IMI Kalman filter feature.

   Region of interest
      The region over which the IMI will optimize mean emissions at up to 0.125°×0.15625° (≈12km\ :sup:`2`) resolution. Specified by rectilinear latitude/longitude bounds or a shapefile.

   Inversion domain
      The region of interest and a surrounding buffer region. The inversion domain is always rectilinear in shape.
   
   State vector
      The collection of 2D elements to be optimized in the inversion. This includes 2D emission elements (at up to 0.125°×0.15625° resolution) within the region of interest and buffer elements. It may also contain OH and/or boundary conditions.

   Buffer elements
      The 2D emission elements that make the buffer region. Default number is 8.

   Prior emission estimate
      Best estimate of emissions before performing the inversion, based on bottom-up inventories.

   Posterior emission estimate
      Best estimate of emissions after performing the inversion.

   Averaging kernel sensitivity
      Estimates how sensitive the posterior solution for a given state vector (emission) element is to observations as opposed to the prior estimate. An emission element with averaging kernel sensitivity 0 is not quantified by the observations at all, and the inversion returns the prior value for that element. An emission element with averaging kernel sensitivity 1 is fully quantified by the observations, and the inversion results for that element are independent of the prior estimate. 

   Degrees of freedom for signal (DOFS)
      The sum of the averaging kernel sensitivities for the region of interest. Measures the information content of the observations towards optimizing the state vector; represents the number of independent pieces of information on the state vector that the observations can quantify.
