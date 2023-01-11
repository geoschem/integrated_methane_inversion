Using the IMI Clustering Options
================================

Why use the clustering options?
-------------------------------
The main computational cost of the IMI is running the perturbation simulations necessary to 
construct the jacobian. This requires running a (GEOS-Chem) Jacobian simulation for each 
state vector element. The default state vector that is generated with the IMI has state 
vector elements in native resolution, meaning each element corresponds with a GEOS-Chem grid 
cell (.25 degree or .5 degree resolution). However, if your state vector has a sufficiently 
large number of elements this can limit the feasibility of running the IMI -- either due to
prohibitively high AWS costs or compute time. Clustering your state vector elements reduces 
the number of state vector elements by aggregating elements together. 

Using the IMI clustering config options
---------------------------------------
To enable the IMI clustering options in the imi config file set 
``ReducedDimensionStateVector: true``. This enables the clustering component of the IMI. 
Once enabled the IMI uses your specified ``ClusteringPairs`` to aggregate state vector elements 
within your domain of interest (excluding buffer elements). eg:

::

    ReducedDimensionStateVector: true
    ClusteringPairs:
      - [1, 15]
      - [2, 24]


Each clustering pair consists of the the aggregation level and the number of cells you are 
allocating with the aggregation level. In the above example, the user is requesting 15 native 
resolution state vector elements and 24 state vector elements to be aggregated with another 
element. Any additional elements that have not been allocated are then aggregated into a 
single element. Using the above clustering pairs, if the domain of interest has 63
elements in the original state vector, 15 of the elements would maintain the original resolution 
and 48 of the elements would be aggregated into 24 2-gridcell elements. If the original state 
vector has 75 elements in the domain of interest, then the remaining 12 unallocated elements are
aggregated into a single element, netting a new state vector with 40 elements in the domain of 
interest. If the specified cluster pairs request a greater number of elements to aggregate than 
are available in the original state vector the IMI will print an error to the output log.

Additionally, the ``ForcedNativeResolutionElements`` is a configuration option that allows you to
specify areas that you would like to maintain high resolution for. Eg:

::
    
    ForcedNativeResolutionElements:
      - [31.5, -104]


For instance, if the user suspects a location to be an emission hotspot they can specify the 
lat/lon coordinates as in the example above and the clustering algorithm will ensure that the
native resolution element is preserved during the aggregation. In order for the IMI to 
preserve the element, you must have enough native resolution pixels specified in your 
``ClusteringPairs``.

Note: The IMI preserves the original state vector file as NativeStateVector.nc in your run directory.

IMI clustering scheme
---------------------
The IMI clustering algorithm uses simplified version of the k-means method described 
`in Nesser et al., 2021 <https://doi.org/10.5194/amt-14-5521-2021>`_ to maintain native 
resolution in areas with high information content (high prior emissions, high observation 
density), while aggregating cells with low information content.

Reducing computational cost while maintaining inversion quality
---------------------------------------------------------------
While clustering is an effective method for alleviating computational constraints for 
running inversions at high resolution for large regions, it can introduce aggregation error
and degrade the quality of your inversion 
(`Turner and Jacob., 2014 <https://doi.org/10.5194/acp-15-7039-2015>`_ ). 
Therefore, it is important to weigh the computational benefits of reducing your state vector
against the inversion quality loss. This can be done by iteratively tuning the cluster
pairings and running the `IMI preview <../advanced/imi-preview.html>`__.IMI preview to assess 
the estimated DOFS. Ideally, you should find a middle groud where the estimated DOFS and 
computation cost is at a acceptable level before proceeding with the inversion.

