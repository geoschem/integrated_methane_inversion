IMI FAQs
========

This page documents frequently asked questions about the IMI.

What is the IMI?
~~~~~~~~~~~~~~~~~
* The IMI is an open-source software tool for quantifying methane emissions at up to 0.25°x0.3125° (≈25-km) and weekly resolution using satellite observations from the `TROPOspheric Monitoring Instrument <https://www.tropomi.eu/>`_ (TROPOMI), a prior estimate of emissions (e.g., a bottom-up emission inventory), and the `GEOS-Chem chemical transport model <https://geoschem.github.io/index.html>`_.
* It uses an analytical Bayesian inversion method that starts from the prior estimate of emissions and improves it with information from the satellite observations to produce a posterior estimate with full error characterization.
* It can be applied to any region and period of interest.

How do I access the IMI?
~~~~~~~~~~~~~~~~~~~~~~~~~
* There are several ways to access the IMI:
    * Use the free IMI product on the `Amazon Web Services (AWS) Marketplace <https://aws.amazon.com/marketplace/pp/prodview-hkuxx4h2vpjba>`_.
    * Download the `IMI source code <https://github.com/geoschem/integrated_methane_inversion>`_ and run it locally.
    * Use the `Integral Earth web user interface <https://integralearth.github.io/>`_ for the IMI.

How do I cite the IMI?
~~~~~~~~~~~~~~~~~~~~~~~
* You can cite the IMI with the corresponding research paper for your application: `IMI 1.0 paper <https://doi.org/10.5194/gmd-15-5787-2022>`_ and/or `IMI 2.0 paper <https://doi.org/10.5194/egusphere-2024-2700>`_.

How much does it cost to use the IMI?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Running the IMI on AWS incurs fees for using AWS compute resources. 
* The typical cost for a 500x500 km\ :sup:`2` domain (e.g., the Permian Basin) for 1 month is about $20. 
* Cost scales with duration and domain size. Costs for larger domains can be effectively mitigated using the smart clustering capability available through the IMI.

Where is the IMI documented?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The IMI is documented at `imi.readthedocs.io <../index.html>`_.

What is the IMI Preview?
~~~~~~~~~~~~~~~~~~~~~~~~~
* The `IMI Preview <../getting-started/imi-preview.html>`_ is a feature for evaluating an IMI configuration without actually running an inversion. With the IMI Preview you will:
* Visualize the TROPOMI observations, bottom-up emission inventories, and point source data to be used in the inversion. 
* Estimate the information content (degrees of freedom for signal, called DOFS) of the inversion.
* Estimate the USD cost of running the inversion on AWS.
* The IMI Preview has no significant costs, and we strongly recommend using it to ensure that the proposed IMI configuration will lead to a successful inversion.

Does the IMI support continuous emission monitoring?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Yes. Users can continuously monitor emissions for a region of interest using the `IMI Kalman filter <../advanced/kalman-filter-mode.html>`_ feature.

Does the IMI support use of custom prior emission inventories?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Yes. Users can override the default IMI prior emission inventories with their own by following `these instructions <../advanced/custom-prior-emissions-hemco.html>`_.

Can I get information on individual point sources from the IMI?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* No. The IMI gives you total emissions in 25x25 km\ :sup:`2` grid cells, including all sources. However, it scrapes point source observations from other databases to improve the prior estimate (optional) and to put the total emission results in context.

Can I get sectoral information from the IMI?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Yes. The IMI uses prior sectoral information to allocate the posterior emission estimates to specific sectors.

Can I detect offshore emissions with the IMI?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Yes. Although the satellite observations are mainly over land, the IMI has the option to use glint observations over water. Offshore emissions can also be quantified from observations of the plume advected over land by onshore flow.

What can I get from the IMI that's different from point source data providers (GHGSat, Carbon Mapper, Kayrros...)?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* The IMI is a completely different product. It provides total gridded continuous emissions, not snapshot emissions from specific point sources as from point source data providers. 
* The two are complementary. The IMI information is most useful for emission reporting, for understanding contributions from different sectors, for monitoring emission trends, and for quantifying long-term averages.


