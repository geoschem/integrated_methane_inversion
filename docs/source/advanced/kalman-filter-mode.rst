
================================
The IMI Kalman Filter mode
================================

What is a Kalman Filter Inversion?
----------------------------------
A Kalman filter is a mathematical algorithm, developed by Rudolph Kalman, that estimates the state of a system by combining measurements and predictions while considering uncertainties. It operates recursively, continuously updating its estimate of the system state based on new measurements.

Kalman filters can be applied in atmospheric inversions by dividing an inversion period is into smaller time intervals, such as weekly chunks. An inversion is sequentially run for each interval, estimating the emissions for that specific period based on measurements and predictions. The resulting optimized emissions are then used as prior emissions for the next interval, allowing the prior emissions of each successive week to be informed by the previous weeks.

.. image:: img/kalman_filter.png
    :width: 500px
    :align: center
    :alt: Kalman filter diagram

Why use Kalman Filter Mode?
---------------------------
This approach enables tracking of how emissions change over time and provides insights into their distribution throughout the inversion period. By using the Kalman filter mode in the inversion, users can calculate intermediate emissions at the desired update frequency, such as weekly, revealing the temporal evolution of emissions.

How to use the Kalman Filter mode
=================================
The IMI Kalman Mode can be applied simpy by updating the ``KalmanMode`` config variable to ``true``. This will enable the Kalman filter mode using the specified update frequency, nudge factor, and first period.
Example Kalman filter config variables:

::

    ## Kalman filter options
    KalmanMode: true
    UpdateFreqDays: 7
    NudgeFactor: 0.1
    FirstPeriod: 1
      

The update frequency (``UpdateFreqDays``) is the number of days to for each inversion interval in the kalman filter.

FirstPeriod
-----------
The ``FirstPeriod`` config variable allows a user to select which chunked interval they would like the Kalman Filter to start on. This is most useful if you have a number of periods succeed eg. 5 out of 8 inversion periods succeed, and you would like to start the Kalman Filter on the 6th period. The ``FirstPeriod`` variable is set to 1 by default, which means the Kalman Filter will start on the first inversion period. If you would like to start the Kalman Filter on the 6th period, you would set ``FirstPeriod`` to 6. The ``FirstPeriod`` variable is a convenience variable, and is not required to run the Kalman Filter mode. If you do not specify a ``FirstPeriod``, the Kalman Filter will start on the first inversion period by default.

The NudgeFactor
---------------
A true Kalman Filter would use the posterior emissions from the previous interval as the prior emissions for the next interval. However, in practice, a direct substitution of the posterior emissions as the prior in the subsequent interval can lead to some emission elements getting locked at very low values. Retaining some information from the prior emissions can help to avoid this issue (`Varon et al., 2023 <https://acp.copernicus.org/preprints/acp-2022-749/>`_ ). The Kalman filter mode in the IMI allows users to specify a nudge factor, which is the fraction of the original emissions inventory that is retained in the prior for the next iteration. The rest of the emissions (1 - ``NudgeFactor``) come from the posterior emissions of the previous iteration.



Notes:
- Can run template setup and spinup run seperately, but DoJacobian, DoInversion, and DoPosterior must all be set to true or false at the same time
- Error should be printed if this is not the case TODO: in sanitize config file
- Preview is a prerequisite to running the inversion
- update to run for final shortened period