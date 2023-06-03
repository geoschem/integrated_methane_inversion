Using the IMI Kalman Filter mode
================================

What is a Kalman Filter Inversion?
----------------------------------
A Kalman filter is a mathematical algorithm, developed by Rudolph Kalman, that estimates the state of a system by combining measurements and predictions while considering uncertainties. It operates recursively, continuously updating its estimate of the system state based on new measurements.

Kalman filters can be applied in atmospheric inversions by dividing an inversion period is into smaller time intervals, such as weekly chunks. An inversion is sequentially run for each interval, estimating the emissions for that specific period based on measurements and predictions. The resulting optimized emissions are then used as prior emissions for the next interval, allowing the prior emissions of each successive week to be informed by the previous weeks.

Why use Kalman Filter Mode?
---------------------------
This approach enables tracking of how emissions change over time and provides insights into their distribution throughout the inversion period. By using the Kalman filter mode in the inversion, users can calculate intermediate emissions at the desired update frequency, such as weekly, revealing the temporal evolution of emissions.

How to use the Kalman Filter mode
---------------------------------
Example Kalman filter config variables:

::

    ## Kalman filter options
    KalmanMode: true
    UpdateFreqDays: 7
    NudgeFactor: 0.1
    FirstPeriod: 1
      

Notes:
- Can run template setup and spinup run seperately, but DoJacobian, DoInversion, and DoPosterior must all be set to true or false at the same time
- Error should be printed if this is not the case TODO: in sanitize config file
- Preview is a prerequisite to running the inversion
- update to run for final shortened period