.. _manual-running-label:

Manually running the workflow
=============================

Make sure to edit any commands containing ``${VARIABLE_NAME}`` based on how you configured the workflow.

1. From `/home/ubuntu/setup_CH4/`, type ``./setup_ch4_inversion.sh``.
2. Navigate to ``${MY_PATH}/${RUN_NAME}`` as as specified in `setup_ch4_inversion.sh`.
3. If you set DO_SPINUP to true, submit the spinup simulation by navigating to the SpinupRun directory and typing ``sbatch ${RUN_NAME}_Spinup.run``.
4. Submit the Jacobian runs by navigating to the `jacobian_runs` directory.
 4a. If running just the base simulation (run 0000, e.g. for prior or
         posterior), you can submit your run using the `${RUN_NAME}_0000.run` script
         in your `run_dirs/${RUN_NAME}_0000` run directory. You will first need to
         uncomment the ``#SBATCH`` lines in that file and make any necessary changes
         to those options in order to submit to the queue. Then, at the command
         line, type ``sbatch ${RUN_NAME}_0000.run`` at the command line.
 4b. If running the full set of runs to construct the Jacobian, type ``./submit_jacobian_simulations_array.sh``.

5. Run the inversion by navigating to `inversion/` and typing ``./run_inversion.sh``. The ``FETCHTROPOMI`` variable in ``run_inversion.sh`` determines whether TROPOMI data files will be automatically 
downloaded for the designated inversion timeframe. Note that unlike the GEOS-Chem data download, **TROPOMI data files will be downloaded even if they already exist on your instance**. Set ``FETCHTROPOMI`` to ``false`` 
if you're rerunning for a certain time period to avoid unnecessary file downloads. 
