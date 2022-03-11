Running the IMI with tmux
=========================

The IMI can be run with `tmux <https://man7.org/linux/man-pages/man1/tmux.1.html>`_ as an alternative to sbatch. Like sbatch, tmux
allows you to run a program on your EC2 instance, disconnect, and then reconnect later to check progress. 

Because of the way the IMI is parallelized, using tmux can grant a small to moderate speed-up.

.. note::
    Before running the IMI with tmux, make sure the ``UseSlurm`` option in the :doc:`configuration file <../getting-started/imi-config-file>` 
    is set to ``false``.

Using tmux
----------
tmux comes preinstalled on the AMI. To start tmux run the following:

    $ tmux 

This enters a tmux shell. From there you can run the inversion script:
    
    $ ./run_imi.sh > imi_output.log
    
This will start the workflow. To keep it running in the background, press ``ctrl-b``. 
Then press ``d`` (without holding ``ctrl``) to detach the tmux shell and get back to the original terminal.
At this point you can disconnect the from ssh and the IMI will continue to run in the background.

To check back in on the IMI, ssh back onto the EC2 instance and run the following to attach to the active tmux session:
    
    $ tmux attach-session -t 0
