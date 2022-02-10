Running the IMI with tmux
=========================

Ordinarily, when ssh'ed into an ec2 instance, any commands or processes you start will terminate once you disconnect from ssh. However, `tmux <https://man7.org/linux/man-pages/man1/tmux.1.html>`_ allows you to run programs in the background, disconnect, and reconnect later to check the progress of your inversion.

Using tmux
----------
tmux comes preinstalled on the ami used for running the IMI workflow on the cloud. To start tmux run the following:

    $ tmux 

This enters a tmux shell. From there you can run the inversion script:
    
    $ ./run_imi.sh
    
This will start the workflow. To keep it running in the background, press ``ctrl-b`` then press ``d`` (without holding ``ctrl``) to detach the tmux shell and get back to the original terminal. At this point you can disconnect the from ssh and the workflow will continue to run in the background.


To check back in on the workflow, ssh back onto the EC2 instance and run the following to attach to the tmux session running the workflow:
    
    $ tmux attach-session -t 0
