Setting up Jupyter on EC2
=========================

The IMI relies on Jupyter notebooks to visualize the results of inversions. However, in order to 
view and run Jupyter notebooks you will need to set up a jupyter server on your EC2 instance and 
securely access it from your local browser. You can do this using a few different methods:

Using an SSL certificate (AWS Recommended)
------------------------------------------

Follow `these short instructions <https://docs.aws.amazon.com/dlami/latest/devguide/setup-jupyter.html>`_ to 
set up and connect to a jupyter notebook server on AWS using a self-signed SSL certificate. 
Note: If you are using Git-BASH, or similar software with ssh, you can follow the 
`Configure a Linux or macOS Client <https://docs.aws.amazon.com/dlami/latest/devguide/setup-jupyter-configure-client-linux.html>`_ 
section, which provides a simpler setup than the Windows instructions.

Using an automatically generated authentication token
-----------------------------------------------------
If the above AWS recommended method is causing trouble you can also use the following method to create and connect to 
a jupyter server using an authentication token. The authentication token is a randomly generated hash code appended to 
the jupyter server url. The token, similar to a password, verifies you have permission to access the server.

To set up a jupyter notebook server on your ec2 instance, run the following command on your **remote/EC2** terminal::

  $ jupyter notebook --no-browser --port 8080

This will start a jupyter server on port 8080 and will print out a link with an authentication token, eg::
  
  $ jupyter notebook --no-browser --port 8080
    ....
    http://localhost:8080/?token=7a7ae708966c68e631bc76ba9eae7b1d287e4747cf7072e7

Then in a new **local** terminal (or GIT-Bash) window run the following command::

  $ ssh -NL 8080:localhost:8080 -i /path/to/private_key

This creates an ssh tunnel from your ec2 instance to your local computer over port 8080, which will allow you to view your jupyter 
notebooks from your browser. Go to the link outputted from your remote serve command above (eg. http://localhost:8080/?token=7a7ae708966c68e631bc76ba9eae7b1d287e4747cf7072e7).
