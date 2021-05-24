AMI Specifications
==================

TODO: Move this to its own page on the RTD site.


The latest Amazon Machine Image (AMI) for the UMI Workflow contains the following software libraries:

- GNU Compiler Collection 8.2.0
- NetCDF-Fortran 4.5.3
- Slurm 20.11.0.1
- Python 3.8.7
- GEOS-Chem Classic 13.0.2 TODO: Update to 13.0.2



.. _quick-start-label:

Quick start guide for new users
===============================


1. Sign up for an Amazon Web Services (AWS) account
---------------------------------------------------

If you do not already have an AWS account (whether personal or through your institution), you'll need to sign up for one.
Go to http://aws.amazon.com and click on "Create an AWS Account" in the upper-right corner:

.. figure:: img/create_aws_account.png
  :target: https://aws.amazon.com
  :width: 400 px

You will need to enter some basic personal information and a credit card number. Running the UMI Workflow is relatively inexpensive (usually on the order of 10s of dollars).
The exact cost primarily depends on the length of your simulations and how long you leave your Cloud instance running beyond the steps of the UMI Workflow.

.. note::
  Students can check out subsidized educational credits at https://aws.amazon.com/education/awseducate/.
  

2. Add S3 user permissions so you can download input data
---------------------------------------------------------

Most input data for the UMI Workflow are stored on the AWS Cloud, but are not included in your instance by default. Instead, relevant data
for your customized simulation are fetched automatically during the workflow. These automatically fetched fields include data GEOS-Chem meteorology and chemistry input fields,
as well as TROPOMI methane fields. To enable this data retrieval, you need to grant S3 download permissions to a user in your AWS account.


The easiest way to enable S3-to-EC2 downloading (and uploading) is to grant S3 access to all EC2 (Elastic Compute Cloud, AWS's basic computing node service) 
instances that are launched on your account.
TODO: Transfer these instructions to this page. Instructions on how to do so are available at 
https://cloud-gc.readthedocs.io/en/latest/chapter03_advanced-tutorial/iam-role.html#create-a-new-iam-role.



3. Launch an instance with the UMI Workflow pre-installed
---------------------------------------------------------

Once you've setup S3 permissions on your AWS account, login to the AWS console and click on EC2.

.. figure:: img/main_console.png
  :width: 600 px

In the EC2 console, you can see your current selected region in the top right.
Choosing a region closer to your physical location will improve your network connectivity, but may result in increased costs compared to using the region
where GEOS-Chem data are hosted (US East (N.Virginia)).

.. figure:: img/region_list.png
  :width: 300 px

.. _choose_ami-label:

In the EC2 console, click on "AMIs" (Amazon Machine Images) under "IMAGES" on the left navigation bar. Then select "Public images" and search for ``TODO:AMI_ID`` or ``TODO:AMI_NAME``.
This image contains the latest version of the UMI Workflow.

.. figure:: img/search_ami.png

An AMI fully specifies the software side of your virtual system, including the operating system, software libraries, and default data files. 
Now it's time to specify the hardware for running your system. Hardware choices differ primarily in CPU and RAM counts. 

You can select from a large number of instance types at the "Step 2: Choose an Instance Type" screen. The UMI Workflow will run more quickly with a higher number of CPUs. 
TODO: choose ideal computational node (this one may be unnecessarily powerful as it is built for inter-node connection). Choose the c5n.9xlarge instance type, which includes 36 CPU cores and 96GB of RAM. 

.. figure:: img/choose_instance_type.png

.. _skip-ec2-config-label:

**Then, just click on "Review and Launch".** You don't need to touch other options this time. This brings you to "Step 7: Review Instance Launch". Click on the Launch button again.

.. _keypair-label:

When you first use EC2, you will be asked to create and download a file called a "Key Pair". It is equivalent to the password you enter to ``ssh`` to your local server.

Give your "Key Pair" a name, click on "Download Key Pair", and finally click on "Launch Instances". In the future, you can simply select "Choose an existing Key Pair", select your previously created Key Pair, and launch.

.. figure:: img/key_pair.png
  :width: 500 px


Once launched, you can monitor the instance in the EC2-Instance console as shown below. Within < 1min of initialization, "Instance State" should become "running" (refresh the page if the status stays as "pending"):

.. figure:: img/running_instance.png

You now have your own system running on the cloud! Note that you will be charged every hour that you leave this instance running, so make sure to do the 
:ref:`final tutorial step: shutdown the server <shutdown-label>` if you need to pause your work to avoid being charged continuously.

.. _login_ec2-label:

4. Login to your instance
------------------------------

Select your instance, click on the "Connect" button (shown in the above figure) near the blue "Launch Instance" button, then you should see this instruction page:

.. figure:: img/connect_instruction.png
  :width: 500 px

- On Mac or Linux, use the ``ssh -i ...`` command under "Example" to connect to the server in the terminal. Some minor changes are needed:

  (1) ``cd`` to the directory where your Key Pair is stored (people often put the key in ``~/.ssh/`` but any directory is fine.)
  (2) Use ``chmod 400 your-key-name.pem`` to change the key pair's permission (also mentioned in the above figure; only need to do this the first time you login).
  (3) Change the user name in that command from ``root`` to ``ubuntu``, so the full command will be like ``ssh -i "your-key-name.pem" ubuntu@xxx.amazonaws.com``

- On Windows, you can install `Git-BASH <https://gitforwindows.org>`_ to emulate a Linux terminal. Simply accept all default options during installation, as the goal here is just to use Bash, not Git. 
Alternatively, you can use `MobaXterm <http://angus.readthedocs.io/en/2016/amazon/log-in-with-mobaxterm-win.html>`_, `Putty <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/putty.html>`_, 
`Linux Subsystem <https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/WSL.html>`_ or `PowerShell with OpenSSH <https://blogs.msdn.microsoft.com/powershell/2017/12/15/using-the-openssh-beta-in-windows-10-fall-creators-update-and-windows-server-1709/>`_. 
The Git-BASH solution should be the most painless, but these other options should work as well.


Once you've followed the above instructions, you should see a "Welcome to Ubuntu" message indicating you've logged into your new EC2 instance.


5. Configure and run the UMI Workflow
-------------------------------------

Navigate to the UMI Workflow setup directory::

  $ cd ~/setup_CH4

Open the ``setup_ch4_inversion.sh`` script::

  $ emacs setup_ch4_inversion.sh
  
TODO: Add instructions from README.MD to describe different customization options and using the workflow.

There are some AWS-specific configuration options that you should leave alone. ``DATA_PATH`` points to an existing folder on your instance containing a small quantity of GEOS-Chem input data.
More data will be downloaded to this folder. You can point ``BC_FILES`` to any folder where you would like to download boundary condition files, but you should leave ``BC_DRYRUN`` set to ``true``
if you plan on using default boundary conditions files so that this script automatically downloads missing BC files.
``SPINUP_DRYRUN`` and ``PROD_DRYRUN`` determine whether the workflow scripts will attempt to automatically download other required GEOS-Chem data. 
You can always leave these options set to ``true`` because the workflow scripts will only download required data that is not already located in ``DATA_PATH``. 

TODO: More workflow instructions here


The ``FETCHTROPOMI`` variable in ``run_inversion.sh`` determines whether TROPOMI data files will be automatically downloaded for the designated inversion timeframe.
Note that unlike the GEOS-Chem data download, **TROPOMI data files will be downloaded even if they already exist on your instance**. Set ``FETCHTROPOMI`` to ``false`` 
if you're rerunning for a certain time period to avoid unnecessary file downloads. 

TODO: More workflow instructions here


6. Analyze output data with Python
----------------------------------

TODO: Fill in this section once all scripts are ready


.. _shutdown-label:

7. Shut down the instance
-------------------------

If you're done using your instance for awhile or don't plan on using it again, you should either shutdown or terminate your instance. 
Shutting down or terminating your instance will minimize or completely stop, respectively, new charges to your account.


Right-click on the instance in your console to get this menu:

.. image:: img/terminate.png

You have two options now, "Stop" to shutdown or "Terminate" to completely delete your instance:

- "Stop" will make the system inactive. You won't be charged for CPU time, but you will be charged a negligible disk storage fee.
  You can restart the server at any time and all files will be preserved. When an instance is stopped, you can also change its hardware type (right click on the instance -> "Instance Settings" -> "Change Instance Type") 
- "Terminate" will completely remove that instance so you won't be charged for it any further.
  Unless you save your system as an AMI or transfer the data to other storage services, you will lose all your data and software.

TODO: Add section (here or elsewhere) on exporting data to S3. Also add information somewhere about modifying instance type when not running UMI to save money
