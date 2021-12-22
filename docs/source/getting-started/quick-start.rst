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

You will need to enter some basic personal information and a credit card number. Running the IMI Workflow is relatively inexpensive (usually on the order of 10s of dollars).
The exact cost primarily depends on the length of your simulations and how long you leave your Cloud instance running beyond the steps of the IMI Workflow.

.. note::
  Students can check out subsidized educational credits at https://aws.amazon.com/education/awseducate/.
  

2. Add S3 user permissions so you can download input data
---------------------------------------------------------

Most input data for the IMI Workflow are stored on the AWS Cloud, but are not included in your instance by default. Instead, relevant data
for your customized simulation are fetched automatically during the workflow. These automatically fetched fields include GEOS-Chem meteorology data and chemistry input fields,
as well as TROPOMI methane fields. To enable this data retrieval, you need to grant S3 download permissions to a user in your AWS account.


The easiest way to enable S3-to-EC2 downloading (and uploading) is to grant S3 access to an IAM role, which when designated to an EC2 instance (Elastic Compute Cloud, AWS's basic computing node service), will give that EC2 instance full access to S3. Instructions on how to create an IAM role with full s3 access is available at `this link to the GEOS-Chem Documentation <https://cloud-gc.readthedocs.io/en/latest/chapter03_advanced-tutorial/iam-role.html#create-a-new-iam-role>`. For more information on `IAM Roles check out the AWS Documentation <https://docs.aws.amazon.com/IAM/latest/UserGuide/id_roles.html>`.



3. Launch an instance with the IMI Workflow pre-installed
---------------------------------------------------------

Once you've setup S3 permissions on your AWS account, login to the AWS console and click on EC2.

.. figure:: img/main_console.png
  :width: 600 px

In the EC2 console, you can see your current selected region in the top right.
Choosing a region closer to your physical location will improve your network connectivity, but may result in increased costs compared to using the region where GEOS-Chem data are hosted (us-east-1 (N.Virginia)).

.. figure:: img/region_list.png
  :width: 300 px

.. _choose_ami-label:

In the EC2 console, click on "AMIs" (Amazon Machine Images) under "IMAGES" on the left navigation bar. Then select "Public images" and search for ``TODO:AMI_ID`` or ``TODO:AMI_NAME``.
This image contains the latest version of the IMI Workflow and has all the necessary software dependencies preinstalled.

.. figure:: img/search_ami.png

An AMI fully specifies the software side of your virtual system, including the operating system, software libraries, and default data files. 
Now it's time to specify the hardware for running your system. Hardware choices differ primarily in CPU and RAM counts. 

You can select from a large number of instance types at the "Step 2: Choose an Instance Type" screen. The IMI Workflow will run more quickly with a higher number of CPUs. 
TODO: choose ideal computational node. Choose the c5.9xlarge instance type, which includes 36 CPU cores and 72GB of RAM. Depending on your use case you may choose to use a different instance type that provides more/less cores and memory.

.. figure:: img/choose_instance_type.png

.. _skip-ec2-config-label:


Proceed to Step 3 and select the ``IAM Role`` you created in :ref:`Step 2 <2. Add S3 user permissions so you can download input data>`. All other config settings in step 3 can be left as the defaults.

.. figure:: img/assign_iam_to_ec2.png

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


5. Configure and run the IMI Workflow
-------------------------------------

Navigate to the IMI Workflow setup directory::

  $ cd ~/setup_CH4

Open the ``setup_ch4_inversion.sh`` script::

  $ emacs setup_ch4_inversion.sh


This script contains many settings you can modify according to your scientific needs. You should not need to modify
certain settings if you are using an instance generated from the AMI to run the workflow. These settings are italicized in the list below.

Setup settings
~~~~~~~~~~~~~~

- CreateClusterFile: Create a netCDF cluster file containing the gridboxes or regions where emissions will be perturbed. If this is set to false, the ClusterFile option must be specified further down.
- SetupTemplateRundir: Copy run directory files from GEOS-Chem and replace text in input files according to settings in setup_ch4_inversion.sh. The template run directory is needed to set up the spinup, 
  Jacobian, and posterior run directories. It only needs to be generated once.
- SetupSpinupRun: Create a run directory to spinup a new restart file representative of your model setup.
- SetupJacobianRunDirectories: Setup run directories for each of your perturbation clusters. The output from these simulations will be used to construct the Jacobian.
- SetupInversion: Copy scripts used to post-process GEOS-Chem data, build the Jacobian, and run the inversion.
- SetupPosteriorRun: Create a run directory to submit a posterior run.

Environment files
~~~~~~~~~~~~~~~~~

- *NCOEnv: Bash script to load NCO software package*
- *GCCEnv: Bash script to load software packages needed to compileand run  GEOS-Chem Classic*
- *CondaEnv: Name of conda environment containing python packages needed to execute the scripts in CH4_TROPOMI_INV. See the example in envs/Harvard-Cannon/ch4_inv.yml.*

File paths
~~~~~~~~~~
- RunName: Specify a name for your simulations.
- MyPath: Set it to the file path where you want to setup the CH4 inversion run directories.
- *DataPath: Path to non-emissions data that will replate {DATA_ROOT} token in input.geos and HEMCO_Config.rc. The default is path is ``/home/ubuntu/ExtData/``.*
  *The emissions data path is set in HEMCO_Config.rc and by default is set to ``/home/ubuntu/ExtData/HEMCO/``.*
- ClusterFile: Path to netCDF file containing clusters to perturb.
- UseBCsForRestart: Logical to determine whether boundary condition files should be used in place of a restart file.
- RestartFile: Path to the initial restart file. Sample restart files may be found in TODO.
- BCfiles: Path to boundary condition files to be used for nested grid simulations. This will replace the path in HEMCO_Config.rc.

Data download settings
~~~~~~~~~~~~~~~~~~~~~~

Your instance created from the AMI does not include much of the required input data for even a very short simulation, so it is recommended you set each of these to true
to automatically fetch required input data. Note that for each of these settings (except RestartDownload) files will not be redownloaded if they already exist at the requested path on your instance, 
so it is safe to leave these set to "true" in successive runs of the workflow. 

- SpinupDryrun: Set to "true" to automatically fetch all required GEOS-Chem meteorology and emissions input files for the spinup run from Amazon S3 (your account will incur small charges).
- ProdDryrun: Set to "true" to automatically fetch all required GEOS-Chem meteorology and emissions input files for the production runs from Amazon S3 (your account will incur small charges).
- RestartDownload: Set to "true" to automatically fetch a default restart file from Amazon S3 (your account will incur a small charge).
- BCDryrun: Set to "true" to automatically fetch default required boundary conditions files from Amazon S3 (your account will incur small charges).

Grid settings
~~~~~~~~~~~~~
- Res: Options are "4x5", "2x2.5", "0.5x0.625", and "0.25x0.3125".
- Met: Options are "merra2" or "geosfp".
- LonMin: Minimum longitude edge of your domain.
- LonMax: Maximum longitude edge of your domain. e.g. For global simulations use -180.0, 180.0. For nested NA simulations use -140.0, -40.0 (0.5x0.625); -130.0, -60.0 (0.25x0.3125).
  For nested Asia simulations use   60.0, 150.0 (0.5x0.625); 70.0, 140.0 (0.25x0.3125)
- LatMin: Minimum latitude edge of your domain.
- LatMax: Maximum latitude edge of your domain. e.g. For global simulations use -90.0, 90.0. For nested NA simulations use 10.0, 70.0 (0.5x0.625); 9.75, 60.0 (0.25x0.3125).
  For nested Asia simulations use  15.0, 55.0 (0.5x0.625); -11.0, 55.0 (0.25x0.3125).
- HalfPolar: Set to "T" for global simulations to use half polar boxes. Set to "F" for nested grid simulations.
- Levs: Set to 47 to use the reduced 47-level grid recommended for CH4 simulations.
- NestedGrid: Set to "F" for global simulations or "T" for nested simulations.
- Region: Set to "" for global or "NA", "AS", "CH", "EU" for default domains
- Buffer: Set to "0 0 0 0" for global simulations. For nested simulations, the recommendation is to use "3 3 3 3" to use 3 grid cells along the nested-grid domain for your buffer zone.

Jacobian settings
~~~~~~~~~~~~~~~~~
- PerturbValue: Perturbation value to apply to clusters in analytical inversion.

Additional settings (change options in input.geos)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- GOSAT: Set to true to use GOSAT observation operator.
- TCCON: Set to true to use TCCON observation operator.
- UseEmisSF: Set to true to use emission scale factors from a previous analytical inversion.
- UseSeparateWetlandsSF: Set to true to use separate scale factors for wetland and nonwetland emissions.
- UseOHSF: Set to true to use OH scale factors from a previous analytical inversion.
- PLANEFLIGHT: Set to true to use the planeflight diagnostic.
- HourlyCH4: Set to true to save out hourly CH4 concentrations and pressure edges.


Save and close setup_ch4_inversion.sh when you're done editing configuration settings.

TODO: Add instructions for using automated workflow

`Click here <manual-running>`__ for instructions on manually running each step of the workflow (an alternative to using the automated workflow run script).

TODO: Add instructions from README.MD to describe different customization options and using the workflow.


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

TODO: Add section (here or elsewhere) on exporting data to S3. Also add information somewhere about modifying instance type when not running IMI to save money
