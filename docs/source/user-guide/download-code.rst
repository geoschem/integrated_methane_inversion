.. _get-code-steps:

#####################
Download instructions
#####################

Follow these directions to download the IMI source code. If you are running on AWS, the AMI already includes the IMI source code so you can skip these steps.

=============
Clone the IMI
=============

Before cloning or running the IMI, make sure that you meet the `minimum software requirements <https://geos-chem.readthedocs.io/en/stable/getting-started/system-req-soft.html>`_. 
The IMI uses Git for downloading and maintaining versioning of its code base.

To download the latest stable IMI version, type:

.. code-block:: console

    $ git clone https://github.com/geoschem/integrated_methane_inversion.git

This will clone the IMI code into a local folder named ``integrated methane_inversion``.

Next, download and initialize the submodules for the IMI with:

.. code-block:: console

    $ git submodule update --init --recursive

.. tip:: 
   To download IMI source code into a folder named something other than :file:`integrated_methane_inversion`, supply the name of the
   folder at the end of the :command:`git clone` command.  For example:
   
   .. code-block:: console

      $ git clone https://github.com/geoschem/integrated_methane_inversion.git my-imi-dir

Navigate to the IMI folder and view the contents:

.. code-block:: console

    $ cd integrated_methane_inversion
    $ ls
    config.yml  docs/  envs/  LICENSE.md  README.md  resources/  run_imi.sh*  src/

.. tip::

   To use an older IMI version (e.g. 2.0.0), follow
   these additional steps:

   .. code-block:: console

      $ git checkout tags/2.0.0                  # Points HEAD to the tag "2.0.0"
      $ git branch version_2.0.0                 # Creates a new branch at tag "2.0.0"
      $ git checkout version_2.0.0               # Checks out the version_2.0.0 branch

   You can do this for any tag in the version history.   For a list of
   all tags, type:

   .. code-block:: console

      $ git tag

   If you have any unsaved  changes, make sure you commit those to a
   branch prior to updating versions.

.. _get-code-steps-branch:

====================================
Create a branch in IMI for your work
====================================

When you first clone the IMI code, it will will be in a **detached HEAD state**. In
other words, the code is checked out but a branch is not
created. Adding new code to a detached HEAD state is very dangerous
and should be avoided. You should instead make a branch it the same
point as the detached HEAD, and then add your own modifications into
that branch.

Navigate to your IMI code repository and type:

.. code-block:: console

    $ git branch

You will see output similar to this:

.. code-block:: text

    *(HEAD detached at xxxxxxxx)
    main

where ``xxxxxxxx`` denotes the hash of the commit at which the code
has been checked out.

At ths point, you may now create a branch in which to store your own
modifications to the GEOS-Chem science codebase.  Type:

.. code-block:: console

   $ git branch feature/my-updates
   $ git checkout feature/my-updates

Instead of :file:`feature/my-updates`, you may choose a name that reflects
the nature of your updates.

.. note::

   This naming convention adheres to the `Github Flow
   <https://guides.github.com/introduction/flow/>`_
   conventions (i.e. new feature branches start with
   :file:`feature/`, bug fix branches start with :file:`bugfix/`, etc.

If you now type:

.. code-block:: console

   $ git branch

You will see that we are checked out onto the branch that you just
created and are no longer in detached HEAD state.

.. code-block:: text

   * feature/my-updates
   main

At this point, you may proceed to add your modifications into the IMI and follow the :doc:`guidelines for contributing <../help-and-reference/CONTRIBUTING>` any mature updates back into the IMI.

.. _get-code-steps-info:

=============
Forward model
=============

The IMI uses `GEOS-Chem <https://geoschem.github.io/>`_ as the forward model to
perturb emissions during the jacobian run phase of the IMI. The first time you execute the IMI via the
:file:`run_imi.sh` script, it will download (via ``git clone``) the expected version of GEOS-Chem as
``GCClassic/`` or ``GCHP/`` (the high-performance version of GEOS-Chem) within the IMI code directory so
that you have:

.. code-block:: console

    $ cd integrated_methane_inversion
    $ ls
    config.yml  envs/       LICENSE.md  resources/   src/
    docs/       GCClassic/  README.md   run_imi.sh*

See the model documentation for more information:

* `GCClassic documentation <https://geos-chem.readthedocs.io/en/stable>`_
* `GCHP documentation <https://gchp.readthedocs.io/en/stable/>`_

