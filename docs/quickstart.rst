Quick Start
*********************

Minimal details to start using the program.

Running The Simulator
=====================

To run the program from the command line simply run the *simulator.py* file in the project top level directory, e.g. ``python PATH_TO_PROJECT\simulator.py``.  Running with the command line flag ``-h`` will show the arguments that can be given.  These are:

* ``-h`` Show help message and exit
* ``-c`` Name of configuration file to use.  Default is *config.ini*
* ``-k`` Generate keyfile and exit if this flag is present
* ``-d`` Generate default configuration file and exit if this flag is present.  Also optionally takes an argument of the file name to use.

Configuration File
==================

All parameters for an epidemic simulation are stored within a configuration file.  This is a ``.ini`` format document giving *key = value* pairs for any required parameter keys on separate lines.

The ``-d`` flag above generates a default configuration file that can be used to run the simulator.  This only specifies the required keys, and uses their default values.

The ``-k`` flag above generates a keyfile.  The keyfile gives details on all possible configuration options.  It gives the name of the key, default value, whether it is required or not, and a comment explaining what the key does.

For more details on the configuration options, see :ref:`Configuration Options <Configuration-Options>` and the keyfile.

Host Position File
==================

One of the required keys in the configuration file is *HostPosFile*.  This gives the name of the file containing the host locations.  The style of this is the same as the *Webidemics* input.  The first line is the total number of hosts, with each subsequent line the x and y positions separated with a space.  For example::

  5
  0 0
  1 0
  0 1
  1 1
  0.5 0.5

gives the following host structure::

  x x
   x
  x x

The default name for this file is *hosts.txt*.  More details see :ref:`Host Files <Host-Files>`.

Host Initial Conditions File
=============================

The final required input file specifies the host initial conditions.  Similarly to the host position file, the first line gives the number of hosts.

.. note::
  The number of hosts should match the number of hosts in the host position file.

Subsequent lines then give the initial state of the corresponding host in the host position file.  These states should be single characters matching those present in the model i.e. S/E/C/D/I/R.

More details see :ref:`Host Files <Host-Files>`.

Running within Python
===========================

As well as from the command line as shown above, the simulator can also be used from within a python environment.  The most basic usage is through the main function of the module::

  import IndividualSimulator

  all_results = IndividualSimulator.main(config_file="config.ini")

This is the same as running ``python simulator.py -c config.ini`` as shown above.  Similar to the ``-d`` and ``-k`` options on the command line, the main function has optional arguments ``keyFile`` and ``defaultConfig``, which are by default ``False`` adn ``None`` respectively.
