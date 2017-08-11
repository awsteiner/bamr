User's Guide
============

Installation
------------
    
The \bm executable requires the installation of `GSL
<http://www.gnu.org/software/gsl>`_ (versions 1.16 and later), `HDF5
<http://www.hdfgroup.org>`_ (versions 1.8.14 and later), the most
current version of `O2scl <http://web.utk.edu/~asteine1/o2scl>`_ from
the github repository, and MPI (tested with ``openmpi-1.8.6``). After
these four packages are successfully installed, you will need to edit
the ``makefile`` and then compile ``bamr`` before execution.

Download
--------

The most recent release version can be obtained from either
of the following::

  svn checkout svn://svn.code.sf.net/p/bamr/code/trunk bamr
  git clone https://github.com/awsteiner/bamr.git


Basic Usage
-----------
    
The basic usage is something like::

  ./bamr -model twop -run default.in -mcmc
  
to perform a one day run with model ``twop`` with the input
file in ``default.in``. 

There are several variables which can be modified with the
``set`` command, e.g.::

  ./bamr -model twop -set max_time 43200 -run default.in -mcmc

which runs for 12 hours instead of the default 24 hours. 

An example of an MPI invocation is::

  mpirun -np 4 ./bamr -set model twop -run default.in -mcmc &

which runs with four processors on the current machine.

Also try::

  ./bamr -help

which outputs some additional information on the 
relevant functions and parameters. 

Data Files
----------

The input data files are HDF5 files (typically named with a
``.o2`` extension) which contain one ``o2scl::table3d``
object giving the probability density of one neutron star
observation as a slice in that table. The command
``add-data``, which adds a neutron star data set, takes four 
arguments (and a fifth optional argument):

- The ASCII name of the neutron star for the output chains
- The input data file name
- The name of the ``o2scl::table3d`` slice in which the data 
  is stored
- The initial guess for the neutron star's mass
- The fifth, optional, argument is the name of the \ref
  ``o2scl::table3d`` object in the data file

It is assumed that the "x" grid in the data file refers to the radius
and the "y" grid refers to the mass. The data need not be normalized.
By default ``bamr`` renormalizes the data so that the maximum
probability is 1. If the parameter ``norm_max`` is set to false, then
the data is renormalized to have a total integral of 1.

Output Files
------------

Output is stored in HDF files with a prefix given by the
argument to the ``mcmc`` command, one set of files
for each processor. Output includes files with the 
following suffixes (where X is the processor index):

- ``_X_out``: Main output file containing full Markov chain(s)
  and most of the parameter values
- ``_X_scr``: Running output of entire simulation

If the executable is run directly (without ``mpirun``)
then X is always zero.

Representations of the EOS and the :math:`M-R` curve are stored on
grids for each Monte Carlo point. The number of points in the grid
is specified by the parameter ``grid_size``. The energy
density grid is specified by the limits ``e_low`` and
``e_high``, the baryon density grid is specified by the
limits ``nb_low`` and ``nb_high``, and the mass grid is
specified by the limits ``m_low`` and ``m_high``.
A default run stores pressure as a function of energy density
(in columns with prefix ``P_``), radius as a function
of mass (in columns with prefix ``R_``), central pressure
as a function of mass (in columns with prefix ``PM_``),
pressure as a function of baryon density (in columns with 
prefix ``Pnb_``), and energy per baryon as a function
of baryon density (in columns with prefix ``EoA_``).

Some Details
------------

The basic functionality is provided in the class
:cpp:class:`bamr::bamr_class` and :cpp:class:`bamr::model` . All of
the "models" (EOS parameterizations) are children of
:cpp:class:`bamr::model` .

If the initial guess has no probability, then the code will fail.
The easiest fix is just to change the initial guess.

In order to make the output more efficient, the table representing
the full Markov chain is divided up into tables with 10,000 rows
each, named ``markov_chain0``, ``markov_chain1``, and so on. The
total number of tables is stored in the HD5 output in the
integer ``n_chains``.

Different models have different optimal MC step sizes. The step size
for each parameter is chosen to be the difference betwen the high and
low limiting values divided by the value
``o2scl::mcmc_para_base::step_fac`` (which can be changed from the
command line using ``-set step_fac <value>``. Increasing or decreasing
this value may give better results.

The EOS results are stored in a table in
:cpp:member:`bamr::model_data::eos` and the TOV results are stored in
:cpp:member:`bamr::model_data::mvsr` .

Object Creation
---------------

The ``main()`` function creates a :cpp:class:`bamr::mcmc_bamr`
instance, setting the number of OpenMP threads as an argument to the
constructor. The ``mcmc_bamr`` constructor creates a
:cpp:class:`bamr::bamr_class` instance for each OpenMP thread, stored
in :cpp:member:`bamr::mcmc_bamr::bc_arr` Then, the ``mcmc_bamr``
constructor creates an instance of the :cpp:class:`bamr::settings` and
:cpp:class:`bamr::ns_data` classes for the shared pointer objects
:cpp:member:`bamr::mcmc_bamr::set` and
:cpp:member:`bamr::mcmc_bamr::nsd` respectively. These two objects are
shared between all of the :cpp:class:`bamr::bamr_class` instances The
``main()`` function then calls
:cpp:func:`bamr::mcmc_bamr::setup_cli()` and then calls the
``o2scl::cli::run_auto()`` function for the object
:cpp:member:`bamr::mcmc_bamr::cl` at which point the execution
proceeds according to the command-line options specified by the user.

Typically the user first specifies the model with ``-model``, at which
point :cpp:func:`bamr::mcmc_bamr::set_model()` creates one
descendant of :cpp:class:`bamr::model` and placed in
:cpp:member:`bamr::bamr_class::mod`. Then the user begins the Monte
Carlo with ``-mcmc`` which calls
:cpp:func:`bamr::mcmc_bamr::mcmc_func()` which in turn calls
``o2scl::mcmc_para_table::mcmc()``. Internally, this latter function
creates several (two times the number of OpenMP threads times the
number of Monte Carlo walkers) copies of :cpp:class:`bamr::model_data`
for use by the Monte Carlo.
     
Crust Model
-----------
    
The crust is computed in ``o2scl::eos_tov_interp`` using the crust EOS
from ``o2scl::eos_tov_interp::default_low_dens_eos()`` . In
O2scl, by default, the transition pressure in this function
is assumed to be the largest pressure in the crust which for the
default crust EOS is the pressure corresponding to a baryon density of
0.08 :math:`\mathrm{fm}^{-3}`. In order to ensure that the transition is
more smooth and that the core EOS (specified by the user) is used for
all pressures above a baryon density of 0.08 :math:`\mathrm{fm}^{-3}`, the
default \o2 procedure is modified. The transition pressure is
decreased by 20 percent and a pressure width of 20 percent is supplied
to ``o2scl::eos_tov_interp::set_transition()`` .
    
Model ``qstar`` from :cpp:class:`bamr::quark_star` is typically run
with ``use_crust`` set to false.

EOS Model
---------

Several EOS models are provided. New models (i.e. new
children of the :cpp:class:`bamr::model` class) must perform several tasks

- The function :cpp:func:`bamr::model::compute_eos()` should use the
  parameters to compute the EOS
- The energy density should be stored in a column named
  ``ed`` and the pressure in ``pr`` with the correct units
  set for each column (currently only ``1/fm^4`` is supported).
- If the model provides the symmetry energy and its density
  derivative, it should be stored as constants named ``"S"``
  and ``"L"`` in the table (in :math:`1/\mathrm{fm}` ).
- Causality is automatically checked in
  :cpp:func:`bamr::model::compute_star()`, but the
  :cpp:func:`bamr::model::compute_eos()` function should check that the
  pressure is not decreasing.
- Finally, it is recommended to set the interpolation type in the
  ``o2scl::table_units`` object to linear interpolation.

Post-processing Code
--------------------

The ``process`` program contains some code to analyze the ``bamr``
output files. Principally, it attempts to remove autocorrelations
from the data by separating the data into blocks, and using the
fluctuations in the mean of each block to determine the
uncertainty in the mean. For example::

  ./process -set xscale 197.33 -hist L out.o2 bamr_0_out 

creates a new file called ``out.o2`` with an object of type 
``o2scl::table`` which represents a probability distribution for $ L $
from the ``bamr`` output file named ``bamr_0_out``. The option ``-set
xscale 197.33`` ensures that the new table converts from
:math:`\mathrm{fm}^{-1}` to :math:`\mathrm{MeV}`. The five columns are
``reps``, ``avgs``, ``errs``, ``plus`` and ``minus``, which give:

- 1. The central bin values for :math:`L` (in :math:`\mathrm{MeV}`)
- 2. The probability of the specified value
- 3. The uncertainty in the probabiliy of the specified value
- 4. Column 2 plus column 3
- 5. Column 2 minus column 3

Recent Change Log
-----------------

Summer 2017: Completely reworked for OpenMP/MPI with new O2scl
class ``mcmc_para_base`` .

April 2015: Added ``process.cpp`` and created new functions \ref
bamr::model::setup_params and \ref bamr::model::remove_params() .
Added several new models.

Acknowledgements
----------------

I would like to thank Paulo Bedaque, Ed Brown, Farrukh Fattoyev, Chris
Fryer, Stefano Gandolfi, Jim Lattimer, Joonas Nattila and Will Newton
for their collaboration on these projects.

.. toctree::
   :maxdepth: 2
   :caption: Contents


  
