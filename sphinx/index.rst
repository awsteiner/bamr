bamr documentation
==================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
    
This User's Guide describes the open-source MPI implementation of
a Bayesian analysis of mass and radius data to determine the mass
versus radius curve and the equation of state of dense matter.
This package will principally be useful for those physicists and
astrophysicists who are already familiar with C++ and are
interested in modifying this code for their own use. This 
documentation was generated from git commit

.. include:: static/rev.txt
   :literal:

This code was originally supported by Chandra grant
`TM1-12003X
<http://cxc.harvard.edu/target_lists/cycle12/theory12.html#12400566>`_.

This is a beta version. Use at your own risk.

Currently, bamr is hosted as a `git repository
<http://www.github.com/awsteiner/bamr>`_ . It is also listed at the
`Astrophysics Source Code Library (ASCL) <http://ascl.net/1408.020>`_
and at the `ADS <http://adsabs.harvard.edu/abs/2014ascl.soft08020S>`_
. This HTML documentation is hosted `here
<http://web.utk.edu/~asteine1/bamr>`_ . While you are not required to
do so, please consider citing the ASCL entry following the method
described `here <http://ascl.net/wordpress/?page_id=351>`_ , and the
relevant references in the bibliography below.

If you are considering using this code for your research, I
encourage you to contact me so that I can help you with the
details and so that you can let me know if and how this code is
useful to you. Nevertheless, you are not required to contact me
and I will be improving documentation and updating this code as
time permits.

.. note:: The bamr source code is undergoing some major 
	  upgrades which have not yet been documented. 

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
  
to perform a one day run with model \c twop with the input
file in \c default.in. 

There are several variables which can be modified with the
\c set command, e.g.::

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
``.o2`` extension) which contain one \ref o2scl::table3d
object giving the probability density of one neutron star
observation as a slice in that table. The command
``add-data``, which adds a neutron star data set, takes four 
arguments (and a fifth optional argument):
- The ASCII name of the neutron star for the output chains
- The input data file name
- The name of the \ref o2scl::table3d slice in which the data 
is stored
- The initial guess for the neutron star's mass
- The fifth, optional, argument is the name of the \ref
o2scl::table3d object in the data file

It is assumed that the "x" grid in the data file refers to the
radius and the "y" grid refers to the mass. The data need not
be normalized. By default \c bm renormalizes the data so that
the maximum probability is 1. If the parameter ``norm_max``
is set to false, then the data is renormalized to have a total
integral of 1.

Output Files
------------

Output is stored in HDF files with a prefix given by the
argument to the \c mcmc command, one set of files
for each processor. Output includes files with the 
following suffixes (where X is the processor index):
- \c _X_out: Main output file containing full Markov chain(s)
and most of the parameter values
- \c _X_scr: Running output of entire simulation

If the executable is run directly (without ``mpirun``)
then X is always zero.

Representations of the EOS and the \f$ M-R \f$ curve are stored on
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

\hline
\section detail_sect Some Details

The basic functionality is provided in the class \ref
bamr::bamr_class and \ref bamr::model . All of the "models" (EOS
parameterizations) are children of \ref bamr::model .

If the initial guess has no probability, then the code will fail.
The easiest fix is just to change the initial guess.

In order to make the output more efficient, the table representing
the full Markov chain is divided up into tables with 10,000 rows
each, named \c markov_chain0, \c markov_chain1, and so on. The
total number of tables is stored in the HD5 output in the
integer ``n_chains``.

Different models have different optimal MC step sizes. The step
size for each parameter is chosen to be the difference betwen the
high and low limiting values divided by the value \c
o2scl::mcmc_base::step_fac (which can be changed from the command
line using "-set step_fac <value>". Increasing or decreasing this
value may give better results.

The EOS results are stored in a table in \ref
bamr::model_data::eos and the TOV results are stored in \ref
bamr::model_data::mvsr .

\hline
\section crust_sect Crust Model
    
The crust is computed in \ref o2scl::eos_tov_interp using the
crust EOS from \ref o2scl::eos_tov_interp::default_low_dens_eos()
. In \o2, by default, the transition pressure in this function is
assumed to be the largest pressure in the crust which for the
default crust EOS is the pressure corresponding to a baryon
density of 0.08 $ \mathrm{fm}^{-3} $. In order to ensure that
the transition is more smooth and that the core EOS (specified by
the user) is used for all pressures above a baryon density of 0.08
$ \mathrm{fm}^{-3} $, the default \o2 procedure is modified.
The transition pressure is decreased by 20 percent and a pressure
width of 20 percent is supplied to
o2scl::eos_tov_interp::set_transition() .
    
Model ``qstar`` from \ref bamr::quark_star is typically run
with ``use_crust`` set to false.

EOS Model
---------

Several EOS models are provided. New models (i.e. new
children of the \ref bamr::model class) must perform several tasks
- The function \ref bamr::model::compute_eos() should use the
parameters to compute the EOS
- The energy density should be stored in a column named
``ed`` and the pressure in ``pr`` with the correct units
set for each column (currently only ``1/fm^4`` is supported).
- If the model provides the symmetry energy and its density
derivative, it should be stored as constants named ``"S"``
and ``"L"`` in the table (in $ 1/\mathrm{fm} $ ).
- Causality is automatically checked in bamr::compute_star(), but
the \ref bamr::model::compute_eos() function should check that the
pressure is not decreasing.
- Finally, it is recommended to set the interpolation type in the
\ref o2scl::table_units object to linear interpolation.

Post-processing Code
--------------------

The ``process`` program contains some code to analyze the \c bamr
output files. Principally, it attempts to remove autocorrelations
from the data by separating the data into blocks, and using the
fluctuations in the mean of each block to determine the
uncertainty in the mean. For example, 

\code
./process -set xscale 197.33 -hist L out.o2 bamr_0_out 
\endcode

creates a new file called \c out.o2 with an object of type \ref
o2scl::table which represents a probability distribution for $ L
$ from the \c bamr output file named \c bamr_0_out . The option
``-set xscale 197.33`` ensures that the new table converts
from $ \mathrm{fm}^{-1} $ to $ \mathrm{Mev} $ . The five
columns are \c reps, \c avgs, \c errs, \c plus and \c minus, which
give:

- The central bin values for $ L $ (in $ \mathrm{MeV} $)
- The probability of the specified value
- The uncertainty in the probabiliy of the specified value
- Column 2 plus column 3
- Column 2 minus column 3

\hline
\section changelog_sect Recent Change Log

April 2015: Added process.cpp and created new functions \ref
bamr::model::setup_params and \ref bamr::model::remove_params() .
Added several new models.

\hline
\section ack_sect Acknowledgements

I would like to thank Paulo Bedaque, Ed Brown, Farrukh Fattoyev,
Chris Fryer, Stefano Gandolfi, Jim Lattimer, Joonas
N&auml;ttil&auml; and Will Newton for their collaboration on these
projects.

\hline
\section ref_sect Bibliography

Some of the references which contain links should direct you to
the work referred to using its DOI identifer.

Class Documentation
===================
	     
bamr_class
----------
	     
.. doxygenclass:: bamr::bamr_class
   :members:
   :protected-members:
   :undoc-members:

model_data
----------
      
.. doxygenclass:: bamr::model_data
   :members:
   :protected-members:
   :undoc-members:

settings
--------
      
.. doxygenclass:: bamr::settings
   :members:
   :protected-members:
   :undoc-members:

ns_data
-------
      
.. doxygenclass:: bamr::ns_data
   :members:
   :protected-members:
   :undoc-members:

model
-----
      
.. doxygenclass:: bamr::model
   :members:
   :protected-members:
   :undoc-members:

two_polytropes
--------------
      
.. doxygenclass:: bamr::two_polytropes
   :members:
   :protected-members:
   :undoc-members:

alt_polytropes
--------------
      
.. doxygenclass:: bamr::alt_polytropes
   :members:
   :protected-members:
   :undoc-members:

fixed_pressure
--------------
      
.. doxygenclass:: bamr::fixed_pressure
   :members:
   :protected-members:
   :undoc-members:

generic_quarks
--------------
      
.. doxygenclass:: bamr::generic_quarks
   :members:
   :protected-members:
   :undoc-members:

quark_star
----------
      
.. doxygenclass:: bamr::quark_star
   :members:
   :protected-members:
   :undoc-members:

qmc_neut
--------
      
.. doxygenclass:: bamr::qmc_neut
   :members:
   :protected-members:
   :undoc-members:

qmc_fixp
--------
      
.. doxygenclass:: bamr::qmc_fixp
   :members:
   :protected-members:
   :undoc-members:

qmc_threep
----------
      
.. doxygenclass:: bamr::qmc_threep
   :members:
   :protected-members:
   :undoc-members:

qmc_twolines
------------
      
.. doxygenclass:: bamr::qmc_twolines
   :members:
   :protected-members:
   :undoc-members:

nstar_cold2
-----------
      
.. doxygenclass:: bamr::nstar_cold2
   :members:
   :protected-members:
   :undoc-members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
