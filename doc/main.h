/** \file main.h
    \brief File containing user's guide documentation
*/
/** \mainpage 

    \section ug_section User's Guide
    
    This User's Guide describes the open-source MPI implementation of
    a Bayesian analysis of mass and radius data to determine the mass
    versus radius curve and the equation of state of dense matter.
    This package will principally be useful for those physicists and
    astrophysicists who are already familiar with C++ and are
    interested in modifying this code for their own use. This 
    documentation was generated from git commit
    \htmlinclude rev.txt

    This code was originally supported by Chandra grant
    <a href="http://cxc.harvard.edu/target_lists/cycle12/theory12.html#12400566">TM1-12003X</a>.

    This is a beta version. Use at your own risk.

    Currently, \bm is hosted as a git repository at
    http://www.github.com/awsteiner/bamr . The corresponding entry at
    the Astrophysics Source Code Library (ASCL) is
    http://ascl.net/1408.020 and at the ADS is
    http://adsabs.harvard.edu/abs/2014ascl.soft08020S . This HTML
    documentation is hosted at http://web.utk.edu/~asteine1/bamr .
    While you are not required to do so, please consider citing the
    ASCL entry following the method described at
    http://ascl.net/wordpress/?page_id=351 , and the relevant
    references in the bibliography below.

    If you are considering using this code for your research, I
    encourage you to contact me so that I can help you with the
    details and so that you can let me know if and how this code is
    useful to you. Nevertheless, you are not required to contact me
    and I will be improving documentation and updating this code as
    time permits.
    
    \hline
    \section contents_sect Contents

    - \ref install_sect
    - \ref dl_sect
    - \ref usage_sect
    - \ref infile_sect 
    - \ref outfile_sect
    - \ref detail_sect
    - \ref crust_sect
    - \ref model_sect
    - \ref func_stack_sect
    - \ref postproc_sect
    - \ref changelog_sect
    - \ref ack_sect
    - \ref ref_sect
    - \ref license_page

    \hline
    \section install_sect Installation
    
    The \bm executable requires the installation of 
    <a href="http://www.gnu.org/software/gsl">GSL</a> 
    (versions 1.16 and later),     
    <a href="http://www.hdfgroup.org">HDF5</a> 
    (versions 1.8.14 and later), 
    \htmlonly
    the most current version of
    <a href="http://web.utk.edu/~asteine1/o2scl">O<span style='position:
    relative; top: 0.3em; font-size: 0.8em'>2</span>scl</a>
    \endhtmlonly
    from the github repository, and MPI (tested with openmpi-1.8.6).
    After these four packages are successfully installed, you will
    need to edit the \c makefile and then compile \bm before
    execution.

    \hline
    \section dl_sect Download

    The most recent release version can be obtained from either
    of the following:
    \verbatim
    svn checkout svn://svn.code.sf.net/p/bamr/code/trunk bamr
    git clone https://github.com/awsteiner/bamr.git
    \endverbatim

    \hline
    \section usage_sect Basic Usage
    
    The basic usage is something like
    \verbatim
    ./bamr -model twop -run default.in -mcmc
    \endverbatim
    to perform a one day run with model \c twop with the input
    file in \c default.in. 

    There are several variables which can be modified with the
    \c set command, e.g. 
    \verbatim
    ./bamr -model twop -set max_time 43200 -run default.in -mcmc
    \endverbatim
    which runs for 12 hours instead of the default 24 hours. 

    An example of an MPI invocation is
    \verbatim
    mpirun -np 4 ./bamr -set model twop -run default.in -mcmc &
    \endverbatim
    which runs with four processors on the current machine.

    Also try
    \verbatim
    ./bamr -help
    \endverbatim
    which outputs some additional information on the 
    relevant functions and parameters. 

    \hline
    \section infile_sect Data Files

    The input data files are HDF5 files (typically named with a
    <tt>.o2</tt> extension) which contain one \ref o2scl::table3d
    object giving the probability density of one neutron star
    observation as a slice in that table. The command
    <tt>add-data</tt>, which adds a neutron star data set, takes four 
    arguments (and a fifth optional argument):
    - The ASCII name of the neutron star for the output chains
    - The input data file name
    - The name of the \ref o2scl::table3d slice in which the data 
    is stored
    - The initial guess for the neutron star's mass
    - The fifth argument is the name of the \ref o2scl::table3d object
    in the data file
    
    It is assumed that the "x" grid in the data file refers to the
    radius and the "y" grid refers to the mass. The data need not
    be normalized. By default \c bm renormalizes the data so that
    the maximum probability is 1. If the parameter <tt>norm_max</tt>
    is set to false, then the data is renormalized to have a total
    integral of 1.

    \hline
    \section outfile_sect Output Files

    Output is stored in HDF files with a prefix given by the
    argument to the \c mcmc command, one set of files
    for each processor. Output includes files with the 
    following suffixes (where X is the processor index):
    - \c _X_out: Main output file containing full Markov chain(s)
    and most of the parameter values
    - \c _X_scr: Running output of entire simulation

    If the executable is run directly (without <tt>mpirun</tt>)
    then X is always zero.

    Representations of the EOS and the \f$ M-R \f$ curve are stored on
    grids for each Monte Carlo point. The number of points in the grid
    is specified by the parameter <tt>grid_size</tt>. The energy
    density grid is specified by the limits <tt>e_low</tt> and
    <tt>e_high</tt>, the baryon density grid is specified by the
    limits <tt>nb_low</tt> and <tt>nb_high</tt>, and the mass grid is
    specified by the limits <tt>m_low</tt> and <tt>m_high</tt>.
    A default run stores pressure as a function of energy density
    (in columns with prefix <tt>P_</tt>), radius as a function
    of mass (in columns with prefix <tt>R_</tt>), central pressure
    as a function of mass (in columns with prefix <tt>PM_</tt>),
    pressure as a function of baryon density (in columns with 
    prefix <tt>Pnb_</tt>), and energy per baryon as a function
    of baryon density (in columns with prefix <tt>EoA_</tt>).

    \hline
    \section detail_sect Some Details

    The basic functionality is provided in the \ref bamr::bamr_class
    and each Monte Carlo point is an object of type \ref bamr::entry.
    All of the "models" (EOS parameterizations) are children of \ref
    bamr::model class.

    If the initial guess has no probability, then the code will fail.
    This is indicated by the line \c "Initial weight zero." in
    the \c _scr file. The easiest fix is just to change the initial 
    guess.

    In order to make the output more efficient, the table representing
    the full Markov chain is divided up into tables with 10,000 rows
    each, named \c markov_chain0, \c markov_chain1, and so on. The
    total number of tables is stored in the integer <tt>n_chains</tt>.

    Different models have different optimal MC step sizes. The step
    size for each parameter is chosen to be the difference betwen the
    high and low limiting values divided by the value \c step_fac .
    Increasing or decreasing this value may give better results.

    The EOS results are stored in a table in the \ref bamr::model::cns
    data member and the TOV results are stored in a table in the \ref
    o2scl::tov_solve object. In order to prevent needless copying of
    the EOS and TOV tables back and forth, the code keeps two
    instances of the model objects (\ref bamr::bamr_class::modp and
    \ref bamr::bamr_class::modp2) and two copies of the \ref
    o2scl::tov_solve objects (\ref bamr::bamr_class::ts and \ref
    bamr::bamr_class::ts2). The class \ref bamr::bamr_class flips back
    and forth between these two result sets depending on whether or
    not the MCMC algorithm gives an acceptance or a rejection.

    In order to ensure that several processors do not try to
    access the same input file at the same time, MPI
    calls are used to force each processor to take turns.

    \hline
    \section crust_sect Crust Model
    
    The crust is computed in \ref o2scl::eos_tov_interp using the
    crust EOS from \ref o2scl::eos_tov_interp::default_low_dens_eos()
    . In \o2, by default, the transition pressure in this function is
    assumed to be the largest pressure in the crust which for the
    default crust EOS is the pressure corresponding to a baryon
    density of 0.08 \f$ \mathrm{fm}^{-3} \f$. In order to ensure that
    the transition is more smooth and that the core EOS (specified by
    the user) is used for all pressures above a baryon density of 0.08
    \f$ \mathrm{fm}^{-3} \f$, the default \o2 procedure is modified.
    The transition pressure is decreased by 20 percent and a pressure
    width of 20 percent is supplied to
    o2scl::eos_tov_interp::set_transition() .
    
    Model <tt>qstar</tt> from \ref bamr::quark_star is typically run
    with <tt>use_crust</tt> set to false.

    \hline
    \section model_sect EOS Model

    Several EOS models are provided. New models (i.e. new
    children of the \ref bamr::model class) must perform several tasks

    - The function \ref bamr::model::compute_eos() should use the
    parameters in the \ref bamr::entry argument to compute the EOS and
    store it in the object returned by \ref
    o2scl::nstar_cold::get_eos_results().

    - The energy density should be stored in a column named
    <tt>ed</tt> and the pressure in <tt>pr</tt> with the correct units
    set for each column (currently only <tt>1/fm^4</tt> is supported).

    - If \ref bamr::bamr_class::baryon_density is true and the EOS
    model did not already compute the baryon density in a column named
    <tt>"nb"</tt>, then \ref bamr::model::compute_eos() should return
    one baryon density and energy density in \ref
    bamr::model::baryon_density_point().

    - If the model provides the symmetry energy and its density
    derivative, it should be stored as constants named <tt>"S"</tt>
    and <tt>"L"</tt> in the table (in \f$ 1/\mathrm{fm} \f$ ).

    - Causality is automatically checked in bamr::compute_star(), but
    the \ref bamr::model::compute_eos() function should check that the
    pressure is not decreasing.

    - Finally, it is recommended to set the interpolation type in the
    \ref o2scl::table_units object to linear interpolation.

    \hline
    \section func_stack_sect Partial Function Call Stack

    The top-level functions in the call stack are:
    - \ref bamr::bamr_class::run()
      - \ref bamr::bamr_class::setup_cli()
      - Command <tt>"model"</tt>: \ref bamr::bamr_class::set_model()
      - Command <tt>"add-data"</tt>: \ref bamr::bamr_class::add_data()
      - Command <tt>"first-point"</tt>: 
        \ref bamr::bamr_class::set_first_point()
      - Command <tt>"mcmc"</tt>: \ref bamr::bamr_class::mcmc()
        - \ref bamr::bamr_class::mcmc_init()
        - \ref bamr::bamr_class::load_mc()
        - \ref bamr::bamr_class::init_grids_table()
          - \ref bamr::bamr_class::table_names_units()
        - Run initial point:
          - \ref bamr::bamr_class::compute_weight() (see below)
          - \ref bamr::bamr_class::add_measurement()
            - \ref bamr::bamr_class::fill_line()
          - \ref bamr::bamr_class::output_best()
        - Main MCMC loop: 
          - If at least one source: \ref bamr::bamr_class::select_mass()
          - \ref bamr::bamr_class::compute_weight() (see below)
          - \ref bamr::bamr_class::make_step()
          - \ref bamr::bamr_class::add_measurement()
            - \ref bamr::bamr_class::fill_line()
          - \ref bamr::bamr_class::output_best()
          - \ref bamr::bamr_class::update_files()
            - If first file update: \ref bamr::bamr_class::first_update()
      - Done with <tt>"mcmc"</tt> command. 

    The operation of \ref bamr::bamr_class::compute_weight() can
    be summarized with:
    - \ref bamr::bamr_class::compute_weight()
      - \ref bamr::bamr_class::compute_star()
        - If the model has an EOS: 
          - \ref bamr::model::compute_eos() to compute the EOS
          - Check pressure is increasing everywhere
          - Compute baryon density if necessary
          - Call \ref bamr::bamr_class::prepare_eos()
          - Compute crust if necessary
          - \ref o2scl::tov_solve::mvsr() to compute the mass-radius curve
          - Check maximum mass
          - If some masses are too large: \ref bamr::bamr_class::select_mass() 
        - Otherwise if there's no EOS: \ref bamr::model::compute_mr()
        - Test for causality

    \hline
    \section postproc_sect Postprocessing code

    The \c process program contains some code to analyze the \c bamr
    output files. Principally, it attempts to remove autocorrelations
    from the data by separating the data into blocks, and using the
    fluctuations in the mean of each block to determine the
    uncertainty in the mean. For example, 

    \code
    ./process -set xscale 197.33 -hist L out.o2 bamr_0_out 
    \endcode

    creates a new file called \c out.o2 with an object of type \ref
    o2scl::table which represents a probability distribution for \f$ L
    \f$ from the \c bamr output file named \c bamr_0_out . The option
    <tt>-set xscale 197.33</tt> ensures that the new table converts
    from \f$ \mathrm{fm}^{-1} \f$ to \f$ \mathrm{Mev} \f$ . The five
    columns are \c reps, \c avgs, \c errs, \c plus and \c minus, which
    give:

    - The central bin values for \f$ L \f$ (in \f$ \mathrm{MeV} \f$)
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
    
    \anchor Bedaque15sv Bedaque15sv:
    <a href="http://dx.doi.org/10.1103/PhysRevLett.114.031103">
    P. Bedaque and A.W. Steiner</a>,
    Phys. Rev. Lett. \b 114 (2015) 031103.

    \anchor Fryer15tf Fryer15tf:
    <a href="http://dx.doi.org/10.1088/0004-637X/812/1/24">
    C.L. Fryer, K. Belczynski, E. Ramirez-Ruiz, S. Rosswog, G. Shen, and 
    A.W. Steiner</a>,
    Astrophys. J. \b 812 (2015) 1.

    \anchor Gandolfi12mm Gandolfi12mm:
    <a href="http://dx.doi.org/10.1103/PhysRevC.85.032801">
    S. Gandolfi, J. Carlson, and S. Reddy</a>
    Phys. Rev. C \b 85 (2012) 032801.

    \anchor Lattimer14co Lattimer14co:
    <a href="http://dx.doi.org/10.1140/epja/i2014-14040-y">
    J.M. Lattimer and A.W. Steiner</a>,
    Eur. Phys. J. A \b 50 (2014) 40.

    \anchor Lattimer14ns Lattimer14ns:
    <a href="http://dx.doi.org/10.1088/0004-637X/784/2/123">
    J.M. Lattimer and A.W. Steiner</a>,
    Astrophys. J. \b 784 (2014) 123.

    \anchor Nattila15eo Nattila15eo:
    <a href="http://arxiv.org/abs/1509.06561">
    J. N&auml;ttil&auml;, A.W. Steiner, J.J.E. Kajava, V.F. Suleimanov, and
    J. Poutanen</a>,
    arxiv.org/1509.06561.

    \anchor Steiner10te Steiner10te:
    <a href="http://dx.doi.org/10.1088/0004-637X/722/1/33">
    A.W. Steiner, J.M. Lattimer, E.F. Brown</a>,
    Astrophys. J. \b 722 (2010) 33.

    \anchor Steiner12cn Steiner12cn:
    <a href="http://dx.doi.org/10.1103/PhysRevLett.108.081102">
    A.W. Steiner and S. Gandolfi</a>,
    Phys. Rev. Lett. \b 108 (2012) 081102.

    \anchor Steiner13tn Steiner13tn:
    <a href="http://dx.doi.org/10.1088/2041-8205/765/1/L5">
    A.W. Steiner, J.M. Lattimer, E.F. Brown</a>,
    Astrophys. J. Lett. \b 765 (2013) 5.

    \anchor Steiner15un Steiner15un:
    <a href="http://dx.doi.org/10.1103/PhysRevC.91.015804">
    A.W. Steiner, S. Gandolfi, F.J. Fattoyev, and W.G. Newton</a>,
    Phys. Rev. C \b 91 (2015) 015804.

    \anchor Steiner16ns Steiner16ns:
    <a href="http://arxiv.org/abs/1510.07515">
    A.W. Steiner, J.M. Lattimer, and E.F. Brown</a>,
    Eur. Phys. J. A (in press).

    \page license_page Licensing
   
    All code is licensed under version 3 of the GPL as provided in the
    files \c COPYING and in \c doc/gpl_license.txt.

    \verbinclude gpl_license.txt

    This documentation is provided under the GNU Free Documentation
    License, as given below and provided in \c
    doc/fdl_license.txt. 
    
    \verbinclude fdl_license.txt
*/
