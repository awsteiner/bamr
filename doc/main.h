/** \file main.h
    \brief File containing user's guide documentation
*/
/** \mainpage User's Guide
    
    This document describes the open-source MPI implementation of a
    Bayesian analysis of mass and radius data to determine the mass
    versus radius curve and the equation of state of dense matter.
    This package will principally be useful for those physicists and
    astrophysicists who are already familiar with C++ and are
    interested in modifying this code for their own use.

    This implementation was originally supported by Chandra grant
    TM1-12003X.

    This is a beta version. Use at your own risk.

    Currently, \bm is dual-hosted as an SVN respostory at
    http://www.sourceforge.net/projects/bamr and a git repository at
    http://www.github.com/awsteiner/bamr . This documentation (when it
    corresponds to a release) is hosted at http://bamr.sourceforge.net
    .

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
    - \ref ack_sect
    - \ref ref_sect
    - \ref license_page

    \hline
    \section install_sect Installation
    
    The \bm executable requires the installation of 
    <a href="http://www.gnu.org/software/gsl">GSL</a> 
    (versions 1.15 and later),     
    <a href="http://www.hdfgroup.org">HDF5</a> 
    (versions 1.8.4 and later), 
    \htmlonly
    <a href="http://o2scl.sourceforge.net">O<span style='position: relative; top: 0.3em; font-size: 0.8em'>2</span>scl</a> (version 0.913 only), 
    \endhtmlonly
    \latexonly
    O$_2$scl
    \endlatexonly
    and MPI (tested with openmpi-1.4.2). After these four packages are
    successfully installed, you will need to edit the \c makefile and
    then compile \bm before execution.

    \hline
    \section dl_sect Download

    The most recent .tar.gz file can be obtained from
    http://sourceforge.net/projects/bamr/files/latest/download
    or you may check out the most recent release version
    from the SVN repository using
    \verbatim
    svn checkout svn://svn.code.sf.net/p/bamr/code/trunk bamr
    \endverbatim

    \hline
    \section usage_sect Basic Usage
    
    The basic usage is something like
    \verbatim
    mpirun -np 4 ./bamr -model twop -run default.in -mcmc run1
    \endverbatim
    to perform a one day run with model \c twop with the input
    file in \c default.in. 

    There are several variables which can be modified with the
    \c set command, e.g. 
    \verbatim
    ./bamr -model twop -set max_time 43200 -run default.in -mcmc run2
    \endverbatim

    An example of an MPI invocation is
    \verbatim
    mpirun -np 4 ./bamr -set model twop -run default.in -mcmc run3 &
    \endverbatim
    which runs with four processors on the current machine.

    Also try
    \verbatim
    ./bamr -help
    \endverbatim
    which outputs some additional information on the 
    relevant functions and parameters. 

    \hline
    \section infile_sect Data files

    The data files are HDF5 files (typically named with a <tt>.o2</tt>
    extension) which contain one \ref table3d object giving the
    probability density of a neutron star observation as a slice in
    that table.

    \hline
    \section outfile_sect Output Files

    Output is stored in HDF files with a prefix given by the
    argument to the \c mcmc command, one set of files
    for each processor. Output includes files with the 
    following suffixes (where X is the processor index):
    - \c _X_best: "Best" point found during run
    - \c _X_out: Ensemble of histograms of EOS relative to SLy4
    - \c _X_scr: Running output of entire simulation

    If the executable is run directly (without <tt>mpirun</tt>)
    then X is always zero.

    \hline
    \section detail_sect Some Details

    The basic functionality is provided in the \ref bamr class and
    each Monte Carlo point is an object of type \ref entry. All of the
    "models" (EOS parameterizations) are children of \ref model class.

    If the initial guess has no probability, then the code will fail.
    This is indicated by the line \c "Initial weight zero." in
    the \c _scr file. The easiest fix is just to change the initial 
    guess.

    \comment
    While the current code available is designed to minimize crashes,
    they do occasionally occur. The TOV solver cannot handle equations
    of state which diverge or are strongly discontinuous. Errors 
    of the form
    \verbatim
    terminate called after throwing an instance of
    'o2scl::exc_overflow_error'
    what():  Error ezerodiv in file tridiag.c at line 117.
    matrix must be positive definite
    \endverbatim
    usually come from the interpolation routine that the TOV solver is
    using to interpolate the EOS. (This particular error is actually a
    C++ exception thrown inside GSL error since O2scl has replaced the
    default GSL error handler.) The best way to debug this is to
    examine the \c _scr file to find the last parameter set which the
    TOV solver tried to evaluate and examine the EOS. One can often do
    this by using the potentially difficult parameter set as the
    initial guess in the input file (i.e. \c default.in) 
    and setting \c debug_eos to \c 1 (true) to force the
    program to output the EOS which it is sending to the TOV solver.
    \endcomment

    \hline
    \section model_sect EOS Model

    Some EOS models are already provided. New models (i.e. new
    children of the \ref model class) must perform several tasks

    - The function \ref model::compute_eos() should use the parameters
    in the \ref entry argument to compute the EOS and store it in the
    object returned by \ref cold_nstar::get_eos_results().

    - The energy density should be stored in a column named
    <tt>ed</tt> and the pressure in <tt>pr</tt> with the correct units
    set for each column (currently only <tt>1/fm^4</tt> is supported).

    - If \ref bamr::baryon_density is true and the EOS model did not
    already compute the baryon density in a column named <tt>"nb"</tt>,
    then \ref model::compute_eos() should return one baryon density
    and energy density in \ref model::baryon_density_point().

    - If the model provides the symmetry energy and its density
    derivative, it should be stored as constants named <tt>"S"</tt>
    and <tt>"L"</tt> in the table (in \f$ 1/\mathrm{fm} \f$ ).

    - Causality is automatically checked in bamr::compute_star(), but
    the \ref model::compute_eos() function should check that the
    pressure is not decreasing. 

    - Finally, it is recommended to set the interpolation type in the
    \ref table_units object to linear interpolation.

    \comment
    \hline
    \section plot_sect Plotting

    Plotting requires the installation of 
    <a href="http://root.cern.ch">ROOT</a>, but only requires the most
    basic functionality of ROOT so all of the extra packages which are
    available can be disabled. The relevant \c makefile targets
    are \c plot and \c plot2d.
    \endcomment

    \hline
    \section ack_sect Acknowledgements

    I would like to thank Ed Brown, Stefano Gandolfi, and Jim 
    Lattimer for their collaboration on these projects.

    \hline
    \section ref_sect Bibliography

    Some of the references which contain links should direct you to
    the work referred to using its DOI identifer.
    
    \anchor Lattimer13 Lattimer13:
    \htmlonly
    <a href="http://arxiv.org/abs/">
    J.M. Lattimer and A.W. Steiner</a>,
    \endhtmlonly
    \latexonly
    \endlatexonly
    Astrophys. J., submitted (2013).

    \anchor Steiner10 Steiner10:
    \htmlonly
    <a href="http://dx.doi.org/10.1088/0004-637X/722/1/33">
    A.W. Steiner, J.M. Lattimer, E.F. Brown</a>,
    \endhtmlonly
    \latexonly
    \href{http://dx.doi.org/10.1088/0004-637X/722/1/33}{
    A.~W. Steiner, J.~M. Lattimer, E.~F. Brown},
    \endlatexonly
    Astrophys. J. \b 722 (2010) 33.

    \anchor Steiner12 Steiner12:
    \htmlonly
    <a href="http://dx.doi.org/10.1103/PhysRevLett.108.081102">
    A.W. Steiner and S. Gandolfi</a>,
    \endhtmlonly
    \latexonly
    \href{http://dx.doi.org/10.1103/PhysRevLett.108.081102}{
    A.~W. Steiner and S. Gandolfi},
    \endlatexonly
    Phys. Rev. Lett. \b 108 (2012) 081102.

    \anchor Steiner13 Steiner13:
    \htmlonly
    <a href="http://dx.doi.org/10.1088/2041-8205/765/1/L5">
    A.W. Steiner, J.M. Lattimer, E.F. Brown</a>,
    \endhtmlonly
    \latexonly
    \href{http://dx.doi.org/10.1088/2041-8205/765/1/L5}{
    A.~W. Steiner, J.~M. Lattimer, E.~F. Brown},
    \endlatexonly
    Astrophys. J. Lett. \b 765 (2013) 5.

    \page license_page Licensing
   
    All code is licensed under version 3 of the GPL as provided in the
    files \c COPYING and in \c doc/gpl_license.txt.

    \verbinclude gpl_license.txt

    This documentation is provided under the GNU Free Documentation
    License, as given below and provided in \c
    doc/fdl_license.txt. 
    
    \verbinclude fdl_license.txt
*/
