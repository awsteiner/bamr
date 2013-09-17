/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2013, Andrew W. Steiner
  
  This file is part of Bamr.
  
  Bamr is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.
  
  Bamr is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with Bamr. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/
#ifndef MODELS_H
#define MODELS_H

#include <iostream>

#include <o2scl/cold_nstar.h>
#include <o2scl/schematic_eos.h>
#include <o2scl/gsl_root_brent.h>
#include <o2scl/prob_dens_func.h>

#include "misc.h"
#include "entry.h"

#ifndef DOXYGEN
namespace o2scl {
#endif

  /** \brief Base class for an EOS parameterization
   */
  class model {

  public:

    /// TOV solver and storage for the EOS table
    cold_nstar2 cns;

    model() {
      cns.nb_start=0.01;
    }

    virtual ~model() {}

    /** \brief Set the lower boundaries for all the parameters,
	masses, and radii
    */
    virtual void low_limits(entry &e)=0;

    /** \brief Set the upper boundaries for all the parameters,
	masses, and radii
    */
    virtual void high_limits(entry &e)=0;

    /// Return the name of parameter with index \c i
    virtual std::string param_name(size_t i)=0;

    /** \brief Compute the EOS corresponding to parameters in 
	\c e and put output in \c tab_eos
    */
    virtual void compute_eos(entry &e, bool &fail, std::ofstream &scr_out)=0;

    /// Compute the M-R curve directly
    virtual void compute_mr(entry &e, std::ofstream &scr_out,
			    o2_shared_ptr<table_units<> >::type tab_mvsr,
			    bool &success) {
      return;
    }

    /** \brief A point to calibrate the baryon density with
     */
    virtual void baryon_density_point(double &n1, double &e1) {
      n1=0.0;
      e1=0.0;
      return;
    }

  };

  /** \brief Two polytropes

      Based on the model from \ref Steiner10.
  */
  class two_polytropes : public model {

  protected:

    /// Low-density EOS
    schematic_eos se;

    /// Neutron for \ref se
    fermion neut;

    /// Proton for \ref se
    fermion prot;
    
    /// The fiducial baryon density
    double nb_n1;

    /// The fiducial energy density
    double nb_e1;

  public:

    /** \brief A point to calibrate the baryon density with

        This just returns \ref nb_n1 and \ref nb_e1, which
	are computed in \ref compute_eos().
    */
    virtual void baryon_density_point(double &n1, double &e1) {
      n1=nb_n1;
      e1=nb_e1;
      return;
    }

    /// Create a model object
    two_polytropes();

    virtual ~two_polytropes() {}

    /** \brief Set the lower boundaries for all the parameters,
	masses, and radii
    */
    virtual void low_limits(entry &e);

    /** \brief Set the upper boundaries for all the parameters,
	masses, and radii
    */
    virtual void high_limits(entry &e);

    /// Return the name of parameter with index \c i
    virtual std::string param_name(size_t i);

    /** \brief Compute the EOS corresponding to parameters in 
	\c e and put output in \c tab_eos
    */
    virtual void compute_eos(entry &e, bool &fail, std::ofstream &scr_out);

  };

  /** \brief Alternate polytropes

      As in \ref two_polytropes, but in terms of the exponents instead
      of the polytropic indices.
  */
  class alt_polytropes : public two_polytropes {

  public:

    virtual ~alt_polytropes() {}

    /** \brief Set the lower boundaries for all the parameters,
	masses, and radii
    */
    virtual void low_limits(entry &e);

    /** \brief Set the upper boundaries for all the parameters,
	masses, and radii
    */
    virtual void high_limits(entry &e);

    /// Return the name of parameter with index \c i
    virtual std::string param_name(size_t i);

    /** \brief Compute the EOS corresponding to parameters in 
	\c e and put output in \c tab_eos
    */
    virtual void compute_eos(entry &e, bool &fail, std::ofstream &scr_out);
  
  };

  /** \brief Fixed pressures
    
      Instead of polytropes, linearly interpolate pressures on a fixed
      grid of energy densities.
  */
  class fixed_pressure : public two_polytropes {

  public:

    virtual ~fixed_pressure() {}

    /** \brief Set the lower boundaries for all the parameters,
	masses, and radii
    */
    virtual void low_limits(entry &e);

    /** \brief Set the upper boundaries for all the parameters,
	masses, and radii
    */
    virtual void high_limits(entry &e);

    /// Return the name of parameter with index \c i
    virtual std::string param_name(size_t i);

    /** \brief Compute the EOS corresponding to parameters in 
	\c e and put output in \c tab_eos
    */
    virtual void compute_eos(entry &e, bool &fail, std::ofstream &scr_out);  

  };

  /// A strange quark star model
  class quark_star : public two_polytropes {
  
  public:

    /// The bag constant
    double B;

    /// The paramter controlling non-perturbative corrections to \f$ \mu^4 \f$
    double c;

    /// The gap
    double Delta;

    /// The strange quark mass
    double ms;

    /// The solver to find the chemical potential for zero pressure
    gsl_mroot_hybrids<> gmh;

    /// Desc
    gsl_root_brent<> grb;

    quark_star() {
    }

    virtual ~quark_star() {}
  
    /// Compute the pressure as a function of the chemical potential
    int pressure(size_t nv, const ubvector &x, ubvector &y);

    /// Compute the pressure as a function of the chemical potential
    double pressure2(double mu);

    /** \brief Set the lower boundaries for all the parameters,
	masses, and radii
    */
    virtual void low_limits(entry &e);
  
    /** \brief Set the upper boundaries for all the parameters,
	masses, and radii
    */
    virtual void high_limits(entry &e);
  
    /// Return the name of parameter with index \c i
    virtual std::string param_name(size_t i);
  
    /** \brief Compute the EOS corresponding to parameters in 
	\c e and put output in \c tab_eos
    */
    virtual void compute_eos(entry &e, bool &fail, std::ofstream &scr_out);

  };

  /** \brief Generic quark model

      Alford et al. 2005 parameterizes quark matter with
      \f[
      P = \frac{3 b_4}{4 \pi^2} \mu^4 - \frac{3 b_2}{4 \pi^2} \mu^2 -B 
      \f]
      where \f$ \mu \f$ is the quark chemical potential. QCD corrections 
      can be parameterized by expressing \f$ b_4 \equiv 1-c \f$ , 
      and values of \f$ c \f$ up to 0.4 (or maybe even larger) are
      reasonable (see discussion after Eq. 4 in Alford et al. (2005)).

      The parameter \f$ b_2 = m_s^2 - 4 \Delta^2 \f$ for CFL quark
      matter, and can thus be positive or negative. A largest possible
      range might be somewhere between \f$ (400~\mathrm{MeV})^2 \f$,
      which corresponds to the situation where the gap is zero and the
      strange quarks receive significant contributions from chiral
      symmetry breaking, to \f$ (150~\mathrm{MeV})^2-4
      (200~\mathrm{MeV})^2 \f$ which corresponds to a bare strange
      quark with a large gap. In units of \f$ \mathrm{fm}^{-1} \f$ ,
      this corresponds to a range of about \f$ -3.5 \f$ to \f$
      4~\mathrm{fm}^{-2} \f$ . In Alford et al. (2010), they choose a
      significantly smaller range, from \f$ -1 \f$ to \f$
      1~\mathrm{fm}^{-2} \f$.
      
      Simplifying the parameterization to 
      \f[
      P = a_4 \mu^4 +a_2 \mu^2 - B 
      \f]
      gives the following ranges
      \f[
      a_4 = 0.045~\mathrm{to}~0.08
      \f]
      and 
      \f[
      a_2 = -0.3~\mathrm{to}~0.3~\mathrm{fm}^{-2}
      \f]
      for the "largest possible range" described above or
      \f[
      a_2 = -0.08~\mathrm{to}~0.08~\mathrm{fm}^{-2}
      \f]
      for the range used by Alford et al. (2010). 
      
      The energy density is
      \f[
      \varepsilon = B + a_2 \mu^2 + 3 a_4 \mu^4
      \f]
    
      Note that 
      \f{eqnarray*}
      \frac{dP}{d \mu} &=& 2 a_2 \mu + 4 a_4 \mu^3 \nonumber \\ 
      \frac{d\varepsilon}{d \mu} &=& 2 a_2 \mu + 12 a_4 \mu^3
      \f}
  */
  class generic_quarks : public two_polytropes {
  
  public:
  
    virtual ~generic_quarks() {}

    /** \brief Set the lower boundaries for all the parameters,
	masses, and radii
    */
    virtual void low_limits(entry &e);

    /** \brief Set the upper boundaries for all the parameters,
	masses, and radii
    */
    virtual void high_limits(entry &e);

    /// Return the name of parameter with index \c i
    virtual std::string param_name(size_t i);

    /** \brief Compute the EOS corresponding to parameters in 
	\c e and put output in \c tab_eos
    */
    virtual void compute_eos(entry &e, bool &fail, std::ofstream &scr_out);

  };

#ifndef DOXYGEN
}
#endif

#endif
