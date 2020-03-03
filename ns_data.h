/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2019, Andrew W. Steiner
  
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
/** \file ns_data.h
    \brief Definition of \ref ns_data
*/
#ifndef NS_DATA_H
#define NS_DATA_H

#include <iostream>

#include <o2scl/table3d.h>

#include "settings.h"

namespace bamr {
  
  /** \brief Neutron star data object

      This class is designed so that multiple OpenMP threads can use
      the same const instance of this class (i.e. so long as they do
      not change the member data).

      \future Maybe it would be better to restructure this 
      object rather than having many vectors of the same size.
   */
  class ns_data {

  public:

    ns_data() {
      n_sources=0;
    }      

    virtual ~ns_data() {
    }
    
    /// \name Input neutron star data
    //@{
    /// Input probability distributions
    std::vector<o2scl::table3d> source_tables;

    /// Alternate input probability distributions
    std::vector<o2scl::table3d> source_tables_alt;

    /// The names for each source
    std::vector<std::string> source_names;

    /// The names of the table in the data file
    std::vector<std::string> table_names;

    /// File names for each source
    std::vector<std::string> source_fnames;

    /// Alternate file names for each source
    std::vector<std::string> source_fnames_alt;

    /// Slice names for each source
    std::vector<std::string> slice_names;

    /// The initial set of neutron star masses
    std::vector<double> init_mass_fracs;

    /** \brief The number of sources
     */
    size_t n_sources;

    /** \brief Add a data distribution to the list
     */
    virtual int add_data(std::vector<std::string> &sv, bool itive_com);

    /** \brief Add a data distribution to the list
     */
    virtual int add_data_alt(std::vector<std::string> &sv, bool itive_com);

    /** \brief Load input probability distributions
	
	Ensure all MPI ranks read all files while ensuring that
	no two ranks simultaneously read the same file.
	
	Let MPI_Size be denoted n. If n is 1, no messages are sent and
	the single rank always reads the files in order. When n is
	larger than 1, rank 0 still reads the files in order. Rank n-1
	begins by reading the second file, rank n-2 begins by reading
	the third file, and so on, until there are no more ranks left
	which are not reading files or there are no more files left.
	After reading, all ranks which read files send a message to
	the next rank (in sequence) prompting them to proceed with the
	next file. Once each rank has read all of the files, it does
	not send or receive any further messages (the number of
	messages sent for n>1 is n times the number of files times).
	If \ref settings::mpi_load_debug is set to 1, then the various
	MPI messages are copied to scr_out, no files are actually
	read, and exit() is called after all messages have been sent
	and received by all ranks.
    */
    virtual void load_mc(std::ofstream &scr_out, int mpi_nprocs,
			 int mpi_rank, std::shared_ptr<settings> set);
    //@}

  };
  
}

#endif
