// $Id$

#ifndef NUTO_DEBUG_H
#define NUTO_DEBUG_H

#ifdef DEBUG

#include <iostream>

#include <boost/typeof/typeof.hpp>
#include <boost/foreach.hpp>

//! @author Daniel Arnold, ISM
//! @date April 2010
//! @brief ... some macros for debugging

	#define DBG_OPENER				std::cout << __FILE__ << ":" << __LINE__ << ":\t"

	#define DBG_POSITION			DBG_OPENER << std::endl;

	#define DBG_POSITION_INFO(v)	DBG_OPENER<< "'\"" << (v) << "\"" << std::endl;

	#define DBG_PRINT_VAL(v)		DBG_OPENER << "'" << #v << "' = \"" << (v) << "\"" << std::endl << std::endl;

	#define DBG_PRINT_VEC(v)		DBG_OPENER << "Entries of: \"" << #v << "\""<< std::endl; \
									{std::size_t dbg_i(0); \
									BOOST_FOREACH( BOOST_TYPEOF( *boost::begin(v) ) const& x , v) \
										std::cout << #v << "["<<dbg_i++<<"] : \"" << x<<"\""<< std::endl;} \
									std::cout<<std::endl<<std::endl;



#else //DEBUG
	#define DBG_POSITION ;
	#define DBG_POSITION_INFO(v) ;
	#define DBG_PRINT_VAL(v) ;
	#define DBG_PRINT_VEC(v) ;
#endif //DEBUG


#endif //NUTO_DEBUG_H
