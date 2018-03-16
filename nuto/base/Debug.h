#pragma once

#ifdef NDEBUG
#define DBG_POSITION ;
#define DBG_POSITION_INFO(v) ;
#define DBG_PRINT_VAL(v) ;
#define DBG_PRINT_VEC(v) ;

#else // NDEBUG

#include <iostream>
#include <boost/foreach.hpp>
#include <boost/typeof/typeof.hpp>

#define DBG_OPENER std::cout << __FILE__ << ":" << __LINE__ << ":\t"

#define DBG_POSITION DBG_OPENER << std::endl;

#define DBG_POSITION_INFO(v) DBG_OPENER << "'\"" << (v) << "\"" << std::endl;

#define DBG_PRINT_VAL(v) DBG_OPENER << "'" << #v << "' = \"" << (v) << "\"" << std::endl << std::endl;

#define DBG_PRINT_VEC(v)                                                                                               \
    DBG_OPENER << "Entries of: \"" << #v << "\"" << std::endl;                                                         \
    {                                                                                                                  \
        std::size_t dbg_i(0);                                                                                          \
        BOOST_FOREACH (BOOST_TYPEOF(*boost::begin(v)) const& x, v)                                                     \
            std::cout << #v << "[" << dbg_i++ << "] : \"" << x << "\"" << std::endl;                                   \
    }                                                                                                                  \
    std::cout << std::endl << std::endl;

#endif // NDEBUG
