/*
	Copyright (C) 2012  Andrew Cotter

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


/**
	\file fast_exp.hpp
	\brief FastExp function

	<b>Nicol N. Schraudolph. "A Fast, Compact Approximation of the Exponential
	Function". In <em>Neural Computation</em>, 11(4). 1998.</b>
*/




#ifndef __FAST_EXP_HPP__
#define __FAST_EXP_HPP__

#ifdef __cplusplus




#include <boost/cstdint.hpp>
#include <boost/static_assert.hpp>

#include <limits>

#include <math.h>




namespace _Private {




//============================================================================
//    FastExpHelper class
//============================================================================


// **NOTE: easy enough to do a FastLog, but it isn't any faster than log
template< size_t t_Size >
struct FastExpHelper;


template<>
struct FastExpHelper< 4 > {

	BOOST_STATIC_ASSERT( sizeof( double  ) == 8 );
	BOOST_STATIC_ASSERT( sizeof( int32_t ) == 4 );


	static inline double Exp( double exponent );


private:

	union Union {

		double d;
#if defined( BIG_ENDIAN )
		struct { int32_t lower, upper; };
#elif defined ( LITTLE_ENDIAN )
		struct { int32_t upper, lower; };
#else    // BIG_ENDIAN / LITTLE_ENDIAN
		BOOST_STATIC_ASSERT( false );
#endif    // BIG_ENDIAN / LITTLE_ENDIAN
	};
};


template<>
struct FastExpHelper< 8 > {

	BOOST_STATIC_ASSERT( sizeof( double  ) == 8 );
	BOOST_STATIC_ASSERT( sizeof( int64_t ) == 8 );


	static inline double Exp( double exponent );


private:

	union Union {

		double d;
		int64_t l;
	};
};




//============================================================================
//    FastExpHelper< 4 > inline methods
//============================================================================


double FastExpHelper< 4 >::Exp( double exponent ) {

	static double const scale = ( 1ull << 20 ) / logl( 2.0 );
	static int32_t const shift = roundl( ( 1ull << 20 ) * ( 0x03ff - log2l( 0.5 + ( 1.0 / ( expl( 1.0 ) * logl( 2.0 ) ) ) ) ) );

	// 709 = floor( 2^10 * log( 2 ) )
	if ( exponent < -709 )
		return 0;
	else if ( exponent > 709 )
		return std::numeric_limits< double >::infinity();
	else {

		Union u;
		u.lower = 0;
		u.upper = exponent * scale + 0.5;
		u.upper += shift;
		return u.d;
	}
}




//============================================================================
//    FastExpHelper< 8 > inline methods
//============================================================================


double FastExpHelper< 8 >::Exp( double exponent ) {

	/*
		**FIXME:
			Does this still work theoretically (it does, practically) if shift
			is a double? It's much faster that way on 32-bit architectures
	*/

	static double const scale = ( 1ull << 52 ) / logl( 2.0 );
	static int64_t const shift = roundl( ( 1ull << 52 ) * ( 0x03ff - log2l( 0.5 + ( 1.0 / ( expl( 1.0 ) * logl( 2.0 ) ) ) ) ) );

	// 709 = floor( 2^10 * log( 2 ) )
	if ( exponent < -709 )
		return 0;
	else if ( exponent > 709 )
		return std::numeric_limits< double >::infinity();
	else {

		Union u;
		u.l = exponent * scale + 0.5;
		u.l += shift;
		return u.d;
	}
}




}    // namespace _Private




//============================================================================
//    FastExp function
//============================================================================


inline double FastExp( double exponent ) {

	return _Private::FastExpHelper< sizeof( long ) >::Exp( exponent );
}




#endif    /* __cplusplus */

#endif    /* __FAST_EXP_HPP__ */
