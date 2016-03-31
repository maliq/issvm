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
	\file random_distribution_standard_exponential.hpp
	\brief StandardExponential implementation
*/




#ifndef __RANDOM_DISTRIBUTION_STANDARD_EXPONENTIAL_HPP__
#define __RANDOM_DISTRIBUTION_STANDARD_EXPONENTIAL_HPP__

#ifdef __cplusplus




#include "random_distribution_private_ziggurat.hpp"
#include "helpers.h"

#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <limits>
#include <cmath>




namespace Random {


namespace Distribution {


namespace {




//============================================================================
//    StandardExponentialZiggurat class
//============================================================================


/*
	George Marsaglia, Wai Wan Tsang. "The Ziggurat Method for Generating
	Random Variables"
*/
template<
	typename t_Type,
	unsigned char t_LogLevels = 8
>
struct StandardExponentialZiggurat {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );

	typedef t_Type Type;


	inline StandardExponentialZiggurat();


	template< typename t_Generator >
	inline t_Type const Sample( t_Generator& generator );


private:

	enum {
		LOG_LEVELS = t_LogLevels,
		LEVELS = ( 1 << LOG_LEVELS )
	};
	BOOST_STATIC_ASSERT( LOG_LEVELS >= 1 );
	BOOST_STATIC_ASSERT( LEVELS     >= 2 );


	static inline long double const Density( long double const xx );
	static inline long double const InverseDensity( long double const yy );
	static inline long double const TailArea( long double const xx );


	t_Type m_xx[ LEVELS + 1 ];
	t_Type m_yy[ LEVELS + 1 ];
};




//============================================================================
//    StandardExponentialZiggurat inline methods
//============================================================================


template< typename t_Type, unsigned char t_LogLevels >
StandardExponentialZiggurat< t_Type, t_LogLevels >::StandardExponentialZiggurat() {

	// initialize the exponential ziggurat
	_Private::InitializeZiggurat(
		&m_xx[ 0 ],
		&m_yy[ 0 ],
		LEVELS,
		Density,
		InverseDensity,
		TailArea
	);
}


template< typename t_Type, unsigned char t_LogLevels >
template< typename t_Generator >
t_Type const StandardExponentialZiggurat< t_Type, t_LogLevels >::Sample( t_Generator& generator ) {

	t_Type variate = 0;

	for ( ; ; ) {

		BOOST_STATIC_ASSERT( sizeof( typename t_Generator::DiscreteType ) * 8 >= LOG_LEVELS );
		typename t_Generator::DiscreteType index = ( generator.SampleDiscreteUniform() >> ( sizeof( typename t_Generator::DiscreteType ) * 8 - LOG_LEVELS ) );
		BOOST_ASSERT( index >= 0      );
		BOOST_ASSERT( index <  LEVELS );

		t_Type const xx = generator.SampleStandardUniform() * m_xx[ index + 1 ];
		if ( LIKELY( std::fabs( xx ) < m_xx[ index ] ) ) {    // inside a rectangle

			variate += xx;
			break;    // **NOTE: break
		}
		// nothing below this point needs to be fast--it will run infrequently
		else if ( index < LEVELS - 1 ) {    // rejection sample from a corner

			t_Type const yy = m_yy[ index + 1 ] + generator.SampleStandardUniform() * ( m_yy[ index ] - m_yy[ index + 1 ] );
			if ( yy < std::exp( -xx ) ) {    // **NOTE: no need to use FastExp, since this happens so infrequently

				variate += xx;
				break;    // **NOTE: break
			}

			// **NOTE: continue
		}
		else {    // sample from the tail

			BOOST_ASSERT( index == LEVELS - 1 );
			t_Type const radius = m_xx[ LEVELS - 1 ];

			/*
				We could instead do the following:
					variate = radius - std::log( generator.SampleStandardUniform() );
					break;    // **NOTE: break
				and change "variate +=" to "variate =" everywhere
			*/
			variate += radius;

			// **NOTE: continue
		}
	}

	return variate;
}




//============================================================================
//    StandardExponentialZiggurat static methods
//============================================================================


template< typename t_Type, unsigned char t_LogLevels >
long double const StandardExponentialZiggurat< t_Type, t_LogLevels >::Density( long double const xx ) {

	BOOST_ASSERT( xx >= 0.0l );
	return std::exp( -xx );
}


template< typename t_Type, unsigned char t_LogLevels >
long double const StandardExponentialZiggurat< t_Type, t_LogLevels >::InverseDensity( long double const yy ) {

	BOOST_ASSERT( yy >= 0.0l );
	BOOST_ASSERT( yy <= 1.0l );
	return -std::log( yy );
}


template< typename t_Type, unsigned char t_LogLevels >
long double const StandardExponentialZiggurat< t_Type, t_LogLevels >::TailArea( long double const xx ) {

	BOOST_ASSERT( xx >= 0.0l );
	return std::exp( -xx );
}




}    // anonymous namespace




//============================================================================
//    StandardExponential class
//============================================================================


// To sample from Y~Exp(rate), just let X~Exp() and take Y=X/rate
template< typename t_Type = double >
struct StandardExponential {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );

	typedef t_Type Type;


	inline StandardExponential();


	inline double const PDF( t_Type const value ) const;
	inline double const CDF( t_Type const value ) const;

	template< typename t_Generator >
	inline t_Type const Sample( t_Generator& generator );


private:

	static StandardExponentialZiggurat< t_Type > s_ziggurat;
};


template< typename t_Type >
StandardExponentialZiggurat< t_Type > StandardExponential< t_Type >::s_ziggurat;




//============================================================================
//    StandardExponential inline methods
//============================================================================


template< typename t_Type >
StandardExponential< t_Type >::StandardExponential() {
}


template< typename t_Type >
double const StandardExponential< t_Type >::PDF( t_Type const value ) const {

	return( ( value > 0 ) ? std::exp( -value ) : 0 );
}


template< typename t_Type >
double const StandardExponential< t_Type >::CDF( t_Type const value ) const {

	return( ( value > 0 ) ? ( 1 - std::exp( -value ) ) : 0 );
}


template< typename t_Type >
template< typename t_Generator >
t_Type const StandardExponential< t_Type >::Sample( t_Generator& generator ) {

	return s_ziggurat.Sample( generator );
}




}    // namespace Distribution


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_DISTRIBUTION_STANDARD_EXPONENTIAL_HPP__ */
