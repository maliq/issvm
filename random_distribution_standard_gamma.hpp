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
	\file random_distribution_standard_gamma.hpp
	\brief StandardGamma implementation
*/




#ifndef __RANDOM_DISTRIBUTION_STANDARD_GAMMA_HPP__
#define __RANDOM_DISTRIBUTION_STANDARD_GAMMA_HPP__

#ifdef __cplusplus




#include "random_distribution_standard_gaussian.hpp"
#include "random_distribution_standard_exponential.hpp"
#include "fast_exp.hpp"
#include "helpers.h"

#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <limits>
#include <cmath>




namespace Random {


namespace Distribution {




//============================================================================
//    StandardGamma class
//============================================================================


/*
	George Marsaglia, Wai Wan Tsang. "A Simple Method for Generating Gamma
	Variables"

	To sample from Y~Gamma(shape,scale), just let X~Gamma(shape) and take Y=X*scale
*/
template< typename t_Type = double >
struct StandardGamma {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );

	typedef t_Type Type;


	inline StandardGamma( t_Type const shape );


	inline t_Type const GetShape() const;


	inline double const PDF( t_Type const value ) const;
	inline double const CDF( t_Type const value ) const;

	template< typename t_Generator >
	inline t_Type const Sample( t_Generator& generator );


private:

	t_Type m_shape;
	t_Type m_dd;
	t_Type m_cc;


	static StandardExponential< t_Type > s_standardExponential;
	static StandardGaussian< t_Type > s_standardGaussian;
};


template< typename t_Type >
StandardExponential< t_Type > StandardGamma< t_Type >::s_standardExponential;


template< typename t_Type >
StandardGaussian< t_Type > StandardGamma< t_Type >::s_standardGaussian;




//============================================================================
//    StandardGamma inline methods
//============================================================================


template< typename t_Type >
StandardGamma< t_Type >::StandardGamma( t_Type const shape ) :
	m_shape( shape ),
	m_dd( ( ( m_shape < 1.0 ) ? ( m_shape + 1.0 ) : m_shape ) - ( 1.0 / 3.0 ) ),
	m_cc( ( 1.0 / 3.0 ) / std::sqrt( m_dd ) )
{
	BOOST_ASSERT( m_shape > 0 );
}


template< typename t_Type >
t_Type const StandardGamma< t_Type >::GetShape() const {

	return m_shape;
}


template< typename t_Type >
double const StandardGamma< t_Type >::PDF( t_Type const value ) const {

	/*
		Might be better to use:

		Catherine Loader. "Fast and Accurate Computation of Binomial
		Probabilities"
	*/
	return( ( value > 0 ) ? std::exp( ( m_shape - 1 ) * std::log( value ) - value - lgamma( m_shape ) ) : 0 );
}


template< typename t_Type >
double const StandardGamma< t_Type >::CDF( t_Type const value ) const {

	// **TODO **FIXME
	throw std::runtime_error( "StandardGamma::CDF is unimplemented" );
}


template< typename t_Type >
template< typename t_Generator >
t_Type const StandardGamma< t_Type >::Sample( t_Generator& generator ) {

#ifndef NDEBUG
	t_Type variate = std::numeric_limits< t_Type >::quiet_NaN();
#else    // NDEBUG
	t_Type variate;
#endif    // NDEBUG

	for ( ; ; ) {

		t_Type gaussian;
		t_Type vv;
		do {

			gaussian = s_standardGaussian.Sample( generator );
			vv = gaussian * m_cc + 1.0;

		} while ( vv <= 0.0 );
		vv = Cube( vv );

		typename t_Generator::FloatingType const uniform = generator.SampleStandardUniform();
		t_Type const gaussianSquared = Square( gaussian );
		if ( uniform < 1.0 - 0.0331 * Square( gaussianSquared ) ) {

			variate = m_dd * vv;
			break;    // **NOTE: break
		}
		if ( std::log( uniform ) < 0.5 * gaussianSquared + m_dd * ( 1.0 - vv + std::log( vv ) ) ) {

			variate = m_dd * vv;
			break;    // **NOTE: break
		}
	}

	if ( m_shape < 1.0 )
		variate *= FastExp( s_standardExponential.Sample( generator ) / -m_shape );    // this happens frequently enough that FastExp is probably worth it
		//variate *= std::pow( generator.SampleStandardUniform(), 1.0 / m_shape );

	BOOST_ASSERT( std::isfinite( variate ) );
	return variate;
}




}    // namespace Distribution


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_DISTRIBUTION_STANDARD_GAMMA_HPP__ */
