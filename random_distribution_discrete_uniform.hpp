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
	\file random_distribution_discrete_uniform.hpp
	\brief DiscreteUniform implementation
*/




#ifndef __RANDOM_DISTRIBUTION_DISCRETE_UNIFORM_HPP__
#define __RANDOM_DISTRIBUTION_DISCRETE_UNIFORM_HPP__

#ifdef __cplusplus




#include "helpers.h"

#include <boost/type_traits.hpp>
#include <boost/cstdint.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <limits>




namespace Random {


namespace Distribution {




//============================================================================
//    DiscreteUniform class
//============================================================================


template< typename t_Type = unsigned int >
struct DiscreteUniform {

	BOOST_STATIC_ASSERT( boost::is_integral< t_Type >::value && boost::is_unsigned< t_Type >::value );

	typedef t_Type Type;


	inline DiscreteUniform();                        // uniform on t_Type
	inline DiscreteUniform( t_Type const bound );    // uniform on [0,bound]


	inline t_Type const GetBound() const;


	inline double const PDF( t_Type const value ) const;
	inline double const CDF( t_Type const value ) const;

	template< typename t_Generator >
	inline t_Type const Sample( t_Generator& generator );


private:

	t_Type m_bound;
	unsigned int m_lowBits;
};




//============================================================================
//    DiscreteUniform inline methods
//============================================================================


template< typename t_Type >
DiscreteUniform< t_Type >::DiscreteUniform() :
	m_bound( -1 ),
	m_lowBits( sizeof( t_Type ) * 8 )
{
}


template< typename t_Type >
DiscreteUniform< t_Type >::DiscreteUniform( t_Type const bound ) :
	m_bound( bound ),
	m_lowBits( CountLowBits( bound ) )
{
}


template< typename t_Type >
t_Type const DiscreteUniform< t_Type >::GetBound() const {

	return m_bound;
}


template< typename t_Type >
double const DiscreteUniform< t_Type >::PDF( t_Type const value ) const {

	return( ( value <= m_bound ) ? ( 1 / ( static_cast< double >( m_bound ) + 1 ) ) : 0 );
}


template< typename t_Type >
double const DiscreteUniform< t_Type >::CDF( t_Type const value ) const {

	return( ( value <= m_bound ) ? ( static_cast< double >( value + 1 ) / ( static_cast< double >( m_bound ) + 1 ) ) : 1 );
}


template< typename t_Type >
template< typename t_Generator >
t_Type const DiscreteUniform< t_Type >::Sample( t_Generator& generator ) {

	BOOST_ASSERT( sizeof( typename t_Generator::DiscreteType ) * 8 >= m_lowBits );

	// scaling would result in a non-uniform distribution, so we use rejection sampling
	t_Type result;
	do {
		// **NOTE: we need both shift and mask to handle the case in which m_lowBits = 0 (for which the shift becomes a NOP)
		const unsigned int shift = ( sizeof( typename t_Generator::DiscreteType ) * 8 - m_lowBits );
		const typename t_Generator::DiscreteType mask = ( static_cast< typename t_Generator::DiscreteType >( 1 ) << m_lowBits ) - 1;
		BOOST_ASSERT( CountLowBits( mask ) == m_lowBits );
		result = ( ( generator.SampleDiscreteUniform() >> shift ) & mask );
	} while ( result > m_bound );

	return result;
}




}    // namespace Distribution


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_DISTRIBUTION_DISCRETE_UNIFORM_HPP__ */
