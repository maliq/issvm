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
	\file random_generator_lagged_fibonacci_4.hpp
	\brief LaggedFibonacci4 class
*/




#ifndef __RANDOM_GENERATOR_LAGGED_FIBONACCI_4_HPP__
#define __RANDOM_GENERATOR_LAGGED_FIBONACCI_4_HPP__

#ifdef __cplusplus




#include "random_generator_private_helpers.hpp"
#include "random_generator_private_lagged_fibonacci_4_helper.hpp"
#include "helpers.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <limits>
#include <cmath>




namespace Random {


namespace Generator {




//============================================================================
//    LaggedFibonacci4 class
//============================================================================


template<
	typename t_DiscreteType = unsigned int,
	typename t_FloatingType = double
>
struct LaggedFibonacci4 {

	BOOST_STATIC_ASSERT(
		(
			boost::is_integral< t_DiscreteType >::value &&
			boost::is_unsigned< t_DiscreteType >::value
		) &&
		boost::is_floating_point< t_FloatingType >::value
	);

	typedef t_DiscreteType DiscreteType;
	typedef t_FloatingType FloatingType;

	enum {
		STATE_SIZE = (
			_Private::LaggedFibonacci4Helper< t_DiscreteType >::STATE_SIZE +
			_Private::LaggedFibonacci4Helper< t_FloatingType >::STATE_SIZE
		)
	};


	inline LaggedFibonacci4();
	inline explicit LaggedFibonacci4( LaggedFibonacci4 const& other );


	inline void Seed();    // from /dev/urandom

	template< typename t_Generator >
	inline void Seed( t_Generator& generator );

	template< typename t_Iterator >
	inline void Seed( t_Iterator const& iiBegin, t_Iterator const& iiEnd );


	inline t_DiscreteType const SampleDiscreteUniform();
	inline t_FloatingType const SampleStandardUniform();


private:

	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	_Private::LaggedFibonacci4Helper< t_DiscreteType > m_discreteHelper;
	_Private::LaggedFibonacci4Helper< t_FloatingType > m_floatingHelper;


	friend class boost::serialization::access;
};




//============================================================================
//    LaggedFibonacci4 inline methods
//============================================================================


template< typename t_DiscreteType, typename t_FloatingType >
LaggedFibonacci4< t_DiscreteType, t_FloatingType >::LaggedFibonacci4() {
}


template< typename t_DiscreteType, typename t_FloatingType >
LaggedFibonacci4< t_DiscreteType, t_FloatingType >::LaggedFibonacci4( LaggedFibonacci4 const& other ) :
	m_discreteHelper( other.m_discreteHelper ),
	m_floatingHelper( other.m_floatingHelper )
{
}


template< typename t_DiscreteType, typename t_FloatingType >
void LaggedFibonacci4< t_DiscreteType, t_FloatingType >::Seed() {

	m_discreteHelper.Seed();
	m_floatingHelper.Seed();
}


template< typename t_DiscreteType, typename t_FloatingType >
template< typename t_Generator >
void LaggedFibonacci4< t_DiscreteType, t_FloatingType >::Seed( t_Generator& generator ) {

	m_discreteHelper.Seed( generator );
	m_floatingHelper.Seed( generator );
}


template< typename t_DiscreteType, typename t_FloatingType >
template< typename t_Iterator >
void LaggedFibonacci4< t_DiscreteType, t_FloatingType >::Seed( t_Iterator const& iiBegin, t_Iterator const& iiEnd ) {

	t_Iterator const iiMiddle = iiBegin + _Private::LaggedFibonacci4Helper< t_DiscreteType >::STATE_SIZE;
	BOOST_ASSERT( iiMiddle + _Private::LaggedFibonacci4Helper< t_FloatingType >::STATE_SIZE == iiEnd );

	m_discreteHelper.Seed( iiBegin,  iiMiddle );
	m_floatingHelper.Seed( iiMiddle, iiEnd    );
}


template< typename t_DiscreteType, typename t_FloatingType >
t_DiscreteType const LaggedFibonacci4< t_DiscreteType, t_FloatingType >::SampleDiscreteUniform() {

	return m_discreteHelper();
}


template< typename t_DiscreteType, typename t_FloatingType >
t_FloatingType const LaggedFibonacci4< t_DiscreteType, t_FloatingType >::SampleStandardUniform() {

	return m_floatingHelper();
}


template< typename t_DiscreteType, typename t_FloatingType >
template< typename t_Archive >
void LaggedFibonacci4< t_DiscreteType, t_FloatingType >::serialize( t_Archive& archive, unsigned int const ) {

	archive & m_discreteHelper;
	archive & m_floatingHelper;
}




}    // namespace Generator


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_GENERATOR_LAGGED_FIBONACCI_4_HPP__ */
