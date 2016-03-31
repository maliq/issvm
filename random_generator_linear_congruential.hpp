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
	\file random_generator_linear_congruential.hpp
	\brief LinearCongruential class
*/




#ifndef __RANDOM_GENERATOR_LINEAR_CONGRUENTIAL_HPP__
#define __RANDOM_GENERATOR_LINEAR_CONGRUENTIAL_HPP__

#ifdef __cplusplus




#include "random_generator_private_helpers.hpp"
#include "random_generator_private_linear_congruential_helper.hpp"
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
//    LinearCongruential class
//============================================================================


template<
	typename t_DiscreteType = unsigned int,
	typename t_FloatingType = double
>
struct LinearCongruential {

	BOOST_STATIC_ASSERT(
		(
			boost::is_integral< t_DiscreteType >::value &&
			boost::is_unsigned< t_DiscreteType >::value
		) &&
		boost::is_floating_point< t_FloatingType >::value
	);

	typedef t_DiscreteType DiscreteType;
	typedef t_FloatingType FloatingType;

	enum { STATE_SIZE = _Private::LinearCongruentialHelper< t_DiscreteType >::STATE_SIZE };


	inline LinearCongruential();
	inline explicit LinearCongruential( LinearCongruential const& other );


	inline void Seed();    // from /dev/urandom

	inline void Seed( t_DiscreteType const& seed );

	template< typename t_Generator >
	inline void Seed( t_Generator& generator );

	template< typename t_Iterator >
	inline void Seed( t_Iterator const& iiBegin, t_Iterator const& iiEnd );


	inline t_DiscreteType const SampleDiscreteUniform();
	inline t_FloatingType const SampleStandardUniform();


private:

	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	_Private::LinearCongruentialHelper< t_DiscreteType > m_helper;


	friend class boost::serialization::access;
};




//============================================================================
//    LinearCongruential inline methods
//============================================================================


template< typename t_DiscreteType, typename t_FloatingType >
LinearCongruential< t_DiscreteType, t_FloatingType >::LinearCongruential() {
}


template< typename t_DiscreteType, typename t_FloatingType >
LinearCongruential< t_DiscreteType, t_FloatingType >::LinearCongruential( LinearCongruential const& other ) :
	m_helper( other.m_helper )
{
}


template< typename t_DiscreteType, typename t_FloatingType >
void LinearCongruential< t_DiscreteType, t_FloatingType >::Seed() {

	m_helper.Seed();
}


template< typename t_DiscreteType, typename t_FloatingType >
void LinearCongruential< t_DiscreteType, t_FloatingType >::Seed( t_DiscreteType const& seed ) {

	m_helper.Seed( seed );
}


template< typename t_DiscreteType, typename t_FloatingType >
template< typename t_Generator >
void LinearCongruential< t_DiscreteType, t_FloatingType >::Seed( t_Generator& generator ) {

	m_helper.Seed( generator );
}


template< typename t_DiscreteType, typename t_FloatingType >
template< typename t_Iterator >
void LinearCongruential< t_DiscreteType, t_FloatingType >::Seed( t_Iterator const& iiBegin, t_Iterator const& iiEnd ) {

	m_helper.Seed( iiBegin, iiEnd );
}


template< typename t_DiscreteType, typename t_FloatingType >
t_DiscreteType const LinearCongruential< t_DiscreteType, t_FloatingType >::SampleDiscreteUniform() {

	return m_helper();
}


template< typename t_DiscreteType, typename t_FloatingType >
t_FloatingType const LinearCongruential< t_DiscreteType, t_FloatingType >::SampleStandardUniform() {

	t_FloatingType static const scale = std::pow( 2.0, -static_cast< double >( sizeof( t_DiscreteType ) ) * 8.0 );
	t_FloatingType const result = m_helper() * scale;
	BOOST_ASSERT( ( result >= 0 ) && ( result < 1 ) );
	return result;
}


template< typename t_DiscreteType, typename t_FloatingType >
template< typename t_Archive >
void LinearCongruential< t_DiscreteType, t_FloatingType >::serialize( t_Archive& archive, unsigned int const ) {

	archive & m_helper;
}




}    // namespace Generator


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_GENERATOR_LINEAR_CONGRUENTIAL_HPP__ */
