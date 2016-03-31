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
	\file random_generator_private_linear_congruential_helper.hpp
	\brief LinearCongruentialHelper class
*/




#ifndef __RANDOM_GENERATOR_PRIVATE_LINEAR_CONGRUENTIAL_HELPER_HPP__
#define __RANDOM_GENERATOR_PRIVATE_LINEAR_CONGRUENTIAL_HELPER_HPP__

#ifdef __cplusplus




#include "random_generator_private_helpers.hpp"
#include "helpers.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/type_traits.hpp>
#include <boost/cstdint.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <limits>
#include <sstream>
#include <stdexcept>

#include <errno.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>




namespace Random {


namespace Generator {


namespace _Private {




//============================================================================
//    LinearCongruentialHelper class
//============================================================================


/*
	Pierre L'Ecuyer, Richard Simard. "TestU01: A C Library for Empirical
	Testing of Random Number Helpers"
		- Used for uint32_t specialization
	http://nuclear.llnl.gov/CNP/rng/rngman/node4.html
		- Used for uint64_t specialization
*/
template< typename t_Type >
struct LinearCongruentialHelper {

	BOOST_STATIC_ASSERT(
		boost::is_integral< t_Type >::value &&
		boost::is_unsigned< t_Type >::value
	);

	typedef t_Type Type;

	enum { STATE_SIZE = 1 };


	inline LinearCongruentialHelper();
	inline explicit LinearCongruentialHelper( LinearCongruentialHelper const& other );


	inline void Seed();    // from /dev/urandom

	inline void Seed( t_Type const& seed );

	template< typename t_Generator >
	inline void Seed( t_Generator& generator );

	template< typename t_Iterator >
	inline void Seed( t_Iterator const& iiBegin, t_Iterator const& iiEnd );


	inline t_Type const operator()();


private:

	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	t_Type m_state;


	friend class boost::serialization::access;
};




//============================================================================
//    LinearCongruentialHelper inline methods
//============================================================================


template< typename t_Type >
LinearCongruentialHelper< t_Type >::LinearCongruentialHelper() : m_state( 0 ) {
}


template< typename t_Type >
LinearCongruentialHelper< t_Type >::LinearCongruentialHelper( LinearCongruentialHelper const& other ) : m_state( other.m_state ) {
}


template< typename t_Type >
void LinearCongruentialHelper< t_Type >::Seed() {

	int handle = -1;
	try {

		handle = open( "/dev/urandom", O_RDONLY, S_IREAD );
		if ( handle < 0 ) {

			std::stringstream stream;
			stream << "Unable to open \"/dev/urandom\" for reading: " << strerror( errno );
			throw std::runtime_error( stream.str() );
		}

		int const bytes = read( handle, &m_state, sizeof( m_state ) );
		if ( bytes < 0 ) {

			std::stringstream stream;
			stream << "Unable to read from \"/dev/urandom\": " << strerror( errno );
			throw std::runtime_error( stream.str() );
		}
		else if ( bytes != sizeof( m_state ) ) {

			std::stringstream stream;
			stream << "Read " << bytes << " of " << sizeof( m_state ) << " bytes from \"/dev/urandom\"";
			throw std::runtime_error( stream.str() );
		}

		close( handle );
		handle = -1;
	}
	catch( ... ) {

		if ( handle >= 0 )
			close( handle );
		throw;
	}
	BOOST_ASSERT( handle == -1 );
}


template< typename t_Type >
void LinearCongruentialHelper< t_Type >::Seed( t_Type const& seed ) {

	m_state = seed;
}


template< typename t_Type >
template< typename t_Generator >
void LinearCongruentialHelper< t_Type >::Seed( t_Generator& generator ) {

	m_state = Cast< t_Type >( generator.SampleDiscreteUniform() );
}


template< typename t_Type >
template< typename t_Iterator >
void LinearCongruentialHelper< t_Type >::Seed( t_Iterator const& iiBegin, t_Iterator const& iiEnd ) {

	BOOST_ASSERT( iiBegin + 1 == iiEnd );
	m_state = Cast< t_Type >( *iiBegin );
}


template<>
inline uint32_t const LinearCongruentialHelper< uint32_t >::operator()() {

	m_state = m_state * 69069u + 1u;
	return m_state;
}


template<>
inline uint64_t const LinearCongruentialHelper< uint64_t >::operator()() {

	m_state = m_state * 2862933555777941757ull + 3037000493ull;
	return m_state;
}


template< typename t_Type >
template< typename t_Archive >
void LinearCongruentialHelper< t_Type >::serialize( t_Archive& archive, unsigned int const ) {

	archive & m_state;
}




}    // namespace _Private


}    // namespace Generator


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_GENERATOR_PRIVATE_GENERATOR_PRIVATE_LINEAR_CONGRUENTIAL_HELPER_HPP__ */
