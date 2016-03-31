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
	\file random_generator_private_lagged_fibonacci_4_helper.hpp
	\brief LaggedFibonacci4Helper class
*/




#ifndef __RANDOM_GENERATOR_PRIVATE_LAGGED_FIBONACCI_4_HELPER_HPP__
#define __RANDOM_GENERATOR_PRIVATE_LAGGED_FIBONACCI_4_HELPER_HPP__

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
//    LaggedFibonacci4Helper class
//============================================================================


/*
	Generates floating point numbers uniformly in the range [0,1), or unsigned
	integers uniformly on their range.

	<b>George Marsaglia. "Random numbers for C: The END?".
	<a href="http://groups.google.com/group/sci.crypt/browse_thread/thread/ca8682a4658a124d/">Posted
	to sci.crypt mailing list</a>. 1999.</b>
*/
template< typename t_Type >
struct LaggedFibonacci4Helper {

	BOOST_STATIC_ASSERT(
		boost::is_floating_point< t_Type >::value ||
		(
			boost::is_integral< t_Type >::value &&
			boost::is_unsigned< t_Type >::value
		)
	);

	typedef t_Type Type;


	enum { STATE_SIZE = 256 };


	inline LaggedFibonacci4Helper();
	inline explicit LaggedFibonacci4Helper( LaggedFibonacci4Helper const& other );


	inline void Seed();    // from /dev/urandom

	template< typename t_Generator >
	inline void Seed( t_Generator& generator );

	template< typename t_Iterator >
	inline void Seed( t_Iterator const& iiBegin, t_Iterator const& iiEnd );


	inline t_Type const operator()();


private:

	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	unsigned int m_index;

	t_Type m_buffer[ STATE_SIZE ];


	friend class boost::serialization::access;
};




//============================================================================
//    LaggedFibonacci4Helper inline methods
//============================================================================


template< typename t_Type >
LaggedFibonacci4Helper< t_Type >::LaggedFibonacci4Helper() : m_index( std::numeric_limits< unsigned int >::max() ) {

#ifndef NDEBUG
	std::fill( m_buffer, m_buffer + ARRAYLENGTH( m_buffer ), NaNHelper< t_Type >::NaN() );
#endif    // NDEBUG
}


template< typename t_Type >
LaggedFibonacci4Helper< t_Type >::LaggedFibonacci4Helper( LaggedFibonacci4Helper const& other ) : m_index( other.m_index ) {

	std::copy( other.m_buffer, other.m_buffer + ARRAYLENGTH( m_buffer ), m_buffer );
}


template< typename t_Type >
void LaggedFibonacci4Helper< t_Type >::Seed() {

	typedef typename IntegerHelper< sizeof( t_Type ) >::Type BufferType;
	BufferType buffer[ ARRAYLENGTH( m_buffer ) ];

	int handle = -1;
	try {

		handle = open( "/dev/urandom", O_RDONLY, S_IREAD );
		if ( handle < 0 ) {

			std::stringstream stream;
			stream << "Unable to open \"/dev/urandom\" for reading: " << strerror( errno );
			throw std::runtime_error( stream.str() );
		}

		int const bytes = read( handle, &buffer[ 0 ], sizeof( buffer ) );
		if ( bytes < 0 ) {

			std::stringstream stream;
			stream << "Unable to read from \"/dev/urandom\": " << strerror( errno );
			throw std::runtime_error( stream.str() );
		}
		else if ( bytes != sizeof( buffer ) ) {

			std::stringstream stream;
			stream << "Read " << bytes << " of " << sizeof( buffer ) << " bytes from \"/dev/urandom\"";
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

	Seed( &buffer[ 0 ], &buffer[ 0 ] + ARRAYLENGTH( buffer ) );
}


template< typename t_Type >
template< typename t_Generator >
void LaggedFibonacci4Helper< t_Type >::Seed( t_Generator& generator ) {

	t_Type* ii    = m_buffer;
	t_Type* iiEnd = ii + ARRAYLENGTH( m_buffer );
	for ( ; ii != iiEnd; ++ii )
		*ii = SampleHelper< t_Type >::Sample( generator );

	m_index = 0;
}


template< typename t_Type >
template< typename t_Iterator >
void LaggedFibonacci4Helper< t_Type >::Seed( t_Iterator const& iiBegin, t_Iterator const& iiEnd ) {

	BOOST_ASSERT( ( iiEnd - iiBegin ) == ARRAYLENGTH( m_buffer ) );

	t_Type* jj = m_buffer;
	for ( t_Iterator ii = iiBegin; ii != iiEnd; ++ii )
		*( jj++ ) = Cast< t_Type >( *ii );
	BOOST_ASSERT( jj == m_buffer + ARRAYLENGTH( m_buffer ) );

	m_index = 0;
}


template< typename t_Type >
t_Type const LaggedFibonacci4Helper< t_Type >::operator()() {

	BOOST_ASSERT( m_index <= ARRAYLENGTH( m_buffer ) );
	if ( UNLIKELY( m_index == ARRAYLENGTH( m_buffer ) ) ) {

		// m_buffer length must be a power of 2
		BOOST_ASSERT( ( ARRAYLENGTH( m_buffer ) & ( ARRAYLENGTH( m_buffer ) - 1 ) ) == 0 );

		for ( unsigned int ii = 0; ii < ARRAYLENGTH( m_buffer ); ++ii ) {

			m_buffer[ ii ] = ModulusHelper< t_Type >::Modulus(
				m_buffer[ ii ] +
				m_buffer[ ( ii +  58 ) & ( ARRAYLENGTH( m_buffer ) - 1 ) ] +
				m_buffer[ ( ii + 119 ) & ( ARRAYLENGTH( m_buffer ) - 1 ) ] +
				m_buffer[ ( ii + 178 ) & ( ARRAYLENGTH( m_buffer ) - 1 ) ]
			);
		}

		m_index = 0;
	}

	return m_buffer[ m_index++ ];
}


template< typename t_Type >
template< typename t_Archive >
void LaggedFibonacci4Helper< t_Type >::serialize( t_Archive& archive, unsigned int const ) {

	archive & m_index;
	archive & m_buffer;
}




}    // namespace _Private


}    // namespace Generator


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_GENERATOR_PRIVATE_LAGGED_FIBONACCI_4_HELPER_HPP__ */
