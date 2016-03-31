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
	\file helpers.h
	\brief Helper macros, functions and classes
*/




#ifndef __HELPERS_H__
#define __HELPERS_H__




//============================================================================
//    ARRAYLENGTH macro
//============================================================================


#define ARRAYLENGTH( array )  \
	( sizeof( array ) / sizeof( array[ 0 ] ) )




//============================================================================
//    LIKELY and UNLIKELY macros
//============================================================================


#if defined( __GNUC__ ) && ( __GNUC__ >= 3 )

#define LIKELY( boolean ) __builtin_expect( ( boolean ), 1 )
#define UNLIKELY( boolean ) __builtin_expect( ( boolean ), 0 )

#else    /* defined( __GNUC__ ) && ( __GNUC__ >= 3 ) */

#define LIKELY( boolean ) ( boolean )
#define UNLIKELY( boolean ) ( boolean )

#endif    /* defined( __GNUC__ ) && ( __GNUC__ >= 3 ) */




#ifdef __cplusplus




#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <stdint.h>




//============================================================================
//    PAREN_TYPE macro
//============================================================================


namespace _Private {


template< typename t_Type >
struct ParenHelper;

template< typename t_Type >
struct ParenHelper< void ( t_Type ) > {

	typedef t_Type Type;
};


}    // namespace _Private


/*
	http://lists.boost.org/Archives/boost/2004/08/70107.php
		Post by Paul Mensonides

	Wraps parenthesis around a type, for use in macro parameters which will
	themselves be used as parameters to nested macros. This is just to help you
	get through multiple levels of macro definitions--you'll still need to use
	BOOST_PP_COMMA() to get through the *first* layer.
*/
#define PAREN_TYPE( tt ) _Private::ParenHelper< void ( tt ) >::Type




//============================================================================
//    Square helper function
//============================================================================


template< typename t_Type >
inline t_Type Square( t_Type const& value ) {

	return( value * value );
}




//============================================================================
//    Cube helper function
//============================================================================


template< typename t_Type >
inline t_Type Cube( t_Type const& value ) {

	return( value * value * value );
}




//============================================================================
//    CountBits helper functions
//============================================================================


inline uint8_t const CountBits( uint8_t number ) {

	BOOST_STATIC_ASSERT( sizeof( number ) == 1 );
	number = ( number & 0x55 ) + ( ( number >> 1 ) & 0x55 );
	number = ( number & 0x33 ) + ( ( number >> 2 ) & 0x33 );
	number = ( number & 0x0f ) +   ( number >> 4 );
	return number;
}


inline uint16_t const CountBits( uint16_t number ) {

	BOOST_STATIC_ASSERT( sizeof( number ) == 2 );
	number = ( number & 0x5555 ) + ( ( number >> 1 ) & 0x5555 );
	number = ( number & 0x3333 ) + ( ( number >> 2 ) & 0x3333 );
	number = ( number & 0x0f0f ) + ( ( number >> 4 ) & 0x0f0f );
	number = ( number & 0x00ff ) +   ( number >> 8 );
	return number;
}


inline uint32_t const CountBits( uint32_t number ) {

	BOOST_STATIC_ASSERT( sizeof( number ) == 4 );
	number = ( number & 0x55555555 ) + ( ( number >>  1 ) & 0x55555555 );
	number = ( number & 0x33333333 ) + ( ( number >>  2 ) & 0x33333333 );
	number = ( number & 0x0f0f0f0f ) + ( ( number >>  4 ) & 0x0f0f0f0f );
	number = ( number & 0x00ff00ff ) + ( ( number >>  8 ) & 0x00ff00ff );
	number = ( number & 0x0000ffff ) +   ( number >> 16 );
	return number;
}


inline uint64_t const CountBits( uint64_t number ) {

	BOOST_STATIC_ASSERT( sizeof( number ) == 8 );
	number = ( number & 0x5555555555555555ull ) + ( ( number >>  1 ) & 0x5555555555555555ull );
	number = ( number & 0x3333333333333333ull ) + ( ( number >>  2 ) & 0x3333333333333333ull );
	number = ( number & 0x0f0f0f0f0f0f0f0full ) + ( ( number >>  4 ) & 0x0f0f0f0f0f0f0f0full );
	number = ( number & 0x00ff00ff00ff00ffull ) + ( ( number >>  8 ) & 0x00ff00ff00ff00ffull );
	number = ( number & 0x0000ffff0000ffffull ) + ( ( number >> 16 ) & 0x0000ffff0000ffffull );
	number = ( number & 0x00000000ffffffffull ) +   ( number >> 32 );
	return number;
}




//============================================================================
//    CountLowBits helper functions
//============================================================================


inline unsigned int const CountLowBits( uint8_t number ) {

	BOOST_STATIC_ASSERT( sizeof( number ) == 1 );

	unsigned int bits = 0;

	if ( number & 0xf0 ) {

		bits += 4;
		number >>= 4;
	}

	if ( number & 0x0c ) {

		bits += 2;
		number >>= 2;
	}

	if ( number & 0x02 )
		bits += 2;
	else if ( number & 0x01 )
		++bits;

	return bits;
}


inline unsigned int const CountLowBits( uint16_t number ) {

	BOOST_STATIC_ASSERT( sizeof( number ) == 2 );

	unsigned int bits = 0;

	if ( number & 0xff00 ) {

		bits += 8;
		number >>= 8;
	}

	if ( number & 0x00f0 ) {

		bits += 4;
		number >>= 4;
	}

	if ( number & 0x000c ) {

		bits += 2;
		number >>= 2;
	}

	if ( number & 0x0002 )
		bits += 2;
	else if ( number & 0x0001 )
		++bits;

	return bits;
}


inline unsigned int const CountLowBits( uint32_t number ) {

	BOOST_STATIC_ASSERT( sizeof( number ) == 4 );

	unsigned int bits = 0;

	if ( number & 0xffff0000 ) {

		bits += 16;
		number >>= 16;
	}

	if ( number & 0x0000ff00 ) {

		bits += 8;
		number >>= 8;
	}

	if ( number & 0x000000f0 ) {

		bits += 4;
		number >>= 4;
	}

	if ( number & 0x0000000c ) {

		bits += 2;
		number >>= 2;
	}

	if ( number & 0x00000002 )
		bits += 2;
	else if ( number & 0x00000001 )
		++bits;

	return bits;
}


inline unsigned int const CountLowBits( uint64_t number ) {

	BOOST_STATIC_ASSERT( sizeof( number ) == 8 );

	unsigned int bits = 0;

	if ( number & 0xffffffff00000000ull ) {

		bits += 32;
		number >>= 32;
	}

	if ( number & 0x00000000ffff0000ull ) {

		bits += 16;
		number >>= 16;
	}

	if ( number & 0x000000000000ff00ull ) {

		bits += 8;
		number >>= 8;
	}

	if ( number & 0x00000000000000f0ull ) {

		bits += 4;
		number >>= 4;
	}

	if ( number & 0x000000000000000cull ) {

		bits += 2;
		number >>= 2;
	}

	if ( number & 0x0000000000000002ull )
		bits += 2;
	else if ( number & 0x0000000000000001ull )
		++bits;

	return bits;
}




#endif    /* __cplusplus */




#endif    /* __HELPERS_H__ */
