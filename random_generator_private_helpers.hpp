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
	\file random_generator_private_helpers.hpp
	\brief Random generator private helper functions
*/




#ifndef __RANDOM_GENERATOR_PRIVATE_HELPERS_HPP__
#define __RANDOM_GENERATOR_PRIVATE_HELPERS_HPP__

#ifdef __cplusplus




#include <boost/type_traits.hpp>
#include <boost/cstdint.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <limits>
#include <cmath>




namespace Random {


namespace Generator {


namespace _Private {




//============================================================================
//    NaNHelper class
//============================================================================


template< typename t_Type, bool t_IsFloating = boost::is_floating_point< t_Type >::value >
struct NaNHelper;


template< typename t_Type >
struct NaNHelper< t_Type, false > {

	BOOST_STATIC_ASSERT(
		boost::is_integral< t_Type >::value &&
		boost::is_unsigned< t_Type >::value
	);

	static inline t_Type const NaN() {

		return 0;
	}
};


template< typename t_Type >
struct NaNHelper< t_Type, true > {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );

	static inline t_Type const NaN() {

		return std::numeric_limits< t_Type >::quiet_NaN();
	}
};




//============================================================================
//    SampleHelper class
//============================================================================


template< typename t_Type, bool t_IsFloating = boost::is_floating_point< t_Type >::value >
struct SampleHelper;


template< typename t_Type >
struct SampleHelper< t_Type, false > {

	BOOST_STATIC_ASSERT(
		boost::is_integral< t_Type >::value &&
		boost::is_unsigned< t_Type >::value
	);

	template< typename t_Generator >
	static inline t_Type const Sample( t_Generator& generator ) {

		return generator.SampleDiscreteUniform();
	}
};


template< typename t_Type >
struct SampleHelper< t_Type, true > {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );

	template< typename t_Generator >
	static inline t_Type const Sample( t_Generator& generator ) {

		return generator.SampleStandardUniform();
	}
};




//============================================================================
//    CastHelper class and Cast function
//============================================================================


template<
	typename t_Type,
	typename t_OtherType,
	bool t_IsFloating = boost::is_floating_point< t_Type >::value,
	bool t_OtherIsFloating = boost::is_floating_point< t_OtherType >::value,
	bool t_LessThan = ( sizeof( t_Type ) < sizeof( t_OtherType ) ),
	bool t_GreaterThan = ( sizeof( t_Type ) > sizeof( t_OtherType ) )
>
struct CastHelper;


template< typename t_Type, typename t_OtherType >
struct CastHelper< t_Type, t_OtherType, false, false, false, false > {

	BOOST_STATIC_ASSERT(
		boost::is_integral< t_Type >::value &&
		boost::is_unsigned< t_Type >::value
	);
	BOOST_STATIC_ASSERT(
		boost::is_integral< t_OtherType >::value &&
		boost::is_unsigned< t_OtherType >::value
	);
	BOOST_STATIC_ASSERT( sizeof( t_Type ) == sizeof( t_OtherType ) );

	static inline t_Type const Cast( t_OtherType const& other ) {

		return other;
	}
};


template< typename t_Type, typename t_OtherType >
struct CastHelper< t_Type, t_OtherType, false, false, true, false > {

	BOOST_STATIC_ASSERT(
		boost::is_integral< t_Type >::value &&
		boost::is_unsigned< t_Type >::value
	);
	BOOST_STATIC_ASSERT(
		boost::is_integral< t_OtherType >::value &&
		boost::is_unsigned< t_OtherType >::value
	);
	BOOST_STATIC_ASSERT( sizeof( t_Type ) < sizeof( t_OtherType ) );

	static inline t_Type const Cast( t_OtherType const& other ) {

		return( other >> ( sizeof( t_OtherType ) - sizeof( t_Type ) ) * 8 );
	}
};


template< typename t_Type, typename t_OtherType >
struct CastHelper< t_Type, t_OtherType, false, false, false, true > {

	BOOST_STATIC_ASSERT(
		boost::is_integral< t_Type >::value &&
		boost::is_unsigned< t_Type >::value
	);
	BOOST_STATIC_ASSERT(
		boost::is_integral< t_OtherType >::value &&
		boost::is_unsigned< t_OtherType >::value
	);
	BOOST_STATIC_ASSERT( sizeof( t_Type ) > sizeof( t_OtherType ) );

	static inline t_Type const Cast( t_OtherType const& other ) {

		return( static_cast< t_Type >( other ) << ( sizeof( t_Type ) - sizeof( t_OtherType ) ) * 8 );
	}
};


template< typename t_Type, typename t_OtherType, bool t_LessThan, bool t_GreaterThan >
struct CastHelper< t_Type, t_OtherType, true, false, t_LessThan, t_GreaterThan > {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );
	BOOST_STATIC_ASSERT(
		boost::is_integral< t_OtherType >::value &&
		boost::is_unsigned< t_OtherType >::value
	);

	static inline t_Type const Cast( t_OtherType const& other ) {

		t_Type result = static_cast< t_Type >( other ) / static_cast< t_Type >( std::pow( 2.0, sizeof( t_OtherType ) * 8.0 ) );    // **NOTE: this pow had better be removed at compile-time
		if ( result >= 1 )    // if t_OtherType has more precision than t_Type
			--result;
		BOOST_ASSERT( ( result >= 0 ) && ( result < 1 ) );
		return result;
	}
};


template< typename t_Type, typename t_OtherType, bool t_GreaterThan >
struct CastHelper< t_Type, t_OtherType, true, true, true, t_GreaterThan > {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type      >::value );
	BOOST_STATIC_ASSERT( boost::is_floating_point< t_OtherType >::value );
	BOOST_STATIC_ASSERT( sizeof( t_Type ) < sizeof( t_OtherType ) );

	static inline t_Type const Cast( t_OtherType const& other ) {

		BOOST_ASSERT( ( other >= 0 ) && ( other < 1 ) );
		t_Type result = other;
		if ( result >= 1 )    // if t_OtherType has more precision than t_Type
			--result;
		BOOST_ASSERT( ( result >= 0 ) && ( result < 1 ) );
		return result;
	}
};


template< typename t_Type, typename t_OtherType, bool t_GreaterThan >
struct CastHelper< t_Type, t_OtherType, true, true, false, t_GreaterThan > {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type      >::value );
	BOOST_STATIC_ASSERT( boost::is_floating_point< t_OtherType >::value );
	BOOST_STATIC_ASSERT( sizeof( t_Type ) >= sizeof( t_OtherType ) );

	static inline t_Type const Cast( t_OtherType const& other ) {

		BOOST_ASSERT( ( other >= 0 ) && ( other < 1 ) );
		t_Type result = other;
		BOOST_ASSERT( ( result >= 0 ) && ( result < 1 ) );
		return result;
	}
};


template< typename t_Type, typename t_OtherType, bool t_LessThan, bool t_GreaterThan >
struct CastHelper< t_Type, t_OtherType, false, true, t_LessThan, t_GreaterThan > {

	BOOST_STATIC_ASSERT(
		boost::is_integral< t_Type >::value &&
		boost::is_unsigned< t_Type >::value
	);
	BOOST_STATIC_ASSERT( boost::is_floating_point< t_OtherType >::value );

	static inline t_Type const Cast( t_OtherType const& other ) {

		BOOST_ASSERT( ( other >= 0 ) && ( other < 1 ) );
		t_Type const result = other * static_cast< t_OtherType >( std::pow( 2.0, sizeof( t_Type ) * 8.0 ) );    // **NOTE: this pow had better be removed at compile-time
		return result;
	}
};


template< typename t_Type, typename t_OtherType >
inline t_Type const Cast( t_OtherType const& other ) {

	return CastHelper< t_Type, t_OtherType >::Cast( other );
}




//============================================================================
//    ModulusHelper class
//============================================================================


template< typename t_Type, bool t_IsFloating = boost::is_floating_point< t_Type >::value >
struct ModulusHelper;


template< typename t_Type >
struct ModulusHelper< t_Type, false > {

	BOOST_STATIC_ASSERT(
		boost::is_integral< t_Type >::value &&
		boost::is_unsigned< t_Type >::value
	);

	static inline t_Type const& Modulus( t_Type const& value ) {

		return value;
	}
};


template< typename t_Type >
struct ModulusHelper< t_Type, true > {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );

	static inline t_Type const Modulus( t_Type value ) {

		BOOST_ASSERT( ( value >= 0 ) && ( value < 4 ) );
		if ( value >= 2 )
			value -= 2;
		if ( value >= 1 )
			--value;
		BOOST_ASSERT( ( value >= 0 ) && ( value < 1 ) );

		return value;
	}
};




//============================================================================
//    IntegerHelper class
//============================================================================


template< unsigned int t_Integer >
struct IntegerHelper;

template<> struct IntegerHelper< 1 > { typedef uint8_t  Type; };
template<> struct IntegerHelper< 2 > { typedef uint16_t Type; };
template<> struct IntegerHelper< 4 > { typedef uint32_t Type; };
#ifndef BOOST_NO_INT64_T
template<> struct IntegerHelper< 8 > { typedef uint64_t Type; };
#endif    // BOOST_NO_INT64_T




}    // namespace _Private


}    // namespace Generator


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_GENERATOR_PRIVATE_HELPERS_HPP__ */
