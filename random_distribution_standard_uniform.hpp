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
	\file random_distribution_standard_uniform.hpp
	\brief StandardUniform implementation
*/




#ifndef __RANDOM_DISTRIBUTION_STANDARD_UNIFORM_HPP__
#define __RANDOM_DISTRIBUTION_STANDARD_UNIFORM_HPP__

#ifdef __cplusplus




#include "helpers.h"

#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <limits>




namespace Random {


namespace Distribution {




//============================================================================
//    StandardUniform class
//============================================================================


template< typename t_Type = double >
struct StandardUniform {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );

	typedef t_Type Type;


	inline StandardUniform();


	inline double const PDF( t_Type const value ) const;
	inline double const CDF( t_Type const value ) const;

	template< typename t_Generator >
	inline t_Type const Sample( t_Generator& generator );
};




//============================================================================
//    StandardUniform inline methods
//============================================================================


template< typename t_Type >
StandardUniform< t_Type >::StandardUniform() {
}


template< typename t_Type >
double const StandardUniform< t_Type >::PDF( t_Type const value ) const {

	return( ( ( value >= 0 ) && ( value < 1 ) ) ? 1 : 0 );
}


template< typename t_Type >
double const StandardUniform< t_Type >::CDF( t_Type const value ) const {

	return( ( value <= 0 ) ? 0 : ( ( value >= 1 ) ? 1 : value ) );
}


template< typename t_Type >
template< typename t_Generator >
t_Type const StandardUniform< t_Type >::Sample( t_Generator& generator ) {

	return generator.SampleStandardUniform();
}




}    // namespace Distribution


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_DISTRIBUTION_STANDARD_UNIFORM_HPP__ */
