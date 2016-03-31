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
	\file random_distribution_private_ziggurat.hpp
	\brief Random::Distribution Ziggurat helper functions
*/




#ifndef __RANDOM_DISTRIBUTION_PRIVATE_ZIGGURAT_HPP__
#define __RANDOM_DISTRIBUTION_PRIVATE_ZIGGURAT_HPP__

#ifdef __cplusplus




#include <boost/type_traits.hpp>
#include <boost/cstdint.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <limits>




namespace Random {


namespace Distribution {


namespace _Private {




//============================================================================
//    Ziggurat helper functions
//============================================================================


inline bool const BisectionHelper(
	unsigned int const levels,
	long double const maximum_yy,
	long double const total_area,
	long double const radius,
	long double const ( density        )( long double const ),
	long double const ( inverseDensity )( long double const ),
	long double const ( tailArea       )( long double const )
)
{
	long double xx = radius;
	long double yy = density( xx );

	long double const level_area = tailArea( radius ) + xx * yy;
	if ( ( level_area * levels ) < total_area )
		return true;

	for ( unsigned int ii = 1; ii < levels; ++ii ) {

		yy += level_area / xx;
		if ( yy > maximum_yy )
			return false;
		xx = inverseDensity( yy );
	}

	return true;
}


inline long double const Bisection(
	unsigned int const levels,
	long double const ( density        )( long double const ),
	long double const ( inverseDensity )( long double const ),
	long double const ( tailArea       )( long double const )
)
{
	long double const maximum_yy = density( 0.0l );
	long double const total_area = tailArea( 0.0l );

	long double lower = 1.0l;
	long double upper = lower;
	if ( BisectionHelper( levels, maximum_yy, total_area, lower, density, inverseDensity, tailArea ) ) {

		do {

			upper = lower;
			lower *= 0.5;

		} while ( BisectionHelper( levels, maximum_yy, total_area, lower, density, inverseDensity, tailArea ) );
	}
	else {

		do {

			lower = upper;
			upper *= 2.0;

		} while ( ! BisectionHelper( levels, maximum_yy, total_area, upper, density, inverseDensity, tailArea ) );
	}

	long double difference = upper - lower;
	BOOST_ASSERT( difference >= 0.0 );
	for ( ; ; ) {

		long double const middle = lower + 0.5 * difference;

		if ( BisectionHelper( levels, maximum_yy, total_area, middle, density, inverseDensity, tailArea ) )
			upper = middle;
		else
			lower = middle;

		long double const old_difference = difference;
		difference = upper - lower;
		BOOST_ASSERT( difference >= 0.0 );
		if ( difference == old_difference )
			break;
	}

	return upper;
}


template< typename t_Type >
inline void InitializeZiggurat(
	t_Type* px,
	t_Type* py,
	unsigned int const levels,
	long double const ( density        )( long double const ),
	long double const ( inverseDensity )( long double const ),
	long double const ( tailArea       )( long double const )
)
{
	BOOST_ASSERT( levels >= 2 );

	long double const radius = Bisection( levels, density, inverseDensity, tailArea );

	long double xx = radius;
	long double yy = density( xx );

	long double const level_area = tailArea( radius ) + xx * yy;
	BOOST_ASSERT( ( level_area * levels ) >= tailArea( 0.0l ) );

	/*
		when we sample from the bottom rectangle, we *never* reject. However,
		with some probability (determined here), we will sample from the tail
	*/
	px[ levels ] = level_area / yy;
	py[ levels ] = 0.0l;

	px[ levels - 1 ] = xx;
	py[ levels - 1 ] = yy;
	for ( unsigned int ii = levels - 2; ii > 0; --ii ) {

		yy += level_area / xx;
		xx = inverseDensity( yy );

		px[ ii ] = xx;
		py[ ii ] = yy;
	}
	px[ 0 ] = 0.0l;
	py[ 0 ] = density( 0.0l );
}




}    // namespace _Private


}    // namespace Distribution


}    // namespace Random




#endif    /* __cplusplus */

#endif    /* __RANDOM_DISTRIBUTION_PRIVATE_ZIGGURAT_HPP__ */
