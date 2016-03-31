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
	\file svm_optimizer_classification_private_find_water_level.cpp
	\brief SVM::Optimizer::Classification::Biased::_Private::FindWaterLevel and SVM::Optimizer::Classification::Unbiased::_Private::FindWaterLevel definitions
*/




#include "helpers.h"

#include <boost/assert.hpp>

#include <algorithm>
#include <limits>
#include <utility>




//============================================================================
//    Partition helper function
//============================================================================


static double* const Partition( double* const begin, double* const end ) {

	BOOST_ASSERT( end > begin );

	double* iiLow  = begin;
	double* iiHigh = end;
	--iiHigh;

	double* ii = begin;
	double const pivot = *( ii++ );    // first element is the pivot

	while ( ii <= iiHigh ) {

		if ( *ii < pivot )
			*( iiLow++ ) = *( ii++ );
		else if ( *ii > pivot )
			std::swap( *( iiHigh-- ), *ii );
		else
			++ii;
	}
	BOOST_ASSERT( iiLow <= iiHigh );

	for ( ii = iiLow; ii <= iiHigh; ++ii )
		*ii = pivot;

	// try to split as nearly in half as possible
	double* result = ( begin + ( end - begin ) / 2 );
	if ( result > iiHigh )
		result = iiHigh;
	else if ( result < iiLow )
		result = iiLow;

	return result;
}




namespace SVM {


namespace Optimizer {


namespace Classification {




namespace Unbiased {


namespace _Private {


//============================================================================
//    FindWaterLevel helper function
//============================================================================


double const FindWaterLevel( double* const begin, double* const end, double const total ) {

	BOOST_ASSERT( end > begin );
	BOOST_ASSERT( total >= 0 );

	double kappa = std::numeric_limits< double >::quiet_NaN();
	if ( total == 0 ) {

		kappa = std::numeric_limits< double >::infinity();
		for ( double const* ii = begin; ii != end; ++ii )
			if ( *ii < kappa )
				kappa = *ii;
	}
	else {

		double* pLower = begin;
		double* pUpper = end;

		double lowerMaximum = -std::numeric_limits< double >::infinity();
		double lowerSum     = 0;

		while ( pUpper > pLower ) {

			double middleMaximum = lowerMaximum;
			double middleSum = lowerSum;

			double* pMiddle = Partition( pLower, pUpper ) + 1;
			BOOST_ASSERT( pMiddle >  pLower );
			BOOST_ASSERT( pMiddle <= pUpper );
			for ( double const* ii = pLower; ii < pMiddle; ++ii ) {

				if ( *ii > middleMaximum )
					middleMaximum = *ii;
				middleSum += *ii;
			}

			unsigned int const middleIndex = pMiddle - begin;
			double const middleVolume = middleMaximum * middleIndex - middleSum;

			if ( LIKELY( middleVolume >= total ) )
				pUpper = pMiddle - 1;
			else {

				pLower = pMiddle;
				lowerMaximum = middleMaximum;
				lowerSum = middleSum;
			}
		}
		BOOST_ASSERT( pLower == pUpper );

		unsigned int const lowerIndex = pLower - begin;
		double const lowerVolume = lowerMaximum * lowerIndex - lowerSum;
		BOOST_ASSERT( lowerVolume <= total );

		kappa = ( total - lowerVolume ) / lowerIndex + lowerMaximum;
		BOOST_ASSERT( kappa >= lowerMaximum );
	}

	return kappa;
}


}    // namespace _Private


}    // namespace Unbiased




namespace Biased {


namespace _Private {


//============================================================================
//    FindWaterLevel helper function
//============================================================================


std::pair< double, double > const FindWaterLevel( double* const positiveBegin, double* const positiveEnd, double* const negativeBegin, double* const negativeEnd, double const total ) {

	BOOST_ASSERT( positiveEnd > positiveBegin );
	BOOST_ASSERT( negativeEnd > negativeBegin );
	BOOST_ASSERT( total >= 0 );

	double positiveKappa = std::numeric_limits< double >::quiet_NaN();
	double negativeKappa = std::numeric_limits< double >::quiet_NaN();
	if ( total == 0 ) {

		positiveKappa = std::numeric_limits< double >::infinity();
		for ( double const* ii = positiveBegin; ii != positiveEnd; ++ii )
			if ( *ii < positiveKappa )
				positiveKappa = *ii;

		negativeKappa = std::numeric_limits< double >::infinity();
		for ( double const* ii = negativeBegin; ii != negativeEnd; ++ii )
			if ( *ii < negativeKappa )
				negativeKappa = *ii;
	}
	else {

		double positiveLowerMaximum  = -std::numeric_limits< double >::infinity();
		double positiveMiddleMaximum = -std::numeric_limits< double >::infinity();
		double positiveLowerSum = 0;
		double positiveMiddleSum = 0;

		double* pPositiveLower = positiveBegin;
		double* pPositiveUpper = positiveEnd;
		double* pPositiveMiddle = Partition( pPositiveLower, pPositiveUpper ) + 1;
		for ( double const* ii = pPositiveLower; ii != pPositiveMiddle; ++ii ) {

			if ( *ii > positiveMiddleMaximum )
				positiveMiddleMaximum = *ii;
			positiveMiddleSum += *ii;
		}

		double negativeLowerMaximum  = -std::numeric_limits< double >::infinity();
		double negativeMiddleMaximum = -std::numeric_limits< double >::infinity();
		double negativeLowerSum = 0;
		double negativeMiddleSum = 0;

		double* pNegativeLower = negativeBegin;
		double* pNegativeUpper = negativeEnd;
		double* pNegativeMiddle = Partition( pNegativeLower, pNegativeUpper ) + 1;
		for ( double const* ii = pNegativeLower; ii < pNegativeMiddle; ++ii ) {

			if ( *ii > negativeMiddleMaximum )
				negativeMiddleMaximum = *ii;
			negativeMiddleSum += *ii;
		}

		while ( ( pPositiveUpper > pPositiveLower ) || ( pNegativeUpper > pNegativeLower ) ) {

			int positiveDirection = 0;
			int negativeDirection = 0;

			{	unsigned int const positiveLowerIndex  = pPositiveLower  - positiveBegin;
				unsigned int const positiveMiddleIndex = pPositiveMiddle - positiveBegin;
				unsigned int const positiveUpperIndex  = pPositiveUpper  - positiveBegin;

				unsigned int const negativeLowerIndex  = pNegativeLower  - negativeBegin;
				unsigned int const negativeMiddleIndex = pNegativeMiddle - negativeBegin;
				unsigned int const negativeUpperIndex  = pNegativeUpper  - negativeBegin;

				if ( positiveMiddleIndex > negativeUpperIndex )
					positiveDirection = -1;
				else if ( positiveMiddleIndex < negativeLowerIndex )
					positiveDirection = 1;

				if ( negativeMiddleIndex > positiveUpperIndex )
					negativeDirection = -1;
				else if ( negativeMiddleIndex < positiveLowerIndex )
					negativeDirection = 1;

				if ( ( positiveDirection == 0 ) && ( negativeDirection == 0 ) ) {

					double const positiveMiddleVolume = positiveMiddleMaximum * positiveMiddleIndex - positiveMiddleSum;
					double const negativeMiddleVolume = negativeMiddleMaximum * negativeMiddleIndex - negativeMiddleSum;

					if ( LIKELY( positiveMiddleVolume + negativeMiddleVolume >= total ) ) {

						if ( positiveMiddleIndex > negativeMiddleIndex )
							positiveDirection = -1;
						else if ( negativeMiddleIndex > positiveMiddleIndex )
							negativeDirection = -1;
						else if ( ( positiveUpperIndex - positiveLowerIndex ) > ( negativeUpperIndex - negativeLowerIndex ) )
							positiveDirection = -1;
						else
							negativeDirection = -1;
					}
					else {

						if ( positiveMiddleIndex < negativeMiddleIndex )
							positiveDirection = 1;
						else if ( negativeMiddleIndex < positiveMiddleIndex )
							negativeDirection = 1;
						else if ( ( positiveUpperIndex - positiveLowerIndex ) > ( negativeUpperIndex - negativeLowerIndex ) )
							positiveDirection = 1;
						else
							negativeDirection = 1;
					}
				}
			}

			if ( positiveDirection != 0 ) {

				if ( LIKELY( positiveDirection < 0 ) ) {

					pPositiveUpper = pPositiveMiddle - 1;
					positiveMiddleMaximum = positiveLowerMaximum;
					positiveMiddleSum = positiveLowerSum;
				}
				else {

					pPositiveLower = pPositiveMiddle;
					positiveLowerMaximum = positiveMiddleMaximum;
					positiveLowerSum = positiveMiddleSum;
				}

				if ( pPositiveLower < pPositiveUpper ) {

					pPositiveMiddle = Partition( pPositiveLower, pPositiveUpper ) + 1;
					for ( double const* ii = pPositiveLower; ii < pPositiveMiddle; ++ii ) {

						if ( *ii > positiveMiddleMaximum )
							positiveMiddleMaximum = *ii;
						positiveMiddleSum += *ii;
					}
				}
				else
					pPositiveMiddle = pPositiveLower;
			}

			if ( negativeDirection != 0 ) {

				if ( LIKELY( negativeDirection < 0 ) ) {

					pNegativeUpper = pNegativeMiddle - 1;
					negativeMiddleMaximum = negativeLowerMaximum;
					negativeMiddleSum = negativeLowerSum;
				}
				else {

					pNegativeLower = pNegativeMiddle;
					negativeLowerMaximum = negativeMiddleMaximum;
					negativeLowerSum = negativeMiddleSum;
				}

				if ( pNegativeLower < pNegativeUpper ) {

					pNegativeMiddle = Partition( pNegativeLower, pNegativeUpper ) + 1;
					for ( double const* ii = pNegativeLower; ii < pNegativeMiddle; ++ii ) {

						if ( *ii > negativeMiddleMaximum )
							negativeMiddleMaximum = *ii;
						negativeMiddleSum += *ii;
					}
				}
				else
					pNegativeMiddle = pNegativeLower;
			}
		}
		BOOST_ASSERT( pPositiveLower == pPositiveUpper );
		BOOST_ASSERT( pNegativeLower == pNegativeUpper );
		BOOST_ASSERT( ( pPositiveLower - positiveBegin ) == ( pNegativeLower - negativeBegin ) );

		unsigned int const index = ( pPositiveLower - positiveBegin );
		BOOST_ASSERT( index > 0 );

		double const deltaKappaSum = ( total + positiveLowerSum + negativeLowerSum ) / index - positiveLowerMaximum - negativeLowerMaximum;
		BOOST_ASSERT( deltaKappaSum >= 0 );

		double positiveDeltaKappaBound = deltaKappaSum;
		if ( pPositiveLower != positiveEnd )
			positiveDeltaKappaBound = std::min( positiveDeltaKappaBound, *pPositiveLower - positiveLowerMaximum );

		double negativeDeltaKappaBound = deltaKappaSum;
		if ( pNegativeLower != negativeEnd )
			negativeDeltaKappaBound = std::min( negativeDeltaKappaBound, *pNegativeLower - negativeLowerMaximum );

		positiveKappa = positiveLowerMaximum + 0.5 * ( deltaKappaSum + positiveDeltaKappaBound - negativeDeltaKappaBound );
		negativeKappa = negativeLowerMaximum + 0.5 * ( deltaKappaSum - positiveDeltaKappaBound + negativeDeltaKappaBound );
		BOOST_ASSERT( positiveKappa >= positiveLowerMaximum );
		BOOST_ASSERT( negativeKappa >= negativeLowerMaximum );
		BOOST_ASSERT( ( pPositiveLower == positiveEnd ) || ( positiveKappa <= *pPositiveLower ) );
		BOOST_ASSERT( ( pNegativeLower == negativeEnd ) || ( negativeKappa <= *pNegativeLower ) );
	}

	// kappa =  0.5 * ( positiveKappa + negativeKappa )
	// bias  = -0.5 * ( positiveKappa - negativeKappa )

	return std::pair< double, double >( positiveKappa, negativeKappa );
}


}    // namespace _Private


}    // namespace Biased




}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM
