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
	\file svm_optimizer_classification_biased_sbp.cpp
	\brief SVM::Optimizer::Biased::SBP implementation
*/




#include "svm_optimizer_classification_biased_sbp.hpp"
#include "svm_optimizer_classification_private_find_water_level.hpp"




namespace SVM {


namespace Optimizer {


namespace Classification {


namespace Biased {




//============================================================================
//    SBP methods
//============================================================================


SBP::~SBP() {
}


unsigned int const SBP::TrainingSize() const {

	return m_pKernel->TrainingSize();
}


unsigned int const SBP::ValidationSize() const {

	return( m_pKernel->Size() - m_pKernel->TrainingSize() );
}


void SBP::GetAlphas( double* const begin, double* const end ) const {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	if ( end != begin + trainingSize )
		throw std::runtime_error( "SBP::GetAlphas parameter array has incorrect length" );

	FindTotalKappas();

	double const kappa = 0.5 * ( m_totalPositiveKappa + m_totalNegativeKappa );

	double const scale = 1 / kappa;
	double const* const totalAlphas = m_totalAlphas.Get();
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		begin[ ii ] = totalAlphas[ ii ] * scale;
}


double const SBP::Bias() const {

	FindTotalKappas();

	double const kappa =  0.5 * ( m_totalPositiveKappa + m_totalNegativeKappa );
	double const bias  = -0.5 * ( m_totalPositiveKappa - m_totalNegativeKappa );

	return( bias / kappa );
}


double const SBP::NormSquared() const {

	FindTotalLosses();

	return m_totalNormSquared;
}


void SBP::GetValidationResponses( double* const begin, double* const end ) const {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	BOOST_ASSERT( trainingSize <= totalSize );
	if ( trainingSize >= totalSize )
		throw std::runtime_error( "SBP::GetValidationResponses requires a validation set" );
	if ( end != begin + totalSize - trainingSize )
		throw std::runtime_error( "SBP::GetValidationResponses parameter array has incorrect length" );

	{	FindTotalKappas();

		double const kappa =  0.5 * ( m_totalPositiveKappa + m_totalNegativeKappa );
		double const bias  = -0.5 * ( m_totalPositiveKappa - m_totalNegativeKappa );
		double const scale = 1 / kappa;

		double* pDestination    = begin;
		double* pDestinationEnd = end;
		double const* pResponse = m_totalResponses.Get() + trainingSize;
		double const* pLabel    = m_pKernel->Labels()    + trainingSize;
		for ( ; pDestination != pDestinationEnd; ++pDestination, ++pResponse, ++pLabel ) {

			double const response = ( *pResponse + bias ) * scale;
			*pDestination = ( ( *pLabel > 0 ) ? response : -response );
		}
	}
}


double const SBP::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< float > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SBP::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< double > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SBP::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector< float > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SBP::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector< double > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SBP::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector< float > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SBP::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector< double > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


void SBP::Iterate( Random::Generator::LaggedFibonacci4<>& generator ) {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	unsigned int positiveIndex = trainingSize;
	unsigned int negativeIndex = trainingSize;
	{	unsigned int positiveCount = 0;
		unsigned int positiveLessCount = 0;
		unsigned int negativeCount = 0;
		unsigned int negativeLessCount = 0;
		double const* pResponse    = m_responses.get();
		double const* pResponseEnd = pResponse + trainingSize;
		double const* pLabel = m_pKernel->Labels();
		for ( ; pResponse != pResponseEnd; ++pResponse, ++pLabel ) {

			if ( *pLabel > 0 ) {

				double const response = *pResponse;
				if ( UNLIKELY( response <= m_positiveKappa ) ) {

					++positiveCount;
					if ( LIKELY( response < m_positiveKappa ) )
						++positiveLessCount;
				}
			}
			else {

				double const response = -*pResponse;
				if ( UNLIKELY( response <= m_negativeKappa ) ) {

					++negativeCount;
					if ( LIKELY( response < m_negativeKappa ) )
						++negativeLessCount;
				}
			}
		}
		BOOST_ASSERT( positiveLessCount <= negativeCount );
		BOOST_ASSERT( negativeLessCount <= positiveCount );

		unsigned int const count = std::min( positiveCount, negativeCount );
		BOOST_ASSERT( count >= positiveLessCount );
		BOOST_ASSERT( count >= negativeLessCount );

		BOOST_ASSERT( ( count > 0 ) && ( count <= trainingSize ) );
		unsigned int positiveSample = 1;
		unsigned int negativeSample = 1;
		{	Random::Distribution::DiscreteUniform<> distribution( count - 1 );
			positiveSample = distribution.Sample( generator ) + 1;
			negativeSample = distribution.Sample( generator ) + 1;
		}
		BOOST_ASSERT( ( positiveSample > 0 ) && ( positiveSample <= count ) );
		BOOST_ASSERT( ( negativeSample > 0 ) && ( negativeSample <= count ) );

		unsigned int positiveEqual = count - positiveLessCount;
		unsigned int negativeEqual = count - negativeLessCount;

		pResponse = m_responses.get();
		pLabel = m_pKernel->Labels();
		for ( ; positiveSample + negativeSample > 0; ++pResponse, ++pLabel ) {

			BOOST_ASSERT( pResponse != pResponseEnd );
			if ( *pLabel > 0 ) {

				if ( LIKELY( positiveSample > 0 ) ) {

					double const response = *pResponse;
					if ( UNLIKELY( response <= m_positiveKappa ) ) {

						if ( LIKELY( response < m_positiveKappa ) ) {

							if ( --positiveSample == 0 )
								positiveIndex = ( pResponse - m_responses.get() );
						}
						else if ( positiveEqual > 0 ) {

							if ( --positiveSample == 0 )
								positiveIndex = ( pResponse - m_responses.get() );
							--positiveEqual;
						}
					}
				}
			}
			else {

				if ( LIKELY( negativeSample > 0 ) ) {

					double const response = -*pResponse;
					if ( UNLIKELY( response <= m_negativeKappa ) ) {

						if ( LIKELY( response < m_negativeKappa ) ) {

							if ( --negativeSample == 0 )
								negativeIndex = ( pResponse - m_responses.get() );
						}
						else if ( negativeEqual > 0 ) {

							if ( --negativeSample == 0 )
								negativeIndex = ( pResponse - m_responses.get() );
							--negativeEqual;
						}
					}
				}
			}
		}
	}
	BOOST_ASSERT( positiveIndex < trainingSize );
	BOOST_ASSERT( negativeIndex < trainingSize );

	++m_iterations;
	double const eta = 0.5 / std::sqrt( m_maximumNormSquared * m_iterations );

	BOOST_ASSERT( m_pKernel->Labels()[ positiveIndex ] > 0 );
	m_normSquared += eta * ( 2 * m_responses[ positiveIndex ] + eta * m_normsSquared[ positiveIndex ] );
	m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), positiveIndex, m_alphas[ positiveIndex ] + eta );

	BOOST_ASSERT( ! ( m_pKernel->Labels()[ negativeIndex ] > 0 ) );
	m_normSquared += eta * ( -2 * m_responses[ negativeIndex ] + eta * m_normsSquared[ negativeIndex ] );
	m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), negativeIndex, m_alphas[ negativeIndex ] - eta );

	if ( m_normSquared > 1 ) {

		double const scale = 1 / std::sqrt( m_normSquared );

		double* pDestination    = m_alphas.get();
		double* pDestinationEnd = pDestination + trainingSize;
		for ( ; pDestination != pDestinationEnd; ++pDestination )
			*pDestination *= scale;

		pDestination    = m_responses.get();
		pDestinationEnd = pDestination + totalSize;
		for ( ; pDestination != pDestinationEnd; ++pDestination )
			*pDestination *= scale;

		m_normSquared = 1;
	}

	{	boost::scoped_array< double > responses( new double[ trainingSize ] );
		double* pFront = responses.get();
		double* pBack = pFront + trainingSize - 1;
		double const* pResponse = m_responses.get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pFront <= pBack; ++pResponse, ++pLabel ) {

			if ( *pLabel > 0 )
				*( pFront++ ) = *pResponse;
			else
				*( pBack-- ) = -*pResponse;
		}

		std::pair< double, double > const waterLevel = _Private::FindWaterLevel( responses.get(), pFront, pFront, responses.get() + trainingSize, m_nu * trainingSize );
		m_positiveKappa = waterLevel.first;
		m_negativeKappa = waterLevel.second;
	}

	m_totalAlphas.Add(    m_alphas.get()    );
	m_totalResponses.Add( m_responses.get() );

	m_totalPositiveKappa = std::numeric_limits< double >::quiet_NaN();
	m_totalNegativeKappa = std::numeric_limits< double >::quiet_NaN();

	m_totalNormSquared = std::numeric_limits< double >::quiet_NaN();
	m_totalAverageLoss = std::numeric_limits< double >::quiet_NaN();
}


void SBP::Recalculate() {

	unsigned int const trainingSize = m_pKernel->TrainingSize();

	m_pKernel->RecalculateResponses( m_alphas.get(), m_responses.get() );

	{	double accumulator = 0;
		double const* pAlpha    = m_alphas.get();
		double const* pAlphaEnd = pAlpha + trainingSize;
		double const* pResponse = m_responses.get();
		for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pResponse )
			accumulator += *pAlpha * *pResponse;
		BOOST_ASSERT( accumulator >= 0 );
		m_normSquared = accumulator;
	}

	{	boost::scoped_array< double > responses( new double[ trainingSize ] );
		double* pFront = responses.get();
		double* pBack = pFront + trainingSize - 1;
		double const* pResponse = m_responses.get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pFront <= pBack; ++pResponse, ++pLabel ) {

			if ( *pLabel > 0 )
				*( pFront++ ) = *pResponse;
			else
				*( pBack-- ) = -*pResponse;
		}

		std::pair< double, double > const waterLevel = _Private::FindWaterLevel( responses.get(), pFront, pFront, responses.get() + trainingSize, m_nu * trainingSize );
		m_positiveKappa = waterLevel.first;
		m_negativeKappa = waterLevel.second;
	}

	m_pKernel->RecalculateResponses( m_totalAlphas.Get(), m_totalResponses.Get() );

	m_totalPositiveKappa = std::numeric_limits< double >::quiet_NaN();
	m_totalNegativeKappa = std::numeric_limits< double >::quiet_NaN();

	m_totalNormSquared = std::numeric_limits< double >::quiet_NaN();
	m_totalAverageLoss = std::numeric_limits< double >::quiet_NaN();
}


void SBP::FindTotalKappas() const {

	BOOST_ASSERT( boost::math::isnan( m_totalPositiveKappa ) == boost::math::isnan( m_totalNegativeKappa ) );
	if ( boost::math::isnan( m_totalPositiveKappa ) ) {

		unsigned int const trainingSize = m_pKernel->TrainingSize();

		boost::scoped_array< double > responses( new double[ trainingSize ] );
		double* pFront = responses.get();
		double* pBack = pFront + trainingSize - 1;
		double const* pResponse = m_totalResponses.Get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pFront <= pBack; ++pResponse, ++pLabel ) {

			if ( *pLabel > 0 )
				*( pFront++ ) = *pResponse;
			else
				*( pBack-- ) = -*pResponse;
		}

		std::pair< double, double > const waterLevel = _Private::FindWaterLevel( responses.get(), pFront, pFront, responses.get() + trainingSize, m_nu * trainingSize * m_iterations );
		m_totalPositiveKappa = waterLevel.first;
		m_totalNegativeKappa = waterLevel.second;
	}
	BOOST_ASSERT( ! boost::math::isnan( m_totalPositiveKappa ) );
	BOOST_ASSERT( ! boost::math::isnan( m_totalNegativeKappa ) );
}


void SBP::FindTotalLosses() const {

	BOOST_ASSERT( boost::math::isnan( m_totalNormSquared ) == boost::math::isnan( m_totalAverageLoss ) );
	if ( boost::math::isnan( m_totalNormSquared ) ) {

		unsigned int const trainingSize = m_pKernel->TrainingSize();

		FindTotalKappas();

		double const kappa =  0.5 * ( m_totalPositiveKappa + m_totalNegativeKappa );
		double const bias  = -0.5 * ( m_totalPositiveKappa - m_totalNegativeKappa );

		double normSquared = 0;
		double totalLoss   = 0;
		double const* pAlpha    = m_totalAlphas.Get();
		double const* pAlphaEnd = pAlpha + trainingSize;
		double const* pResponse = m_totalResponses.Get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pResponse, ++pLabel ) {

			double const alpha    = *pAlpha;
			double const response = *pResponse;

			normSquared += alpha * response;

			double const classification = ( ( *pLabel > 0 ) ? response + bias : -( response + bias ) ) / kappa;
			if ( classification < 1 )
				totalLoss += 1 - classification;
		}
		BOOST_ASSERT( normSquared >= 0 );
		BOOST_ASSERT( totalLoss   >= 0 );
		m_totalNormSquared = normSquared / kappa;
		m_totalAverageLoss = totalLoss / trainingSize;
	}
	BOOST_ASSERT( ! boost::math::isnan( m_totalNormSquared ) );
	BOOST_ASSERT( ! boost::math::isnan( m_totalAverageLoss ) );
}




}    // namespace Biased


}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM
