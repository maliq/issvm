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
	\file svm_optimizer_classification_unbiased_sbp.cpp
	\brief SVM::Optimizer::Unbiased::SBP implementation
*/




#include "svm_optimizer_classification_unbiased_sbp.hpp"
#include "svm_optimizer_classification_private_find_water_level.hpp"




namespace SVM {


namespace Optimizer {


namespace Classification {


namespace Unbiased {




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

	FindTotalKappa();

	double const scale = 1 / m_totalKappa;
	double const* const totalAlphas = m_totalAlphas.Get();
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		begin[ ii ] = totalAlphas[ ii ] * scale;
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

	{	FindTotalKappa();

		double const scale = 1 / m_totalKappa;

		double* pDestination    = begin;
		double* pDestinationEnd = end;
		double const* pResponse = m_totalResponses.Get() + trainingSize;
		double const* pLabel    = m_pKernel->Labels()    + trainingSize;
		for ( ; pDestination != pDestinationEnd; ++pDestination, ++pResponse, ++pLabel ) {

			double const response = *pResponse * scale;
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

	unsigned int index = trainingSize;
	{	unsigned int count = 0;
		double const* pResponse    = m_responses.get();
		double const* pResponseEnd = pResponse + trainingSize;
		double const* pLabel = m_pKernel->Labels();
		for ( ; pResponse != pResponseEnd; ++pResponse, ++pLabel ) {

			double const response = ( ( *pLabel > 0 ) ? *pResponse : -*pResponse );
			if ( UNLIKELY( response <= m_kappa ) )
				++count;
		}

		BOOST_ASSERT( ( count > 0 ) && ( count <= trainingSize ) );
		unsigned int sample = 1;
		{	Random::Distribution::DiscreteUniform<> distribution( count - 1 );
			sample = distribution.Sample( generator ) + 1;
		}
		BOOST_ASSERT( ( sample > 0 ) && ( sample <= count ) );

		pResponse = m_responses.get();
		pLabel = m_pKernel->Labels();
		for ( ; sample > 0; ++pResponse, ++pLabel ) {

			BOOST_ASSERT( pResponse != pResponseEnd );
			double const response = ( ( *pLabel > 0 ) ? *pResponse : -*pResponse );
			if ( UNLIKELY( response <= m_kappa ) )
				if ( --sample == 0 )
					index = ( pResponse - m_responses.get() );
		}
	}
	BOOST_ASSERT( index < trainingSize );

	++m_iterations;
	double const eta = 1 / std::sqrt( m_maximumNormSquared * m_iterations );

	if ( m_pKernel->Labels()[ index ] > 0 ) {

		m_normSquared += eta * ( 2 * m_responses[ index ] + eta * m_normsSquared[ index ] );
		m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), index, m_alphas[ index ] + eta );
	}
	else {

		m_normSquared += eta * ( -2 * m_responses[ index ] + eta * m_normsSquared[ index ] );
		m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), index, m_alphas[ index ] - eta );
	}

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
		double* pDestination    = responses.get();
		double* pDestinationEnd = pDestination + trainingSize;
		double const* pResponse = m_responses.get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pDestination != pDestinationEnd; ++pResponse, ++pLabel ) {

			if ( *pLabel > 0 )
				*( pDestination++ ) = *pResponse;
			else
				*( pDestination++ ) = -*pResponse;
		}

		m_kappa = _Private::FindWaterLevel( responses.get(), responses.get() + trainingSize, m_nu * trainingSize );
	}

	m_totalAlphas.Add(    m_alphas.get()    );
	m_totalResponses.Add( m_responses.get() );

	m_totalKappa = std::numeric_limits< double >::quiet_NaN();

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
		double* pDestination    = responses.get();
		double* pDestinationEnd = pDestination + trainingSize;
		double const* pResponse = m_responses.get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pDestination != pDestinationEnd; ++pResponse, ++pLabel ) {

			if ( *pLabel > 0 )
				*( pDestination++ ) = *pResponse;
			else
				*( pDestination++ ) = -*pResponse;
		}

		m_kappa = _Private::FindWaterLevel( responses.get(), responses.get() + trainingSize, m_nu * trainingSize );
	}

	m_pKernel->RecalculateResponses( m_totalAlphas.Get(), m_totalResponses.Get() );

	m_totalKappa = std::numeric_limits< double >::quiet_NaN();

	m_totalNormSquared = std::numeric_limits< double >::quiet_NaN();
	m_totalAverageLoss = std::numeric_limits< double >::quiet_NaN();
}


void SBP::FindTotalKappa() const {

	if ( boost::math::isnan( m_totalKappa ) ) {

		unsigned int const trainingSize = m_pKernel->TrainingSize();

		boost::scoped_array< double > responses( new double[ trainingSize ] );
		double* pDestination    = responses.get();
		double* pDestinationEnd = pDestination + trainingSize;
		double const* pResponse = m_totalResponses.Get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pDestination != pDestinationEnd; ++pResponse, ++pLabel ) {

			if ( *pLabel > 0 )
				*( pDestination++ ) = *pResponse;
			else
				*( pDestination++ ) = -*pResponse;
		}

		m_totalKappa = _Private::FindWaterLevel( responses.get(), responses.get() + trainingSize, m_nu * trainingSize * m_iterations );
	}
	BOOST_ASSERT( ! boost::math::isnan( m_totalKappa ) );
}


void SBP::FindTotalLosses() const {

	BOOST_ASSERT( boost::math::isnan( m_totalNormSquared ) == boost::math::isnan( m_totalAverageLoss ) );
	if ( boost::math::isnan( m_totalNormSquared ) ) {

		unsigned int const trainingSize = m_pKernel->TrainingSize();

		FindTotalKappa();

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

			double const classification = ( ( *pLabel > 0 ) ? response : -response ) / m_totalKappa;
			if ( classification < 1 )
				totalLoss += 1 - classification;
		}
		BOOST_ASSERT( normSquared >= 0 );
		BOOST_ASSERT( totalLoss   >= 0 );
		m_totalNormSquared = normSquared / m_totalKappa;
		m_totalAverageLoss = totalLoss / trainingSize;
	}
	BOOST_ASSERT( ! boost::math::isnan( m_totalNormSquared ) );
	BOOST_ASSERT( ! boost::math::isnan( m_totalAverageLoss ) );
}




}    // namespace Unbiased


}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM
