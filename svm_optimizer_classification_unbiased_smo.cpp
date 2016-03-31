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
	\file svm_optimizer_classification_unbiased_smo.cpp
	\brief SVM::Optimizer::Unbiased::SMO implementation
*/




#include "svm_optimizer_classification_unbiased_smo.hpp"




namespace SVM {


namespace Optimizer {


namespace Classification {


namespace Unbiased {




//============================================================================
//    SMO methods
//============================================================================


SMO::~SMO() {
}


unsigned int const SMO::TrainingSize() const {

	return m_pKernel->TrainingSize();
}


unsigned int const SMO::ValidationSize() const {

	return( m_pKernel->Size() - m_pKernel->TrainingSize() );
}


void SMO::GetAlphas( double* const begin, double* const end ) const {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	if ( end != begin + trainingSize )
		throw std::runtime_error( "SMO::GetAlphas parameter array has incorrect length" );

	std::copy( m_alphas.get(), m_alphas.get() + trainingSize, begin );
}

double const SMO::NormSquared() const {

	FindLosses();

	return m_normSquared;
}


void SMO::GetValidationResponses( double* const begin, double* const end ) const {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	BOOST_ASSERT( trainingSize <= totalSize );
	if ( trainingSize >= totalSize )
		throw std::runtime_error( "SMO::GetValidationResponses requires a validation set" );
	if ( end != begin + totalSize - trainingSize )
		throw std::runtime_error( "SMO::GetValidationResponses parameter array has incorrect length" );

	{	double* pDestination    = begin;
		double* pDestinationEnd = end;
		double const* pResponse = m_responses.get()   + trainingSize;
		double const* pLabel    = m_pKernel->Labels() + trainingSize;
		for ( ; pDestination != pDestinationEnd; ++pDestination, ++pResponse, ++pLabel ) {

			double const response = *pResponse;
			*pDestination = ( ( *pLabel > 0 ) ? response : -response );
		}
	}
}


double const SMO::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< float > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SMO::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< double > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SMO::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector< float > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SMO::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector< double > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SMO::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector< float > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const SMO::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector< double > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


void SMO::Iterate( Random::Generator::LaggedFibonacci4<>& generator ) {

	unsigned int const trainingSize = m_pKernel->TrainingSize();

	double const upperBound = 1.0 / ( m_lambda * trainingSize );

	unsigned int bestIndex = 0;
	double bestAlpha = std::numeric_limits< double >::quiet_NaN();
	double bestScore = -std::numeric_limits< double >::infinity();

	{	double const* pAlpha       = m_alphas.get();
		double const* pAlphaEnd    = pAlpha + trainingSize;
		double const* pResponse    = m_responses.get();
		double const* pNormSquared = m_normsSquared.get();
		double const* pLabel       = m_pKernel->Labels();
		for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pResponse, ++pNormSquared, ++pLabel ) {

			double const alpha       = *pAlpha;
			double const response    = *pResponse;
			double const normSquared = *pNormSquared;
			double const sign = ( ( *pLabel > 0 ) ? 1 : -1 );

			BOOST_ASSERT( normSquared > 0 );
			double delta = ( sign - response ) / normSquared;

			double newAlpha = alpha + delta;
			if ( sign > 0 )
				newAlpha = std::min( upperBound, std::max( 0.0, newAlpha ) );
			else
				newAlpha = std::max( -upperBound, std::min( 0.0, newAlpha ) );
			delta = newAlpha - alpha;
			BOOST_ASSERT( sign * newAlpha >= 0          );
			BOOST_ASSERT( sign * newAlpha <= upperBound );

			double const score = delta * ( ( sign - response ) - 0.5 * delta * normSquared );
			if ( UNLIKELY( score > bestScore ) ) {

				bestIndex = ( pAlpha - m_alphas.get() );
				bestAlpha = newAlpha;
				bestScore = score;
			}
		}
	}

	m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), bestIndex, bestAlpha );

	m_normSquared = std::numeric_limits< double >::quiet_NaN();
	m_averageLoss = std::numeric_limits< double >::quiet_NaN();
	m_alphaSum    = std::numeric_limits< double >::quiet_NaN();
}


void SMO::Recalculate() {

	m_pKernel->RecalculateResponses( m_alphas.get(), m_responses.get() );

	m_normSquared = std::numeric_limits< double >::quiet_NaN();
	m_averageLoss = std::numeric_limits< double >::quiet_NaN();
	m_alphaSum    = std::numeric_limits< double >::quiet_NaN();
}


void SMO::FindLosses() const {

	BOOST_ASSERT( boost::math::isnan( m_normSquared ) == boost::math::isnan( m_averageLoss ) );
	BOOST_ASSERT( boost::math::isnan( m_averageLoss ) == boost::math::isnan( m_alphaSum    ) );
	if ( boost::math::isnan( m_normSquared ) ) {

		unsigned int const trainingSize = m_pKernel->TrainingSize();

		double normSquared = 0;
		double totalLoss   = 0;
		double alphaSum    = 0;
		double const* pAlpha    = m_alphas.get();
		double const* pAlphaEnd = pAlpha + trainingSize;
		double const* pResponse = m_responses.get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pResponse, ++pLabel ) {

			double const alpha    = *pAlpha;
			double const response = *pResponse;

			normSquared += alpha * response;

			double const loss = 1 - ( ( *pLabel > 0 ) ? response : -response );
			if ( loss > 0 )
				totalLoss += loss;

			alphaSum += std::abs( alpha );
		}
		BOOST_ASSERT( normSquared >= 0 );
		BOOST_ASSERT( totalLoss   >= 0 );
		m_normSquared = normSquared;
		m_averageLoss = totalLoss / trainingSize;
		m_alphaSum = alphaSum;
	}
	BOOST_ASSERT( ! boost::math::isnan( m_normSquared ) );
	BOOST_ASSERT( ! boost::math::isnan( m_averageLoss ) );
	BOOST_ASSERT( ! boost::math::isnan( m_alphaSum ) );
}




}    // namespace Unbiased


}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM
