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
	\file svm_optimizer_classification_biased_smo.cpp
	\brief SVM::Optimizer::Biased::SMO implementation
*/




#include "svm_optimizer_classification_biased_smo.hpp"




namespace SVM {


namespace Optimizer {


namespace Classification {


namespace Biased {




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


double const SMO::Bias() const {

	FindLosses();

	return m_bias;
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

	{	FindLosses();

		double* pDestination    = begin;
		double* pDestinationEnd = end;
		double const* pResponse = m_responses.get()   + trainingSize;
		double const* pLabel    = m_pKernel->Labels() + trainingSize;
		for ( ; pDestination != pDestinationEnd; ++pDestination, ++pResponse, ++pLabel ) {

			double const response = *pResponse + m_bias;
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

	unsigned int bestIndex1 = 0;
	double bestScore = -std::numeric_limits< double >::infinity();

	{	double const* pAlpha    = m_alphas.get();
		double const* pAlphaEnd = pAlpha + trainingSize;
		double const* pResponse = m_responses.get();
		double const* pLabel    = m_pKernel->Labels();
		for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pResponse, ++pLabel ) {

			double const alpha    = *pAlpha;
			double const response = *pResponse;
			double const sign = ( ( *pLabel > 0 ) ? 1 : -1 );

			double const gradient = sign - response;

			double score = std::fabs( gradient );
			if ( *pLabel > 0 ) {

				if (
					( ( gradient > 0 ) && ( alpha >= upperBound ) ) ||
					( ( gradient < 0 ) && ( alpha <= 0          ) )
				)
				{
					score = -score;
				}
			}
			else {

				if (
					( ( gradient > 0 ) && ( alpha >= 0           ) ) ||
					( ( gradient < 0 ) && ( alpha <= -upperBound ) )
				)
				{
					score = -score;
				}
			}

			if ( UNLIKELY( score > bestScore ) ) {

				bestIndex1 = ( pAlpha - m_alphas.get() );
				bestScore = score;
			}
		}
	}
	boost::shared_array< double > const row1( m_pKernel->Row( bestIndex1 ) );

	unsigned int bestIndex2 = 0;
	double bestAlpha1 = std::numeric_limits< double >::quiet_NaN();
	double bestAlpha2 = std::numeric_limits< double >::quiet_NaN();
	bestScore = -std::numeric_limits< double >::infinity();

	{	double const k11       = m_normsSquared[ bestIndex1 ];
		double const alpha1    = m_alphas[       bestIndex1 ];
		double const response1 = m_responses[    bestIndex1 ];
		double const sign1 = ( ( m_pKernel->Labels()[ bestIndex1 ] > 0 ) ? 1 : -1 );

		double const* pAlpha       = m_alphas.get();
		double const* pResponse    = m_responses.get();
		double const* pNormSquared = m_normsSquared.get();
		double const* pRow1        = row1.get();
		double const* pLabel       = m_pKernel->Labels();
		for ( unsigned int ii = 0; ii < trainingSize; ++ii, ++pAlpha, ++pResponse, ++pNormSquared, ++pRow1, ++pLabel ) {

			if ( LIKELY( ii != bestIndex1 ) ) {

				double const k12 = *pRow1;
				double const k22 = *pNormSquared;
				double const alpha2    = *pAlpha;
				double const response2 = *pResponse;
				double const sign2 = ( ( *pLabel > 0 ) ? 1 : -1 );

				double const numerator = ( sign1 - response1 ) - ( sign2 - response2 );
				double const denominator = k11 + k22 - 2 * k12;
				BOOST_ASSERT( denominator >= 0 );
				double delta = ( ( std::fabs( numerator ) < denominator * upperBound ) ? ( numerator / denominator ) : ( ( numerator >= 0 ) ? upperBound : -upperBound ) );

				double newAlpha1 = alpha1 + delta;
				if ( sign1 > 0 )
					newAlpha1 = std::min( upperBound, std::max( 0.0, newAlpha1 ) );
				else
					newAlpha1 = std::max( -upperBound, std::min( 0.0, newAlpha1 ) );

				double newAlpha2 = alpha2 - delta;
				if ( sign2 > 0 )
					newAlpha2 = std::min( upperBound, std::max( 0.0, newAlpha2 ) );
				else
					newAlpha2 = std::max( -upperBound, std::min( 0.0, newAlpha2 ) );

				double const newDelta1 = newAlpha1 - alpha1;
				double const newDelta2 = alpha2 - newAlpha2;
				if ( std::fabs( newDelta1 ) < std::fabs( newDelta2 ) ) {

					delta = newDelta1;
					newAlpha2 = alpha2 - delta;
				}
				else if ( std::fabs( newDelta2 ) < std::fabs( newDelta1 ) ) {

					delta = newDelta2;
					newAlpha1 = alpha1 + delta;
				}
				else {

					BOOST_ASSERT( newDelta1 == newDelta2 );
					delta = newDelta1;
				}
				BOOST_ASSERT( sign1 * newAlpha1 >= 0          );
				BOOST_ASSERT( sign1 * newAlpha1 <= upperBound );
				BOOST_ASSERT( sign2 * newAlpha2 >= 0          );
				BOOST_ASSERT( sign2 * newAlpha2 <= upperBound );

				double const score = delta * ( numerator - 0.5 * delta * denominator );
				BOOST_ASSERT( score >= 0 );

				if ( UNLIKELY( score > bestScore ) ) {

					bestIndex2 = ii;
					bestAlpha1 = newAlpha1;
					bestAlpha2 = newAlpha2;
					bestScore = score;
				}
			}
		}
	}

	m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), bestIndex1, bestAlpha1, row1 );
	m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), bestIndex2, bestAlpha2 );

	m_bias        = std::numeric_limits< double >::quiet_NaN();
	m_normSquared = std::numeric_limits< double >::quiet_NaN();
	m_averageLoss = std::numeric_limits< double >::quiet_NaN();
	m_alphaSum    = std::numeric_limits< double >::quiet_NaN();
}


void SMO::Recalculate() {

	m_pKernel->RecalculateResponses( m_alphas.get(), m_responses.get() );

	m_bias        = std::numeric_limits< double >::quiet_NaN();
	m_normSquared = std::numeric_limits< double >::quiet_NaN();
	m_averageLoss = std::numeric_limits< double >::quiet_NaN();
	m_alphaSum    = std::numeric_limits< double >::quiet_NaN();
}


void SMO::FindLosses() const {

	BOOST_ASSERT( boost::math::isnan( m_bias        ) == boost::math::isnan( m_normSquared ) );
	BOOST_ASSERT( boost::math::isnan( m_normSquared ) == boost::math::isnan( m_averageLoss ) );
	BOOST_ASSERT( boost::math::isnan( m_averageLoss ) == boost::math::isnan( m_alphaSum    ) );
	if ( boost::math::isnan( m_bias ) ) {

		unsigned int const trainingSize = m_pKernel->TrainingSize();

		double const upperBound = 1.0 / ( m_lambda * trainingSize );

		double biasNumerator = 0;
		unsigned int biasDenominator = 0;
		double const* pAlpha    = m_alphas.get();
		double const* pAlphaEnd = pAlpha + trainingSize;
		double const* pResponse = m_responses.get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pResponse, ++pLabel ) {

			double const alpha = *pAlpha;
			if ( ( alpha != 0 ) && ( std::abs( alpha ) < upperBound ) ) {

				double const response = *pResponse;

				biasNumerator += ( ( *pLabel > 0 ) ? 1 : -1 ) - response;
				++biasDenominator;
			}
		}
		m_bias = ( ( biasDenominator > 0 ) ? ( biasNumerator / biasDenominator ) : 0 );

		double normSquared = 0;
		double totalLoss   = 0;
		double alphaSum    = 0;
		pAlpha    = m_alphas.get();
		pResponse = m_responses.get();
		pLabel = m_pKernel->Labels();
		for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pResponse, ++pLabel ) {

			double const alpha    = *pAlpha;
			double const response = *pResponse;

			normSquared += alpha * response;

			double const loss = 1 - ( ( *pLabel > 0 ) ? response + m_bias : -( response + m_bias ) );
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
	BOOST_ASSERT( ! boost::math::isnan( m_bias ) );
	BOOST_ASSERT( ! boost::math::isnan( m_normSquared ) );
	BOOST_ASSERT( ! boost::math::isnan( m_averageLoss ) );
	BOOST_ASSERT( ! boost::math::isnan( m_alphaSum ) );
}




}    // namespace Biased


}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM
