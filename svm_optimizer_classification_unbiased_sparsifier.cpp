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
	\file svm_optimizer_classification_unbiased_sparsifier.cpp
	\brief SVM::Optimizer::Unbiased::Sparsifier implementation
*/




#include "svm_optimizer_classification_unbiased_sparsifier.hpp"




namespace SVM {


namespace Optimizer {


namespace Classification {


namespace Unbiased {




//============================================================================
//    Sparsifier methods
//============================================================================


Sparsifier::~Sparsifier() {
}


unsigned int const Sparsifier::TrainingSize() const {

	return m_pKernel->TrainingSize();
}


unsigned int const Sparsifier::ValidationSize() const {

	return( m_pKernel->Size() - m_pKernel->TrainingSize() );
}


void Sparsifier::GetAlphas( double* const begin, double* const end ) const {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	if ( end != begin + trainingSize )
		throw std::runtime_error( "Sparsifier::GetAlphas parameter array has incorrect length" );

	std::copy( m_alphas.get(), m_alphas.get() + trainingSize, begin );
}


double const Sparsifier::NormSquared() const {

	return m_normSquared;
}


void Sparsifier::GetValidationResponses( double* const begin, double* const end ) const {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	BOOST_ASSERT( trainingSize <= totalSize );
	if ( trainingSize >= totalSize )
		throw std::runtime_error( "Sparsifier::GetValidationResponses requires a validation set" );
	if ( end != begin + totalSize - trainingSize )
		throw std::runtime_error( "Sparsifier::GetValidationResponses parameter array has incorrect length" );

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


#if 0    // account for randomness

double const Sparsifier::ValidationError() const {

	unsigned int const size = ValidationSize();
	boost::shared_array< double > responses( new double[ size ] );
	GetValidationResponses( responses.get(), responses.get() + size );

	double error = 0;
	{	double const* ii    = responses.get();
		double const* iiEnd = ii + size;
		for ( ; ii != iiEnd; ++ii ) {

			if ( *ii <= -0.5 )
				++error;
			else if ( *ii < 0.5 )
				error += 0.5 - *ii;
		}
	}

	return( error / size );
}

#endif    // account for randomness


double const Sparsifier::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< float > const& otherVector ) const {

	return EvaluateHelper( generator, otherVector );
}


double const Sparsifier::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< double > const& otherVector ) const {

	return EvaluateHelper( generator, otherVector );
}


double const Sparsifier::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector< float > const& otherVector ) const {

	return EvaluateHelper( generator, otherVector );
}


double const Sparsifier::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector< double > const& otherVector ) const {

	return EvaluateHelper( generator, otherVector );
}


double const Sparsifier::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector< float > const& otherVector ) const {

	return EvaluateHelper( generator, otherVector );
}


double const Sparsifier::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector< double > const& otherVector ) const {

	return EvaluateHelper( generator, otherVector );
}


void Sparsifier::Iterate( Random::Generator::LaggedFibonacci4<>& generator ) {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	FindKappa();

	// find an update pair in the support set which is more than m_targetSuboptimality-suboptimal
	unsigned int index = m_index;
	{	double kappa = m_targetSuboptimality;

		double const* pTargetPrediction    = m_targetPredictions.get();
		double const* pTargetPredictionEnd = pTargetPrediction + trainingSize;
		double const* pAlpha = m_alphas.get();
		double const* pResponse = m_responses.get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pTargetPrediction != pTargetPredictionEnd; ++pTargetPrediction, ++pAlpha, ++pResponse, ++pLabel ) {

			if ( *pAlpha != 0 ) {

				double const targetPrediction = std::min( 1.0, ( ( *pLabel > 0 ) ? *pTargetPrediction : -*pTargetPrediction ) );
				if ( targetPrediction > 0 ) {

					double const response = ( ( *pLabel > 0 ) ? *pResponse : -*pResponse );
					double const delta = targetPrediction - response;
					if ( UNLIKELY( delta > kappa ) ) {

						index = ( pTargetPrediction - m_targetPredictions.get() );
						kappa = delta;
					}
				}
			}
		}
	}

	++m_iterations;
	double const eta = m_eta / m_maximumNormSquared;

	if ( m_pKernel->Labels()[ index ] > 0 ) {

		m_normSquared += eta * ( 2 * m_responses[ index ] + eta * m_normsSquared[ index ] );
		m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), index, m_alphas[ index ] + eta );
	}
	else {

		m_normSquared += eta * ( -2 * m_responses[ index ] + eta * m_normsSquared[ index ] );
		m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), index, m_alphas[ index ] - eta );
	}

	if ( m_normSquared > m_targetNormSquared ) {

		double const scale = std::sqrt( m_targetNormSquared / m_normSquared );

		double* pDestination    = m_alphas.get();
		double* pDestinationEnd = pDestination + trainingSize;
		for ( ; pDestination != pDestinationEnd; ++pDestination )
			*pDestination *= scale;

		pDestination    = m_responses.get();
		pDestinationEnd = pDestination + totalSize;
		for ( ; pDestination != pDestinationEnd; ++pDestination )
			*pDestination *= scale;

		m_normSquared = m_targetNormSquared;
	}

	m_index = trainingSize;
	m_kappa = std::numeric_limits< double >::quiet_NaN();
}


void Sparsifier::Recalculate() {

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

	m_index = trainingSize;
	m_kappa = std::numeric_limits< double >::quiet_NaN();
}


void Sparsifier::FindKappa() const {

	if ( boost::math::isnan( m_kappa ) ) {

        unsigned int const trainingSize = m_pKernel->TrainingSize();

		m_index = trainingSize;
		m_kappa = -std::numeric_limits< double >::infinity();

		double const* pTargetPrediction    = m_targetPredictions.get();
		double const* pTargetPredictionEnd = pTargetPrediction + trainingSize;
		double const* pResponse = m_responses.get();
		double const* pLabel = m_pKernel->Labels();
		for ( ; pTargetPrediction != pTargetPredictionEnd; ++pTargetPrediction, ++pResponse, ++pLabel ) {

			double const targetPrediction = std::min( 1.0, ( ( *pLabel > 0 ) ? *pTargetPrediction : -*pTargetPrediction ) );
			if ( targetPrediction > 0 ) {

				double const response = ( ( *pLabel > 0 ) ? *pResponse : -*pResponse );
				double const delta = targetPrediction - response;
				if ( UNLIKELY( delta > m_kappa ) ) {

					m_index = ( pTargetPrediction - m_targetPredictions.get() );
					m_kappa = delta;
				}
			}
		}
	}
	BOOST_ASSERT( m_index < m_pKernel->TrainingSize() );
	BOOST_ASSERT( ! boost::math::isnan( m_kappa ) );
}




}    // namespace Unbiased


}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM
