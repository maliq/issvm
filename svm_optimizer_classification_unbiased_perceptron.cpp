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
	\file svm_optimizer_classification_unbiased_perceptron.cpp
	\brief SVM::Optimizer::Unbiased::Perceptron implementation
*/




#include "svm_optimizer_classification_unbiased_perceptron.hpp"




namespace SVM {


namespace Optimizer {


namespace Classification {


namespace Unbiased {




//============================================================================
//    Perceptron methods
//============================================================================


Perceptron::~Perceptron() {
}


unsigned int const Perceptron::TrainingSize() const {

	return m_pKernel->TrainingSize();
}


unsigned int const Perceptron::ValidationSize() const {

	return( m_pKernel->Size() - m_pKernel->TrainingSize() );
}


void Perceptron::GetAlphas( double* const begin, double* const end ) const {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	if ( end != begin + trainingSize )
		throw std::runtime_error( "Perceptron::GetAlphas parameter array has incorrect length" );

	double const scale = 1 / std::sqrt( m_normSquared );
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		begin[ ii ] = m_alphas[ ii ] * scale;
}


double const Perceptron::NormSquared() const {

	return 1;
}


void Perceptron::GetValidationResponses( double* const begin, double* const end ) const {

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	BOOST_ASSERT( trainingSize <= totalSize );
	if ( trainingSize >= totalSize )
		throw std::runtime_error( "Perceptron::GetValidationResponses requires a validation set" );
	if ( end != begin + totalSize - trainingSize )
		throw std::runtime_error( "Perceptron::GetValidationResponses parameter array has incorrect length" );

	{	double const scale = 1 / std::sqrt( m_normSquared );

		double* pDestination    = begin;
		double* pDestinationEnd = end;
		double const* pResponse = m_responses.get()   + trainingSize;
		double const* pLabel    = m_pKernel->Labels() + trainingSize;
		for ( ; pDestination != pDestinationEnd; ++pDestination, ++pResponse, ++pLabel ) {

			double const response = *pResponse * scale;
			*pDestination = ( ( *pLabel > 0 ) ? response : -response );
		}
	}
}


double const Perceptron::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< float > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const Perceptron::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< double > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const Perceptron::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector< float > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const Perceptron::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector< double > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const Perceptron::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector< float > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


double const Perceptron::Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector< double > const& otherVector ) const {

	return EvaluateHelper( otherVector );
}


void Perceptron::Iterate( Random::Generator::LaggedFibonacci4<>& generator ) {

	unsigned int const trainingSize = m_pKernel->TrainingSize();

	double const scaledMargin = m_margin * std::sqrt( m_normSquared );

	FindKappa();

	if ( m_kappa < scaledMargin ) {

		unsigned int index = m_index;
		double kappa = scaledMargin;

		{	double const* pAlpha    = m_alphas.get();
			double const* pAlphaEnd = pAlpha + trainingSize;
			double const* pResponse = m_responses.get();
			double const* pLabel = m_pKernel->Labels();
			for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pResponse, ++pLabel ) {

				if ( *pAlpha != 0 ) {

					double const response = ( ( *pLabel > 0 ) ? *pResponse : -*pResponse );
					if ( UNLIKELY( response <= kappa ) ) {

						index = ( pAlpha - m_alphas.get() );
						kappa = response;
					}
				}
			}
		}

		if ( m_pKernel->Labels()[ index ] > 0 ) {

			m_normSquared += 2 * m_responses[ index ] + m_normsSquared[ index ];
			m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), index, m_alphas[ index ] + 1 );
		}
		else {

			m_normSquared += -2 * m_responses[ index ] + m_normsSquared[ index ];
			m_pKernel->SetAlpha( m_alphas.get(), m_responses.get(), index, m_alphas[ index ] - 1 );
		}

		m_index = trainingSize;
		m_kappa = std::numeric_limits< double >::quiet_NaN();
	}
}


void Perceptron::Recalculate() {

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


void Perceptron::FindKappa() const {

	if ( boost::math::isnan( m_kappa ) ) {

        unsigned int const trainingSize = m_pKernel->TrainingSize();

		m_index = trainingSize;
		m_kappa = std::numeric_limits< double >::infinity();

		double const* pResponse    = m_responses.get();
		double const* pResponseEnd = pResponse + trainingSize;
		double const* pLabel = m_pKernel->Labels();
		for ( ; pResponse != pResponseEnd; ++pResponse, ++pLabel ) {

			double const response = ( ( *pLabel > 0 ) ? *pResponse : -*pResponse );
			if ( UNLIKELY( response < m_kappa ) ) {

				m_index = ( pResponse - m_responses.get() );
				m_kappa = response;
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
