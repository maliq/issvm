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
	\file svm_optimizer_classification_biased_perceptron.hpp
	\brief SVM::Optimizer::Biased::Perceptron implementation
*/




#ifndef __SVM_OPTIMIZER_CLASSIFICATION_BIASED_PERCEPTRON_HPP__
#define __SVM_OPTIMIZER_CLASSIFICATION_BIASED_PERCEPTRON_HPP__

#ifdef __cplusplus




#include "svm_optimizer_base.hpp"
#include "svm_kernel.hpp"
#include "random.hpp"
#include "vector.hpp"
#include "array_sum.hpp"
#include "helpers.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/scoped_array.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <sstream>
#include <vector>
#include <stdexcept>




namespace SVM {


namespace Optimizer {


namespace Classification {


namespace Biased {




//============================================================================
//    Perceptron class
//============================================================================


struct Perceptron : public Base {

	inline Perceptron();

	inline Perceptron(
		boost::shared_ptr< Kernel::Base > pKernel,
		std::vector< std::string > const& optimizerParameters
	);

	virtual ~Perceptron();


	virtual unsigned int const TrainingSize()   const;
	virtual unsigned int const ValidationSize() const;

	virtual void GetAlphas( double* const begin, double* const end ) const;
	virtual double const Bias() const;

	virtual double const NormSquared() const;

	virtual void GetValidationResponses( double* const begin, double* const end ) const;

	inline double const Margin() const;


	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< float  > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< double > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector<   float  > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector<   double > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector<  float  > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector<  double > const& otherVector ) const;


	virtual void Iterate( Random::Generator::LaggedFibonacci4<>& generator );

	virtual void Recalculate();


private:

	template< typename t_VectorType >
	inline double const EvaluateHelper( t_VectorType const& vector ) const;

	void FindKappas() const;


	template< typename t_Archive >
	inline void save( t_Archive& archive, unsigned int const ) const;

	template< typename t_Archive >
	inline void load( t_Archive& archive, unsigned int const );

	BOOST_SERIALIZATION_SPLIT_MEMBER()


	boost::shared_ptr< Kernel::Base > m_pKernel;

	double m_margin;

	double m_normSquared;
	unsigned int mutable m_positiveIndex;
	unsigned int mutable m_negativeIndex;
	double mutable m_positiveKappa;
	double mutable m_negativeKappa;

	boost::scoped_array< double > m_alphas;
	boost::scoped_array< double > m_responses;

	boost::scoped_array< double > m_normsSquared;


	friend class boost::serialization::access;
};




//============================================================================
//    Perceptron inline methods
//============================================================================


Perceptron::Perceptron() :
	m_margin( 0 ),
	m_normSquared( 0 ),
	m_positiveIndex( 0 ),
	m_negativeIndex( 0 ),
	m_positiveKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_negativeKappa( std::numeric_limits< double >::quiet_NaN() )
{
}


Perceptron::Perceptron(
	boost::shared_ptr< Kernel::Base > pKernel,
	std::vector< std::string > const& optimizerParameters
) :
	m_pKernel( pKernel ),
	m_margin( 0 ),
	m_normSquared( 0 ),
	m_positiveIndex( m_pKernel->TrainingSize() ),
	m_negativeIndex( m_pKernel->TrainingSize() ),
	m_positiveKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_negativeKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_alphas( new double[ m_pKernel->TrainingSize() ] ),
	m_responses( new double[ m_pKernel->Size() ] ),
	m_normsSquared( new double[ m_pKernel->TrainingSize() ] )
{
	
	std::cout << "I am Biased Perceptron" << std::endl;

	if ( optimizerParameters.size() != 1 ) {

		std::stringstream stream;
		stream << "Biased::Perceptron optimizer expects one parameter (margin), received [";
		{	std::vector< std::string >::const_iterator ii    = optimizerParameters.begin();
			std::vector< std::string >::const_iterator iiEnd = optimizerParameters.end();
			if ( ii != iiEnd ) {

				stream << *ii;
				for ( ++ii; ii != iiEnd; ++ii )
					stream << ',' << *ii;
			}
		}
		stream << ']';
		throw std::runtime_error( stream.str() );
	}
	m_margin = atof( optimizerParameters[ 0 ].c_str() );
	if ( m_margin < 0 )
		throw std::runtime_error( "Biased::Perceptron optimizer parameter margin must be nonnegative" );

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	std::fill( m_alphas.get(),    m_alphas.get()    + trainingSize, 0 );
	std::fill( m_responses.get(), m_responses.get() + totalSize,    0 );

	for ( unsigned int ii = 0; ii < trainingSize; ++ii ) {

		double const normSquared = m_pKernel->KernelInnerProduct( ii, ii );
		m_normsSquared[ ii ] = normSquared;
	}
}


double const Perceptron::Margin() const {

	FindKappas();

	double const kappa = 0.5 * ( m_positiveKappa + m_negativeKappa );
	return( kappa / std::sqrt( m_normSquared ) );
}


template< typename t_VectorType >
double const Perceptron::EvaluateHelper( t_VectorType const& vector ) const {

	FindKappas();

	double const bias = -0.5 * ( m_positiveKappa - m_negativeKappa );
	return( ( m_pKernel->Evaluate( vector, m_alphas.get() ) + bias ) / std::sqrt( m_normSquared ) );
}


template< typename t_Archive >
void Perceptron::save( t_Archive& archive, unsigned int const ) const {

	archive & boost::serialization::base_object< Base >( *this );

	archive & m_pKernel;
	archive & m_margin;
	archive & m_normSquared;

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_alphas[ ii ];

	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & m_responses[ ii ];
}


template< typename t_Archive >
void Perceptron::load( t_Archive& archive, unsigned int const ) {

	archive & boost::serialization::base_object< Base >( *this );

	archive & m_pKernel;
	archive & m_margin;
	archive & m_normSquared;

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	m_alphas.reset( new double[ trainingSize ] );
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_alphas[ ii ];

	m_responses.reset( new double[ totalSize ] );
	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & m_responses[ ii ];

	m_normsSquared.reset( new double[ trainingSize ] );
	for ( unsigned int ii = 0; ii < trainingSize; ++ii ) {

		double const normSquared = m_pKernel->KernelInnerProduct( ii, ii );
		m_normsSquared[ ii ] = normSquared;
	}

	m_positiveIndex = trainingSize;
	m_negativeIndex = trainingSize;
	m_positiveKappa = std::numeric_limits< double >::quiet_NaN();
	m_negativeKappa = std::numeric_limits< double >::quiet_NaN();
}




}    // namespace Biased


}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_OPTIMIZER_CLASSIFICATION_BIASED_PERCEPTRON_HPP__ */
