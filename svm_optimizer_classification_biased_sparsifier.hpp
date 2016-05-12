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
	\file svm_optimizer_classification_biased_sparsifier.hpp
	\brief SVM::Optimizer::Biased::Sparsifier implementation
*/




#ifndef __SVM_OPTIMIZER_CLASSIFICATION_BIASED_SPARSIFIER_HPP__
#define __SVM_OPTIMIZER_CLASSIFICATION_BIASED_SPARSIFIER_HPP__

#ifdef __cplusplus




#include "svm_optimizer_base.hpp"
#include "svm_kernel.hpp"
#include "random.hpp"
#include "vector.hpp"
#include "array_sum.hpp"
#include "data.hpp"
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
//    Sparsifier class
//============================================================================


struct Sparsifier : public Base {

	inline Sparsifier();

	inline Sparsifier(
		boost::shared_ptr< Kernel::Base > pKernel,
		std::vector< std::string > const& optimizerParameters
	);

	virtual ~Sparsifier();


	virtual unsigned int const TrainingSize()   const;
	virtual unsigned int const ValidationSize() const;

	virtual void GetAlphas( double* const begin, double* const end ) const;
	virtual double const Bias() const;

	virtual double const NormSquared() const;

	virtual void GetValidationResponses( double* const begin, double* const end ) const;
	#if 0    // account for randomness
	virtual double const ValidationError() const;
	#endif    // account for randomness

	inline double const Objective() const;


	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< float  > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< double > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector<   float  > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector<   double > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector<  float  > const& otherVector ) const;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector<  double > const& otherVector ) const;


	virtual void Iterate( Random::Generator::LaggedFibonacci4<>& generator );

	virtual void Recalculate();

	double const GetEta() const {
        return m_eta;
    }

    double const getSubOptimality() const{
        return m_targetSuboptimality;
    }

    double const GetNormSquared() const{
        return m_targetNormSquared;
    }



private:

	template< typename t_VectorType >
	inline double const EvaluateHelper( Random::Generator::LaggedFibonacci4<>& generator, t_VectorType const& vector ) const;

	void FindKappas() const;


	template< typename t_Archive >
	inline void save( t_Archive& archive, unsigned int const ) const;

	template< typename t_Archive >
	inline void load( t_Archive& archive, unsigned int const );

	BOOST_SERIALIZATION_SPLIT_MEMBER()


	boost::shared_ptr< Kernel::Base > m_pKernel;

	double m_eta;
	double m_targetSuboptimality;
	double m_targetNormSquared;
	boost::scoped_array< double > m_targetPredictions;
	unsigned int m_iterations;

	double m_normSquared;
	unsigned int mutable m_positiveIndex;
	unsigned int mutable m_negativeIndex;
	double mutable m_positiveKappa;
	double mutable m_negativeKappa;

	boost::scoped_array< double > m_alphas;
	boost::scoped_array< double > m_responses;

	boost::scoped_array< double > m_normsSquared;

	double m_maximumNormSquared;


	friend class boost::serialization::access;
};




//============================================================================
//    Sparsifier inline methods
//============================================================================


Sparsifier::Sparsifier() :
	m_eta( 1 ),
	m_targetSuboptimality( std::numeric_limits< double >::infinity() ),
	m_targetNormSquared( 0 ),
	m_iterations( 0 ),
	m_normSquared( 0 ),
	m_positiveIndex( 0 ),
	m_negativeIndex( 0 ),
	m_positiveKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_negativeKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_maximumNormSquared( 0 )
{
}


Sparsifier::Sparsifier(
	boost::shared_ptr< Kernel::Base > pKernel,
	std::vector< std::string > const& optimizerParameters
) :
	m_pKernel( pKernel ),
	m_eta( 1 ),
	m_targetSuboptimality( std::numeric_limits< double >::infinity() ),
	m_targetNormSquared( 0 ),
	m_targetPredictions( new double[ m_pKernel->TrainingSize() ] ),
	m_iterations( 0 ),
	m_normSquared( 0 ),
	m_positiveIndex( m_pKernel->TrainingSize() ),
	m_negativeIndex( m_pKernel->TrainingSize() ),
	m_positiveKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_negativeKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_alphas( new double[ m_pKernel->TrainingSize() ] ),
	m_responses( new double[ m_pKernel->Size() ] ),
	m_normsSquared( new double[ m_pKernel->TrainingSize() ] ),
	m_maximumNormSquared( 0 )
{
	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	std::cout << "I am Biased Sparsifier" << std::endl;

	if ( ( optimizerParameters.size() < 2 ) || ( optimizerParameters.size() > 4 ) ) {

		std::stringstream stream;
		stream << "Biased::Sparsifier optimizer expects between two and four parameters (prediction_file,norm^2,eta=1,suboptimality=inf), received [";
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
	{	std::vector< double > data;
		LoadVector( optimizerParameters[ 0 ], data );
		if ( data.size() != trainingSize )
			throw std::runtime_error( "Biased::Sparsifier optimizer parameter prediction_file must contain a prediction for every training example" );
		std::copy( data.begin(), data.end(), m_targetPredictions.get() );
	}
	m_targetNormSquared = atof( optimizerParameters[ 1 ].c_str() );
	if ( m_targetNormSquared < 0 )
		throw std::runtime_error( "Biased::Sparsifier optimizer parameter norm^2 must be nonnegative" );
	if ( optimizerParameters.size() > 2 ) {

		m_eta = atof( optimizerParameters[ 2 ].c_str() );
		if ( m_eta <= 0 )
			throw std::runtime_error( "Biased::Sparsifier optimizer parameter eta must be positive" );
		if ( optimizerParameters.size() > 3 )
			m_targetSuboptimality = atof( optimizerParameters[ 3 ].c_str() );
	}

	std::fill( m_alphas.get(),    m_alphas.get()    + trainingSize, 0 );
	std::fill( m_responses.get(), m_responses.get() + totalSize,    0 );

	for ( unsigned int ii = 0; ii < trainingSize; ++ii ) {

		double const normSquared = m_pKernel->KernelInnerProduct( ii, ii );
		if ( normSquared > m_maximumNormSquared )
			m_maximumNormSquared = normSquared;
		m_normsSquared[ ii ] = normSquared;
	}
}


double const Sparsifier::Objective() const {

	FindKappas();

	double const kappa = 0.5 * ( m_positiveKappa + m_negativeKappa );
	return kappa;
}


template< typename t_VectorType >
double const Sparsifier::EvaluateHelper( Random::Generator::LaggedFibonacci4<>& generator, t_VectorType const& vector ) const {

	FindKappas();

	double const bias = 0.5 * ( m_positiveKappa - m_negativeKappa );
#if 0    // include randomness
	Random::Distribution::StandardUniform<> distribution;
	return( m_pKernel->Evaluate( vector, m_alphas.get() ) + bias + distribution.Sample( generator ) - 0.5 );
#else    // exclude randomness
	return( m_pKernel->Evaluate( vector, m_alphas.get() ) + bias );
#endif
}


template< typename t_Archive >
void Sparsifier::save( t_Archive& archive, unsigned int const ) const {

	archive & boost::serialization::base_object< Base >( *this );

	archive & m_pKernel;
	archive & m_eta;
	archive & m_targetSuboptimality;
	archive & m_targetNormSquared;

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_targetPredictions[ ii ];

	archive & m_iterations;
	archive & m_normSquared;

	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_alphas[ ii ];

	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & m_responses[ ii ];
}


template< typename t_Archive >
void Sparsifier::load( t_Archive& archive, unsigned int const ) {

	archive & boost::serialization::base_object< Base >( *this );

	archive & m_pKernel;
	archive & m_eta;
	archive & m_targetSuboptimality;
	archive & m_targetNormSquared;

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	m_targetPredictions.reset( new double[ trainingSize ] );
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_targetPredictions[ ii ];

	archive & m_iterations;
	archive & m_normSquared;

	m_alphas.reset( new double[ trainingSize ] );
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_alphas[ ii ];

	m_responses.reset( new double[ totalSize ] );
	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & m_responses[ ii ];

	m_normsSquared.reset( new double[ trainingSize ] );
	m_maximumNormSquared = 0;
	for ( unsigned int ii = 0; ii < trainingSize; ++ii ) {

		double const normSquared = m_pKernel->KernelInnerProduct( ii, ii );
		if ( normSquared > m_maximumNormSquared )
			m_maximumNormSquared = normSquared;
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

#endif    /* __SVM_OPTIMIZER_CLASSIFICATION_BIASED_SPARSIFIER_HPP__ */
