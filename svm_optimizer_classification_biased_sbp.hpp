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
	\file svm_optimizer_classification_biased_sbp.hpp
	\brief SVM::Optimizer::Biased::SBP implementation
*/




#ifndef __SVM_OPTIMIZER_CLASSIFICATION_BIASED_SBP_HPP__
#define __SVM_OPTIMIZER_CLASSIFICATION_BIASED_SBP_HPP__

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
//    SBP class
//============================================================================


struct SBP : public Base {

	inline SBP();

	inline SBP(
		boost::shared_ptr< Kernel::Base > pKernel,
		std::vector< std::string > const& optimizerParameters
	);

	virtual ~SBP();


	virtual unsigned int const TrainingSize()   const;
	virtual unsigned int const ValidationSize() const;

	virtual void GetAlphas( double* const begin, double* const end ) const;
	virtual double const Bias() const;

	virtual double const NormSquared() const;

	virtual void GetValidationResponses( double* const begin, double* const end ) const;

	inline double const AverageLoss() const;
	inline double const Objective() const;


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

	void FindTotalKappas() const;
	void FindTotalLosses() const;


	template< typename t_Archive >
	inline void save( t_Archive& archive, unsigned int const ) const;

	template< typename t_Archive >
	inline void load( t_Archive& archive, unsigned int const );

	BOOST_SERIALIZATION_SPLIT_MEMBER()


	boost::shared_ptr< Kernel::Base > m_pKernel;

	double m_nu;
	unsigned int m_iterations;

	double m_normSquared;
	double m_positiveKappa;
	double m_negativeKappa;

	boost::scoped_array< double > m_alphas;
	boost::scoped_array< double > m_responses;

	ArraySum< double > m_totalAlphas;
	ArraySum< double > m_totalResponses;

	boost::scoped_array< double > m_normsSquared;

	double m_maximumNormSquared;

	double mutable m_totalPositiveKappa;
	double mutable m_totalNegativeKappa;

	double mutable m_totalNormSquared;
	double mutable m_totalAverageLoss;


	friend class boost::serialization::access;
};




//============================================================================
//    SBP inline methods
//============================================================================


SBP::SBP() :
	m_nu( 0 ),
	m_iterations( 0 ),
	m_normSquared( 0 ),
	m_positiveKappa( 0 ),
	m_negativeKappa( 0 ),
	m_maximumNormSquared( 0 ),
	m_totalPositiveKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_totalNegativeKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_totalNormSquared( std::numeric_limits< double >::quiet_NaN() ),
	m_totalAverageLoss( std::numeric_limits< double >::quiet_NaN() )
{
}


SBP::SBP(
	boost::shared_ptr< Kernel::Base > pKernel,
	std::vector< std::string > const& optimizerParameters
) :
	m_pKernel( pKernel ),
	m_nu( 0 ),
	m_iterations( 0 ),
	m_normSquared( 0 ),
	m_positiveKappa( 0 ),
	m_negativeKappa( 0 ),
	m_alphas( new double[ m_pKernel->TrainingSize() ] ),
	m_responses( new double[ m_pKernel->Size() ] ),
	m_totalAlphas( m_pKernel->TrainingSize() ),
	m_totalResponses( m_pKernel->Size() ),
	m_normsSquared( new double[ m_pKernel->TrainingSize() ] ),
	m_maximumNormSquared( 0 ),
	m_totalPositiveKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_totalNegativeKappa( std::numeric_limits< double >::quiet_NaN() ),
	m_totalNormSquared( std::numeric_limits< double >::quiet_NaN() ),
	m_totalAverageLoss( std::numeric_limits< double >::quiet_NaN() )
{
	
	std::cout << "I am Biased SBP" << std::endl;

	if ( optimizerParameters.size() != 1 ) {

		std::stringstream stream;
		stream << "Biased::SBP optimizer expects one parameter (nu), received [";
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
	m_nu = atof( optimizerParameters[ 0 ].c_str() );
	if ( m_nu < 0 )
		throw std::runtime_error( "Biased::SBP optimizer parameter nu must be nonnegative" );

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	std::fill( m_alphas.get(),    m_alphas.get()    + trainingSize, 0 );
	std::fill( m_responses.get(), m_responses.get() + totalSize,    0 );

	for ( unsigned int ii = 0; ii < trainingSize; ++ii ) {

		double const normSquared = m_pKernel->KernelInnerProduct( ii, ii );
		if ( normSquared > m_maximumNormSquared )
			m_maximumNormSquared = normSquared;
		m_normsSquared[ ii ] = normSquared;
	}
}


double const SBP::AverageLoss() const {

	FindTotalLosses();

	return m_totalAverageLoss;
}


double const SBP::Objective() const {

	FindTotalKappas();

	return( 0.5 * ( m_totalPositiveKappa + m_totalNegativeKappa ) / m_iterations );
}


template< typename t_VectorType >
double const SBP::EvaluateHelper( t_VectorType const& vector ) const {

	FindTotalKappas();

	double const kappa =  0.5 * ( m_totalPositiveKappa + m_totalNegativeKappa );
	double const bias  = -0.5 * ( m_totalPositiveKappa - m_totalNegativeKappa );

	return( ( m_pKernel->Evaluate( vector, m_totalAlphas.Get() ) + bias ) / kappa );
}


template< typename t_Archive >
void SBP::save( t_Archive& archive, unsigned int const ) const {

	archive & boost::serialization::base_object< Base >( *this );

	archive & m_pKernel;
	archive & m_nu;
	archive & m_iterations;
	archive & m_normSquared;
	archive & m_positiveKappa;
	archive & m_negativeKappa;

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_alphas[ ii ];

	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & m_responses[ ii ];

	double const* totalAlphas = m_totalAlphas.Get();
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & totalAlphas[ ii ];

	double const* totalResponses = m_totalResponses.Get();
	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & totalResponses[ ii ];
}


template< typename t_Archive >
void SBP::load( t_Archive& archive, unsigned int const ) {

	archive & boost::serialization::base_object< Base >( *this );

	archive & m_pKernel;
	archive & m_nu;
	archive & m_iterations;
	archive & m_normSquared;
	archive & m_positiveKappa;
	archive & m_negativeKappa;

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	m_alphas.reset( new double[ trainingSize ] );
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_alphas[ ii ];

	m_responses.reset( new double[ totalSize ] );
	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & m_responses[ ii ];

	m_totalAlphas.Reset( trainingSize );
	double* totalAlphas = m_totalAlphas.Get();
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & totalAlphas[ ii ];

	m_totalResponses.Reset( totalSize );
	double* totalResponses = m_totalResponses.Get();
	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & totalResponses[ ii ];

	m_normsSquared.reset( new double[ trainingSize ] );
	m_maximumNormSquared = 0;
	for ( unsigned int ii = 0; ii < trainingSize; ++ii ) {

		double const normSquared = m_pKernel->KernelInnerProduct( ii, ii );
		if ( normSquared > m_maximumNormSquared )
			m_maximumNormSquared = normSquared;
		m_normsSquared[ ii ] = normSquared;
	}

	m_totalPositiveKappa = std::numeric_limits< double >::quiet_NaN();
	m_totalNegativeKappa = std::numeric_limits< double >::quiet_NaN();

	m_totalNormSquared = std::numeric_limits< double >::quiet_NaN();
	m_totalAverageLoss = std::numeric_limits< double >::quiet_NaN();
}




}    // namespace Biased


}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_OPTIMIZER_CLASSIFICATION_BIASED_SBP_HPP__ */
