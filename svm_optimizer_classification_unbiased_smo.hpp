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
	\file svm_optimizer_classification_unbiased_smo.hpp
	\brief SVM::Optimizer::Unbiased::SMO implementation
*/




#ifndef __SVM_OPTIMIZER_CLASSIFICATION_UNBIASED_SMO_HPP__
#define __SVM_OPTIMIZER_CLASSIFICATION_UNBIASED_SMO_HPP__

#ifdef __cplusplus




#include "svm_optimizer_base.hpp"
#include "svm_kernel.hpp"
#include "vector.hpp"

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

#include <iostream>
#include <fstream>


namespace SVM {


namespace Optimizer {


namespace Classification {


namespace Unbiased {




//============================================================================
//    SMO class
//============================================================================


struct SMO : public Base {

	inline SMO();

	inline SMO(
		boost::shared_ptr< Kernel::Base > pKernel,
		std::vector< std::string > const& optimizerParameters
	);

	virtual ~SMO();


	virtual unsigned int const TrainingSize()   const;
	virtual unsigned int const ValidationSize() const;

	virtual void GetAlphas( double* const begin, double* const end ) const;


	virtual void WriteSupport(std::string filename){

		std::cout << "Get Support from Unbiased SMO" << std::endl;
		m_pKernel->writeSupport(m_alphas.get(),0.0,filename);

	}

	virtual double const NormSquared() const;

	virtual void GetValidationResponses( double* const begin, double* const end ) const;

	inline double const AverageLoss() const;
	inline double const AlphaSum() const;

	inline double const PrimalObjective() const;
	inline double const DualObjective() const;


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

	void FindLosses() const;


	template< typename t_Archive >
	inline void save( t_Archive& archive, unsigned int const ) const;

	template< typename t_Archive >
	inline void load( t_Archive& archive, unsigned int const );

	BOOST_SERIALIZATION_SPLIT_MEMBER()


	boost::shared_ptr< Kernel::Base > m_pKernel;

	double m_lambda;

	boost::scoped_array< double > m_alphas;
	boost::scoped_array< double > m_responses;

	boost::scoped_array< double > m_normsSquared;

	double mutable m_normSquared;
	double mutable m_averageLoss;
	double mutable m_alphaSum;


	friend class boost::serialization::access;
};




//============================================================================
//    SMO inline methods
//============================================================================


SMO::SMO() :
	m_lambda( 0 ),
	m_normSquared( std::numeric_limits< double >::quiet_NaN() ),
	m_averageLoss( std::numeric_limits< double >::quiet_NaN() ),
	m_alphaSum(    std::numeric_limits< double >::quiet_NaN() )
{
}


SMO::SMO(
	boost::shared_ptr< Kernel::Base > pKernel,
	std::vector< std::string > const& optimizerParameters
) :
	m_pKernel( pKernel ),
	m_lambda( 0 ),
	m_alphas( new double[ m_pKernel->TrainingSize() ] ),
	m_responses( new double[ m_pKernel->Size() ] ),
	m_normsSquared( new double[ m_pKernel->TrainingSize() ] ),
	m_normSquared( std::numeric_limits< double >::quiet_NaN() ),
	m_averageLoss( std::numeric_limits< double >::quiet_NaN() ),
	m_alphaSum(    std::numeric_limits< double >::quiet_NaN() )
{
	
	std::cout << "I am Unbiased SMO" << std::endl;

	if ( optimizerParameters.size() != 1 ) {

		std::stringstream stream;
		stream << "Unbiased::SMO optimizer expects one parameter (lambda), received [";
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
	m_lambda = atof( optimizerParameters[ 0 ].c_str() );
	if ( m_lambda < 0 )
		throw std::runtime_error( "Unbiased::SMO optimizer parameter lambda must be nonnegative" );

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	std::fill( m_alphas.get(),    m_alphas.get()    + trainingSize, 0 );
	std::fill( m_responses.get(), m_responses.get() + totalSize,    0 );

	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		m_normsSquared[ ii ] = m_pKernel->KernelInnerProduct( ii, ii );
}


double const SMO::AverageLoss() const {

	FindLosses();

	return m_averageLoss;
}


double const SMO::AlphaSum() const {

	FindLosses();

	return m_alphaSum;
}


double const SMO::PrimalObjective() const {

	FindLosses();

	return( 0.5 * m_lambda * m_normSquared + m_averageLoss );
}


double const SMO::DualObjective() const {

	FindLosses();

	return( m_lambda * ( m_alphaSum - 0.5 * m_normSquared ) );
}


template< typename t_VectorType >
double const SMO::EvaluateHelper( t_VectorType const& vector ) const {

	return m_pKernel->Evaluate( vector, m_alphas.get() );
}


template< typename t_Archive >
void SMO::save( t_Archive& archive, unsigned int const ) const {

	archive & boost::serialization::base_object< Base >( *this );

	archive & m_pKernel;
	archive & m_lambda;

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_alphas[ ii ];

	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & m_responses[ ii ];
}


template< typename t_Archive >
void SMO::load( t_Archive& archive, unsigned int const ) {

	archive & boost::serialization::base_object< Base >( *this );

	archive & m_pKernel;
	archive & m_lambda;

	unsigned int const trainingSize = m_pKernel->TrainingSize();
	unsigned int const totalSize    = m_pKernel->Size();

	m_alphas.reset( new double[ trainingSize ] );
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		archive & m_alphas[ ii ];

	m_responses.reset( new double[ totalSize ] );
	for ( unsigned int ii = 0; ii < totalSize; ++ii )
		archive & m_responses[ ii ];

	m_normsSquared.reset( new double[ trainingSize ] );
	for ( unsigned int ii = 0; ii < trainingSize; ++ii )
		m_normsSquared[ ii ] = m_pKernel->KernelInnerProduct( ii, ii );

	m_normSquared = std::numeric_limits< double >::quiet_NaN();
	m_averageLoss = std::numeric_limits< double >::quiet_NaN();
	m_alphaSum    = std::numeric_limits< double >::quiet_NaN();
}




}    // namespace Unbiased


}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_OPTIMIZER_CLASSIFICATION_UNBIASED_SMO_HPP__ */
