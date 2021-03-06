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
	\file svm_kernel_simple.hpp
	\brief SVM::Kernel::Simple implementation
*/




#ifndef __SVM_KERNEL_SIMPLE_HPP__
#define __SVM_KERNEL_SIMPLE_HPP__

#ifdef __cplusplus




#include "svm_kernel_private_data.hpp"
#include "svm_kernel_private_cache.hpp"
#include "svm_kernel_base.hpp"
#include "vector.hpp"
#include "helpers.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/shared_array.hpp>
#include <boost/assert.hpp>

#include <vector>


#include <iostream>
#include <fstream>

namespace SVM {


namespace Kernel {




//============================================================================
//    Simple class
//============================================================================


template< typename t_VectorType, typename t_Traits >
struct Simple : public Base, public _Private::Data< t_VectorType >, public t_Traits {

	inline Simple();

	template< typename t_OriginalVectorType >
	inline Simple(
		std::vector< t_OriginalVectorType >& vectors,
		std::vector< double >& labels,
		unsigned int const trainingSize,
		std::vector< std::string > const& kernelParameters,
		unsigned int const cacheSize = 0
	);

	virtual ~Simple();


	virtual unsigned int const Size() const;
	virtual unsigned int const TrainingSize() const;

	virtual double const* const Labels() const;


	virtual double const KernelInnerProduct(
		unsigned int const index1,
		unsigned int const index2
	) const;

	virtual boost::shared_array< double > const Row(
		unsigned int const index
	);

	virtual void SetAlpha(
		double* const alphas,
		double* const responses,
		unsigned int const index,
		double const alpha
	) const;

	virtual void SetAlpha(
		double* const alphas,
		double* const responses,
		unsigned int const index,
		double const alpha,
		boost::shared_array< double > const& row
	) const;

	virtual void RecalculateResponses(
		double const* const alphas,
		double* const responses
	) const;


	virtual double const Evaluate(
		SparseVector< float > const& otherVector,
		double const* const alphas
	) const;

	virtual double const Evaluate(
		SparseVector< double > const& otherVector,
		double const* const alphas
	) const;

	virtual double const Evaluate(
		SpanVector< float > const& otherVector,
		double const* const alphas
	) const;

	virtual double const Evaluate(
		SpanVector< double > const& otherVector,
		double const* const alphas
	) const;

	virtual double const Evaluate(
		DenseVector< float > const& otherVector,
		double const* const alphas
	) const;

	virtual double const Evaluate(
		DenseVector< double > const& otherVector,
		double const* const alphas
	) const;

	virtual void const writeSupport(
		double const* const alphas, double bias, std::string filename
	) const;

private:

	double const EvaluateHelper(
		DenseVector< typename t_VectorType::Type > const& vector1,
		double const* const alphas
	) const;


	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	unsigned int m_trainingSize;

	_Private::Cache mutable m_cache;


	friend class boost::serialization::access;
};




//============================================================================
//    Simple inline methods
//============================================================================


template< typename t_VectorType, typename t_Traits >
Simple< t_VectorType, t_Traits >::Simple() :
	Base(),
	_Private::Data< t_VectorType >(),
	t_Traits()
{

	std::cout << "I am Kernel Simple ..." << std::endl;
}


template< typename t_VectorType, typename t_Traits >
template< typename t_OriginalVectorType >
Simple< t_VectorType, t_Traits >::Simple(
	std::vector< t_OriginalVectorType >& vectors,
	std::vector< double >& labels,
	unsigned int const trainingSize,
	std::vector< std::string > const& kernelParameters,
	unsigned int const cacheSize
) :
	Base(),
	_Private::Data< t_VectorType >( vectors, labels ),
	t_Traits( kernelParameters ),
	m_trainingSize( trainingSize ),
	m_cache( cacheSize, _Private::Data< t_VectorType >::m_size )
{
	BOOST_ASSERT( trainingSize <= _Private::Data< t_VectorType >::m_size );
	std::cout << "I am Kernel Simple ... " << std::endl;
}


template< typename t_VectorType, typename t_Traits >
Simple< t_VectorType, t_Traits >::~Simple() {
}


template< typename t_VectorType, typename t_Traits >
unsigned int const Simple< t_VectorType, t_Traits >::Size() const {

	return _Private::Data< t_VectorType >::m_size;
}


template< typename t_VectorType, typename t_Traits >
unsigned int const Simple< t_VectorType, t_Traits >::TrainingSize() const {

	return m_trainingSize;
}


template< typename t_VectorType, typename t_Traits >
double const* const Simple< t_VectorType, t_Traits >::Labels() const {

	return _Private::Data< t_VectorType >::m_labels.get();
}


template< typename t_VectorType, typename t_Traits >
double const Simple< t_VectorType, t_Traits >::KernelInnerProduct(
	unsigned int const index1,
	unsigned int const index2
) const
{
	BOOST_ASSERT( index1 < _Private::Data< t_VectorType >::m_size );
	BOOST_ASSERT( index2 < _Private::Data< t_VectorType >::m_size );

	t_VectorType const& vector1( _Private::Data< t_VectorType >::m_vectors[ index1 ] );
	t_VectorType const& vector2( _Private::Data< t_VectorType >::m_vectors[ index2 ] );

	return t_Traits::KernelInnerProduct( vector1, vector2 );
}


template< typename t_VectorType, typename t_Traits >
boost::shared_array< double > const Simple< t_VectorType, t_Traits >::Row( unsigned int const index ) {

	assert( index < TrainingSize() );

	unsigned int const totalSize = _Private::Data< t_VectorType >::m_size;

	std::pair< bool, boost::shared_array< double > > cached = m_cache[ index ];

	if ( LIKELY( ! cached.first ) ) {

		DenseVector< typename t_VectorType::Type > vector1( _Private::Data< t_VectorType >::m_vectors[ index ] );

		double* const cache = cached.second.get();
		t_VectorType const* const trainingVectors = _Private::Data< t_VectorType >::m_vectors.get();
		#pragma omp parallel for schedule( static )
		for ( int ii = 0; ii < static_cast< int >( totalSize ); ++ii )
			cache[ ii ] = t_Traits::KernelInnerProduct( vector1, trainingVectors[ ii ] );
	}

	return cached.second;
}


template< typename t_VectorType, typename t_Traits >
void Simple< t_VectorType, t_Traits >::SetAlpha(
	double* const alphas,
	double* const responses,
	unsigned int const index,
	double const alpha
) const
{
	assert( index < TrainingSize() );

	unsigned int const totalSize = _Private::Data< t_VectorType >::m_size;

	double const deltaAlpha = alpha - alphas[ index ];
	if ( deltaAlpha != 0 ) {

		if ( m_cache.Size() == 0 ) {

			DenseVector< typename t_VectorType::Type > vector1( _Private::Data< t_VectorType >::m_vectors[ index ] );

			t_VectorType const* const trainingVectors = _Private::Data< t_VectorType >::m_vectors.get();
			#pragma omp parallel for schedule( static )
			for ( int ii = 0; ii < static_cast< int >( totalSize ); ++ii )
				responses[ ii ] += deltaAlpha * t_Traits::KernelInnerProduct( vector1, trainingVectors[ ii ] );
		}
		else {

			std::pair< bool, boost::shared_array< double > > cached = m_cache[ index ];

			if ( LIKELY( ! cached.first ) ) {

				DenseVector< typename t_VectorType::Type > vector1( _Private::Data< t_VectorType >::m_vectors[ index ] );

				double* const cache = cached.second.get();
				t_VectorType const* const trainingVectors = _Private::Data< t_VectorType >::m_vectors.get();
				#pragma omp parallel for schedule( static )
				for ( int ii = 0; ii < static_cast< int >( totalSize ); ++ii ) {

					double const kernelInnerProduct = t_Traits::KernelInnerProduct( vector1, trainingVectors[ ii ] );
					cache[ ii ] = kernelInnerProduct;
					responses[ ii ] += deltaAlpha * kernelInnerProduct;
				}
			}
			else {

				double* pResponse    = responses;
				double* pResponseEnd = pResponse + totalSize;
				double const* pCache = cached.second.get();
				for ( ; pResponse != pResponseEnd; ++pResponse, ++pCache ) {

					double const kernelInnerProduct = *pCache;
					*pResponse += deltaAlpha * kernelInnerProduct;
				}
			}
		}
	}

	alphas[ index ] = alpha;
}


template< typename t_VectorType, typename t_Traits >
void Simple< t_VectorType, t_Traits >::SetAlpha(
	double* const alphas,
	double* const responses,
	unsigned int const index,
	double const alpha,
	boost::shared_array< double > const& row
) const
{
	assert( index < TrainingSize() );

	unsigned int const totalSize = _Private::Data< t_VectorType >::m_size;

	double const deltaAlpha = alpha - alphas[ index ];
	if ( deltaAlpha != 0 ) {

		double* pResponse    = responses;
		double* pResponseEnd = pResponse + totalSize;
		double const* pCache = row.get();
		for ( ; pResponse != pResponseEnd; ++pResponse, ++pCache ) {

			double const kernelInnerProduct = *pCache;
			*pResponse += deltaAlpha * kernelInnerProduct;
		}
	}

	alphas[ index ] = alpha;
}


template< typename t_VectorType, typename t_Traits >
void Simple< t_VectorType, t_Traits >::RecalculateResponses(
	double const* const alphas,
	double* const responses
) const
{
	unsigned int const totalSize = _Private::Data< t_VectorType >::m_size;

	t_VectorType const* const trainingVectors = _Private::Data< t_VectorType >::m_vectors.get();
	#pragma omp parallel for schedule( static )
	for ( int ii = 0; ii < static_cast< int >( totalSize ); ++ii )
		responses[ ii ] = EvaluateHelper( trainingVectors[ ii ], alphas );
}


template< typename t_VectorType, typename t_Traits >
double const Simple< t_VectorType, t_Traits >::Evaluate(
	SparseVector< float > const& otherVector,
	double const* const alphas
) const
{
	return EvaluateHelper( otherVector, alphas );
}


template< typename t_VectorType, typename t_Traits >
double const Simple< t_VectorType, t_Traits >::Evaluate(
	SparseVector< double > const& otherVector,
	double const* const alphas
) const
{
	return EvaluateHelper( otherVector, alphas );
}


template< typename t_VectorType, typename t_Traits >
double const Simple< t_VectorType, t_Traits >::Evaluate(
	SpanVector< float > const& otherVector,
	double const* const alphas
) const
{
	return EvaluateHelper( otherVector, alphas );
}


template< typename t_VectorType, typename t_Traits >
double const Simple< t_VectorType, t_Traits >::Evaluate(
	SpanVector< double > const& otherVector,
	double const* const alphas
) const
{
	return EvaluateHelper( otherVector, alphas );
}


template< typename t_VectorType, typename t_Traits >
double const Simple< t_VectorType, t_Traits >::Evaluate(
	DenseVector< float > const& otherVector,
	double const* const alphas
) const
{
	return EvaluateHelper( otherVector, alphas );
}


template< typename t_VectorType, typename t_Traits >
double const Simple< t_VectorType, t_Traits >::Evaluate(
	DenseVector< double > const& otherVector,
	double const* const alphas
) const
{
	return EvaluateHelper( otherVector, alphas );
}

template< typename t_VectorType, typename t_Traits >
void const Simple< t_VectorType, t_Traits >::writeSupport(
	double const* const alphas, double bias, std::string filename
) const
{

	unsigned int const trainingSize = m_trainingSize;
	std::cout << "Write Support from Kernel Vector Data" << std::endl;
	std::cout << filename << std::endl;
	double accumulator = 0;

	double const* pAlpha    = alphas;
	double const* pAlphaEnd = pAlpha + trainingSize;
	t_VectorType const* pTrainingVector = _Private::Data< t_VectorType >::m_vectors.get();

	std::ofstream myfile;
  	myfile.open(filename);
  	
  	for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pTrainingVector) {

		double const alpha = *pAlpha;

		if ( alpha != 0 ){
			myfile << bias << " ";
			myfile << alpha << " ";
			myfile << *pTrainingVector << std::endl;
		}
	}

	myfile.close();

}

template< typename t_VectorType, typename t_Traits >
double const Simple< t_VectorType, t_Traits >::EvaluateHelper(
	DenseVector< typename t_VectorType::Type > const& vector1,
	double const* const alphas
) const
{
	unsigned int const trainingSize = m_trainingSize;

	double accumulator = 0;

	double const* pAlpha    = alphas;
	double const* pAlphaEnd = pAlpha + trainingSize;
	t_VectorType const* pTrainingVector = _Private::Data< t_VectorType >::m_vectors.get();
	for ( ; pAlpha != pAlphaEnd; ++pAlpha, ++pTrainingVector ) {

		double const alpha = *pAlpha;
		if ( alpha != 0 )
			accumulator += alpha * t_Traits::KernelInnerProduct( vector1, *pTrainingVector );
	}

	return accumulator;
}


template< typename t_VectorType, typename t_Traits >
template< typename t_Archive >
void Simple< t_VectorType, t_Traits >::serialize( t_Archive& archive, unsigned int const ) {

	archive & boost::serialization::base_object< Base >( *this );
	archive & boost::serialization::base_object< _Private::Data< t_VectorType > >( *this );
	archive & boost::serialization::base_object< t_Traits >( *this );
	archive & m_trainingSize;
	archive & m_cache;
}




}    // namespace Kernel


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_KERNEL_SIMPLE_HPP__ */
