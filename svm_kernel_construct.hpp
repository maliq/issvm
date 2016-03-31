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
	\file svm_kernel_construct.hpp
	\brief SVM::Kernel::Construct function
*/




#ifndef __SVM_KERNEL_CONSTRUCT_HPP__
#define __SVM_KERNEL_CONSTRUCT_HPP__

#ifdef __cplusplus




#include "svm_kernel_base.hpp"
#include "vector.hpp"

#include <boost/shared_ptr.hpp>

#include <vector>




namespace SVM {


namespace Kernel {




//============================================================================
//    Construct function
//============================================================================


template< typename t_Type, template< typename, typename > class t_KernelType, typename t_Traits, typename t_OriginalVectorType >
inline boost::shared_ptr< Base > Construct(
	std::vector< t_OriginalVectorType >& vectors,
	std::vector< double >& labels,
	unsigned int const trainingSize,
	std::vector< std::string > const& kernelParameters,
	unsigned int const cacheSize = 0
)
{
	uint64_t sparseSize = 0;
	if ( t_Traits::SPARSE_ALLOWED ) {

		typename std::vector< SparseVector< t_Type > >::const_iterator ii    = vectors.begin();
		typename std::vector< SparseVector< t_Type > >::const_iterator iiEnd = vectors.end();
		for ( ; ii != iiEnd; ++ii )
			sparseSize += ii->Size();
	}

	uint64_t spanSize = 0;
	if ( t_Traits::SPAN_ALLOWED ) {

		typename std::vector< SparseVector< t_Type > >::const_iterator ii    = vectors.begin();
		typename std::vector< SparseVector< t_Type > >::const_iterator iiEnd = vectors.end();
		for ( ; ii != iiEnd; ++ii )
			spanSize += SpanVector< t_Type >( *ii ).Size();
	}

	uint64_t denseSize = 0;
	{	typename std::vector< SparseVector< t_Type > >::const_iterator ii    = vectors.begin();
		typename std::vector< SparseVector< t_Type > >::const_iterator iiEnd = vectors.end();
		for ( ; ii != iiEnd; ++ii )
			denseSize += DenseVector< t_Type >( *ii ).Size();
	}

	boost::shared_ptr< Base > result;

	if (
		t_Traits::SPARSE_ALLOWED &&
		( ( ! t_Traits::SPAN_ALLOWED ) || ( sparseSize <= spanSize ) ) &&
		( sparseSize <= denseSize )
	)
	{
		result = boost::shared_ptr< Base >( new t_KernelType< SparseVector< t_Type >, t_Traits >( vectors, labels, trainingSize, kernelParameters, cacheSize ) );
	}
	else if (
		t_Traits::SPAN_ALLOWED &&
		( spanSize <= denseSize )
	)
	{
		result = boost::shared_ptr< Base >( new t_KernelType< SpanVector< t_Type >, t_Traits >( vectors, labels, trainingSize, kernelParameters, cacheSize ) );
	}
	else {

		result = boost::shared_ptr< Base >( new t_KernelType< DenseVector< t_Type >, t_Traits >( vectors, labels, trainingSize, kernelParameters, cacheSize ) );
	}

	return result;
}




}    // namespace Kernel


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_KERNEL_CONSTRUCT_HPP__ */
