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
	\file svm_serialization.hpp
	\brief Exports all SVM classes that require serialization
*/




#ifndef __SVM_SERIALIZATION_HPP__
#define __SVM_SERIALIZATION_HPP__

#ifdef __cplusplus




#include "svm.hpp"

#include <boost/serialization/export.hpp>
#include <boost/serialization/shared_ptr.hpp>




//============================================================================
//    Base serialization registration
//============================================================================


BOOST_CLASS_EXPORT_GUID( SVM::Kernel::Base,    "svm_kernel_base"    )
BOOST_CLASS_EXPORT_GUID( SVM::Optimizer::Base, "svm_optimizer_base" )


BOOST_SERIALIZATION_SHARED_PTR( SVM::Kernel::Base    )
BOOST_SERIALIZATION_SHARED_PTR( SVM::Optimizer::Base )




//============================================================================
//    Linear serialization registration
//============================================================================


BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::Simple< SparseVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ), "svm_kernel_sparse_float_linear"  )
BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::Simple< SparseVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ), "svm_kernel_sparse_double_linear" )

BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::Simple< SpanVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ), "svm_kernel_span_float_linear"  )
BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::Simple< SpanVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ), "svm_kernel_span_double_linear" )

BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::Simple< DenseVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ), "svm_kernel_dense_float_linear"  )
BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::Simple< DenseVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ), "svm_kernel_dense_double_linear" )


BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::Simple< SparseVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ) )
BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::Simple< SparseVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ) )

BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::Simple< SpanVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ) )
BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::Simple< SpanVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ) )

BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::Simple< DenseVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ) )
BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::Simple< DenseVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Linear > ) )




//============================================================================
//    Gaussian serialization registration
//============================================================================


BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::VectorData< SparseVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ), "svm_kernel_sparse_float_gaussian"  )
BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::VectorData< SparseVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ), "svm_kernel_sparse_double_gaussian" )

BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::VectorData< SpanVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ), "svm_kernel_span_float_gaussian"  )
BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::VectorData< SpanVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ), "svm_kernel_span_double_gaussian" )

BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::VectorData< DenseVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ), "svm_kernel_dense_float_gaussian"  )
BOOST_CLASS_EXPORT_GUID( PAREN_TYPE( SVM::Kernel::VectorData< DenseVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ), "svm_kernel_dense_double_gaussian" )


BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::VectorData< SparseVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ) )
BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::VectorData< SparseVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ) )

BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::VectorData< SpanVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ) )
BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::VectorData< SpanVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ) )

BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::VectorData< DenseVector< float  > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ) )
BOOST_SERIALIZATION_SHARED_PTR( PAREN_TYPE( SVM::Kernel::VectorData< DenseVector< double > BOOST_PP_COMMA() SVM::Kernel::Traits::Gaussian > ) )




//============================================================================
//    Optimizer registration
//============================================================================


BOOST_CLASS_EXPORT_GUID( SVM::Optimizer::Classification::Unbiased::SMO,        "svm_optimizer_classification_unbiased_smo"        )
BOOST_CLASS_EXPORT_GUID( SVM::Optimizer::Classification::Unbiased::SBP,        "svm_optimizer_classification_unbiased_sbp"        )
BOOST_CLASS_EXPORT_GUID( SVM::Optimizer::Classification::Unbiased::Perceptron, "svm_optimizer_classification_unbiased_perceptron" )
BOOST_CLASS_EXPORT_GUID( SVM::Optimizer::Classification::Unbiased::Sparsifier, "svm_optimizer_classification_unbiased_sparsifier" )

BOOST_CLASS_EXPORT_GUID( SVM::Optimizer::Classification::Biased::SMO,        "svm_optimizer_classification_biased_smo"        )
BOOST_CLASS_EXPORT_GUID( SVM::Optimizer::Classification::Biased::SBP,        "svm_optimizer_classification_biased_sbp"        )
BOOST_CLASS_EXPORT_GUID( SVM::Optimizer::Classification::Biased::Perceptron, "svm_optimizer_classification_biased_perceptron" )
BOOST_CLASS_EXPORT_GUID( SVM::Optimizer::Classification::Biased::Sparsifier, "svm_optimizer_classification_biased_sparsifier" )


BOOST_SERIALIZATION_SHARED_PTR( SVM::Optimizer::Classification::Unbiased::SMO        )
BOOST_SERIALIZATION_SHARED_PTR( SVM::Optimizer::Classification::Unbiased::SBP        )
BOOST_SERIALIZATION_SHARED_PTR( SVM::Optimizer::Classification::Unbiased::Perceptron )
BOOST_SERIALIZATION_SHARED_PTR( SVM::Optimizer::Classification::Unbiased::Sparsifier )

BOOST_SERIALIZATION_SHARED_PTR( SVM::Optimizer::Classification::Biased::SMO        )
BOOST_SERIALIZATION_SHARED_PTR( SVM::Optimizer::Classification::Biased::SBP        )
BOOST_SERIALIZATION_SHARED_PTR( SVM::Optimizer::Classification::Biased::Perceptron )
BOOST_SERIALIZATION_SHARED_PTR( SVM::Optimizer::Classification::Biased::Sparsifier )




#endif    /* __cplusplus */

#endif    /* __SVM_SERIALIZATION_HPP__ */
