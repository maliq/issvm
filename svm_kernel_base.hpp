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
	\file svm_kernel_base.hpp
	\brief SVM::Kernel::Base implementation
*/




#ifndef __SVM_KERNEL_BASE_HPP__
#define __SVM_KERNEL_BASE_HPP__

#ifdef __cplusplus




#include "vector.hpp"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/shared_array.hpp>




namespace SVM {


namespace Kernel {




//============================================================================
//    Base class
//============================================================================


struct Base {

	inline Base();

	virtual ~Base();


	// Size() is the total number of training+validation elements
	virtual unsigned int const Size() const = 0;
	virtual unsigned int const TrainingSize() const = 0;

	virtual double const* const Labels() const = 0;


	// index1,index2 < Size()
	virtual double const KernelInnerProduct(
		unsigned int const index1,
		unsigned int const index2
	) const = 0;

	// index < TrainingSize()
	virtual boost::shared_array< double > const Row(
		unsigned int const index
	) = 0;

	// alphas[ TrainingSize() ], responses[ Size() ], index < TrainingSize()
	virtual void SetAlpha(
		double* const alphas,
		double* const responses,
		unsigned int const index,
		double const alpha
	) const = 0;

	// alphas[ TrainingSize() ], responses[ Size() ], pRow[ Size() ], index < TrainingSize()
	virtual void SetAlpha(
		double* const alphas,
		double* const responses,
		unsigned int const index,
		double const alpha,
		boost::shared_array< double > const& pRow
	) const = 0;

	// alphas[ TrainingSize() ], responses[ Size() ]
	virtual void RecalculateResponses(
		double const* const alphas,
		double* const responses
	) const = 0;


	// alphas[ TrainingSize() ]
	virtual double const Evaluate(
		SparseVector< float > const& otherVector,
		double const* const alphas
	) const = 0;

	// alphas[ TrainingSize() ]
	virtual double const Evaluate(
		SparseVector< double > const& otherVector,
		double const* const alphas
	) const = 0;

	// alphas[ TrainingSize() ]
	virtual double const Evaluate(
		SpanVector< float > const& otherVector,
		double const* const alphas
	) const = 0;

	// alphas[ TrainingSize() ]
	virtual double const Evaluate(
		SpanVector< double > const& otherVector,
		double const* const alphas
	) const = 0;

	// alphas[ TrainingSize() ]
	virtual double const Evaluate(
		DenseVector< float > const& otherVector,
		double const* const alphas
	) const = 0;

	// alphas[ TrainingSize() ]
	virtual double const Evaluate(
		DenseVector< double > const& otherVector,
		double const* const alphas
	) const = 0;

	virtual void const writeSupport(
		double const* const alphas, double bias, std::string filename
	) const = 0;

private:

	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	friend class boost::serialization::access;
};




//============================================================================
//    Base inline methods
//============================================================================


Base::Base() {
}


template< typename t_Archive >
void Base::serialize( t_Archive& archive, unsigned int const ) {
}




}    // namespace Kernel


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_KERNEL_BASE_HPP__ */
