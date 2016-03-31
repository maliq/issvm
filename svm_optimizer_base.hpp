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
	\file svm_optimizer_base.hpp
	\brief SVM::Optimizer::Base implementation
*/




#ifndef __SVM_OPTIMIZER_BASE_HPP__
#define __SVM_OPTIMIZER_BASE_HPP__

#ifdef __cplusplus




#include "svm_optimizer_base.hpp"
#include "random.hpp"
#include "vector.hpp"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/shared_array.hpp>




namespace SVM {


namespace Optimizer {




//============================================================================
//    Base class
//============================================================================


struct Base {

	inline Base();

	virtual ~Base() = 0;


	virtual unsigned int const TrainingSize()   const = 0;
	virtual unsigned int const ValidationSize() const = 0;

	virtual void GetAlphas( double* const begin, double* const end ) const = 0;

	virtual void WriteSupport(std::string filename);

	virtual double const Bias() const;

	virtual unsigned int const Support() const;
	virtual double const NormSquared() const = 0;

	virtual void GetValidationResponses( double* const begin, double* const end ) const = 0;
	virtual double const ValidationError() const;


	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< float  > const& otherVector ) const = 0;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SparseVector< double > const& otherVector ) const = 0;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector<   float  > const& otherVector ) const = 0;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, SpanVector<   double > const& otherVector ) const = 0;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector<  float  > const& otherVector ) const = 0;
	virtual double const Evaluate( Random::Generator::LaggedFibonacci4<>& generator, DenseVector<  double > const& otherVector ) const = 0;


	virtual void Iterate( Random::Generator::LaggedFibonacci4<>& generator ) = 0;

	virtual void Recalculate() = 0;


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




}    // namespace Optimizer


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_OPTIMIZER_BASE_HPP__ */
