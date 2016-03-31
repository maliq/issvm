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
	\file svm_kernel_traits_linear.hpp
	\brief SVM::Kernel::Traits::Linear implementation
*/




#ifndef __SVM_KERNEL_TRAITS_LINEAR_HPP__
#define __SVM_KERNEL_TRAITS_LINEAR_HPP__

#ifdef __cplusplus




#include "svm_kernel_simple.hpp"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <sstream>
#include <vector>
#include <stdexcept>




namespace SVM {


namespace Kernel {


namespace Traits {




//============================================================================
//    Linear traits class
//============================================================================


struct Linear {

	// allowed vector types
	enum { SPARSE_ALLOWED = true };
	enum { SPAN_ALLOWED   = true };


	inline Linear();

	inline Linear( std::vector< std::string > const& kernelParameters );


protected:

	// (allowed,allowed) or (dense,allowed)
	template< typename t_Vector1Type, typename t_Vector2Type >
	double const KernelInnerProduct( t_Vector1Type const& vector1, t_Vector2Type const& vector2 ) const;


private:

	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	friend class boost::serialization::access;
};




//============================================================================
//    Linear inline methods
//============================================================================


Linear::Linear() {
}


Linear::Linear( std::vector< std::string > const& kernelParameters ) {

	if ( kernelParameters.size() != 0 ) {

		std::stringstream stream;
		stream << "Linear kernel expects no parameters, received [";
		{	std::vector< std::string >::const_iterator ii    = kernelParameters.begin();
			std::vector< std::string >::const_iterator iiEnd = kernelParameters.end();
			if ( ii != iiEnd ) {

				stream << *ii;
				for ( ++ii; ii != iiEnd; ++ii )
					stream << ',' << *ii;
			}
		}
		stream << ']';
		throw std::runtime_error( stream.str() );
	}
}


template< typename t_Vector1Type, typename t_Vector2Type >
double const Linear::KernelInnerProduct( t_Vector1Type const& vector1, t_Vector2Type const& vector2 ) const {

	return vector1.InnerProduct( vector2 );
}


template< typename t_Archive >
inline void Linear::serialize( t_Archive& archive, unsigned int const ) {
}




}    // namespace Traits


}    // namespace Kernel


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_KERNEL_TRAITS_LINEAR_HPP__ */
