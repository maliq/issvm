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
	\file svm_kernel_traits_gaussian.hpp
	\brief SVM::Kernel::Traits::Gaussian implementation
*/




#ifndef __SVM_KERNEL_TRAITS_GAUSSIAN_HPP__
#define __SVM_KERNEL_TRAITS_GAUSSIAN_HPP__

#ifdef __cplusplus




#include "svm_kernel_vector_data.hpp"

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
//    Gaussian triats class
//============================================================================


struct Gaussian {

	// allowed vector types
	enum { SPARSE_ALLOWED = true };
	enum { SPAN_ALLOWED   = true };


	inline Gaussian();

	inline Gaussian( std::vector< std::string > const& kernelParameters );


protected:

	typedef double VectorDataType;

	// parameter can be non-allowed vector types
	template< typename t_OtherVectorType >
	inline VectorDataType const VectorData( t_OtherVectorType const& vector ) const;


	// (allowed,allowed) or (dense,allowed)
	template< typename t_Vector1Type, typename t_Vector2Type >
	double const KernelInnerProduct( t_Vector1Type const& vector1, VectorDataType const& vectorData1, t_Vector2Type const& vector2, VectorDataType const& vectorData2 ) const;


private:

	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	double m_gamma;


	friend class boost::serialization::access;
};




//============================================================================
//    Gaussian inline methods
//============================================================================


Gaussian::Gaussian() :
	m_gamma( 0 )
{
}


Gaussian::Gaussian( std::vector< std::string > const& kernelParameters ) :
	m_gamma( 0 )
{
	if ( kernelParameters.size() != 1 ) {

		std::stringstream stream;
		stream << "Gaussian kernel expects one parameter (gamma), received [";
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
	m_gamma = atof( kernelParameters[ 0 ].c_str() );
	if ( m_gamma <= 0 )
		throw std::runtime_error( "Gaussian kernel parameter gamma must be positive" );
}


template< typename t_OtherVectorType >
typename Gaussian::VectorDataType const Gaussian::VectorData( t_OtherVectorType const& vector ) const {

	return vector.NormSquared();
}


template< typename t_Vector1Type, typename t_Vector2Type >
double const Gaussian::KernelInnerProduct( t_Vector1Type const& vector1, VectorDataType const& vectorData1, t_Vector2Type const& vector2, VectorDataType const& vectorData2 ) const {

	return std::exp( m_gamma * ( 2 * vector1.InnerProduct( vector2 ) - vectorData1 - vectorData2 ) );
}


template< typename t_Archive >
inline void Gaussian::serialize( t_Archive& archive, unsigned int const ) {

	archive & m_gamma;
}




}    // namespace Traits


}    // namespace Kernel


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_KERNEL_TRAITS_GAUSSIAN_HPP__ */
