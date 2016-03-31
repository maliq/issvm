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
	\file svm_kernel_private_data.hpp
	\brief SVM::Kernel::_Private::Data implementation
*/




#ifndef __SVM_KERNEL_PRIVATE_DATA_HPP__
#define __SVM_KERNEL_PRIVATE_DATA_HPP__

#ifdef __cplusplus




#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/scoped_array.hpp>
#include <boost/assert.hpp>

#include <vector>




namespace SVM {


namespace Kernel {


namespace _Private {




//============================================================================
//    Data class
//============================================================================


template< typename t_VectorType >
struct Data {

	inline Data();

	template< typename t_OriginalVectorType >
	inline Data(
		std::vector< t_OriginalVectorType >& vectors,
		std::vector< double >& labels
	);


protected:

	unsigned int m_size;

	boost::scoped_array< t_VectorType > m_vectors;
	boost::scoped_array< double > m_labels;


private:

	template< typename t_Archive >
	inline void save( t_Archive& archive, unsigned int const ) const;

	template< typename t_Archive >
	inline void load( t_Archive& archive, unsigned int const );

	BOOST_SERIALIZATION_SPLIT_MEMBER()


	friend class boost::serialization::access;
};




//============================================================================
//    Data inline methods
//============================================================================


template< typename t_VectorType >
Data< t_VectorType >::Data() :
	m_size( 0 )
{
}


template< typename t_VectorType >
template< typename t_OriginalVectorType >
Data< t_VectorType >::Data(
	std::vector< t_OriginalVectorType >& vectors,
	std::vector< double >& labels
) :
	m_size( vectors.size() ),
	m_vectors( new t_VectorType[ m_size ] ),
	m_labels( new double[ m_size ] )
{
	BOOST_ASSERT( m_size == labels.size() );

	for ( unsigned int ii = 0; ii < m_size; ++ii ) {

		m_vectors[ ii ] = vectors[ ii ];
		vectors[ ii ].Clear();
		m_labels[ ii ] = labels[ ii ];
	}
	vectors.clear();
	labels.clear();
}


template< typename t_VectorType >
template< typename t_Archive >
void Data< t_VectorType >::save( t_Archive& archive, unsigned int const ) const {

	archive & m_size;

	for ( unsigned int ii = 0; ii < m_size; ++ii )
		archive & m_vectors[ ii ];

	for ( unsigned int ii = 0; ii < m_size; ++ii )
		archive & m_labels[ ii ];
}


template< typename t_VectorType >
template< typename t_Archive >
void Data< t_VectorType >::load( t_Archive& archive, unsigned int const ) {

	archive & m_size;

	m_vectors.reset( new t_VectorType[ m_size ] );
	for ( unsigned int ii = 0; ii < m_size; ++ii )
		archive & m_vectors[ ii ];

	m_labels.reset( new double[ m_size ] );
	for ( unsigned int ii = 0; ii < m_size; ++ii )
		archive & m_labels[ ii ];
}




}    // namespace _Private


}    // namespace Kernel


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_KERNEL_PRIVATE_DATA_HPP__ */
