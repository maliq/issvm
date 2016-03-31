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
	\file svm_kernel_private_cache.hpp
	\brief SVM::Kernel::_Private::Cache implementation
*/




#ifndef __SVM_KERNEL_PRIVATE_CACHE_HPP__
#define __SVM_KERNEL_PRIVATE_CACHE_HPP__

#ifdef __cplusplus




#include "helpers.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>
#include <boost/assert.hpp>

#include <list>
#include <map>




namespace SVM {


namespace Kernel {


namespace _Private {




//============================================================================
//    Cache class
//============================================================================


struct Cache {

	inline Cache();

	inline Cache( unsigned int const size, unsigned int const dimension );


	inline unsigned int const Size() const;


	std::pair< bool, boost::shared_array< double > > const operator[]( unsigned int const index );


private:

	template< typename t_Archive >
	inline void save( t_Archive& archive, unsigned int const ) const;

	template< typename t_Archive >
	inline void load( t_Archive& archive, unsigned int const );

	BOOST_SERIALIZATION_SPLIT_MEMBER()


	typedef std::pair< unsigned int, boost::shared_array< double > > Element;

	unsigned int m_size;
	unsigned int m_dimension;

	std::list< Element > m_cache;
	boost::scoped_array< std::list< Element >::iterator > m_table;


	friend class boost::serialization::access;
};




//============================================================================
//    Cache inline methods
//============================================================================


Cache::Cache() :
	m_size( 0 ),
	m_dimension( 0 )
{
}


Cache::Cache( unsigned int const size, unsigned int const dimension ) :
	m_size( size ),
	m_dimension( dimension ),
	m_table( new std::list< Element >::iterator[ m_dimension ] )
{
	BOOST_ASSERT( m_dimension > 0 );

	std::fill( m_table.get(), m_table.get() + m_dimension, m_cache.end() );
}


unsigned int const Cache::Size() const {

	return m_size;
}


template< typename t_Archive >
void Cache::save( t_Archive& archive, unsigned int const ) const {

	archive & m_size;
	archive & m_dimension;
}


template< typename t_Archive >
void Cache::load( t_Archive& archive, unsigned int const ) {

	archive & m_size;
	archive & m_dimension;

	BOOST_ASSERT( m_dimension > 0 );

	m_cache.clear();
	m_table.reset( new std::list< Element >::iterator[ m_dimension ] );
	std::fill( m_table.get(), m_table.get() + m_dimension, m_cache.end() );
}




}    // namespace _Private


}    // namespace Kernel


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_KERNEL_PRIVATE_CACHE_HPP__ */
