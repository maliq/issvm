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
	\file array_sum.hpp
	\brief Implementation of ArraySum class
*/




#ifndef __SVM_ARRAY_SUM_HPP__
#define __SVM_ARRAY_SUM_HPP__

#ifdef __cplusplus




#include "helpers.h"

#include <boost/scoped_array.hpp>
#include <boost/shared_array.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <list>




//============================================================================
//    ArraySum class
//============================================================================


template< typename t_Type >
struct ArraySum {

	inline ArraySum();

	inline ArraySum( unsigned int const size );


	inline void Clear();

	inline void Reset( unsigned int const size );


	inline unsigned int const Size() const;

	inline t_Type const* const Get() const;
	inline t_Type* const Get();



	template< typename t_OtherType >
	void Add( t_OtherType const* other );


private:

	void Collapse();


	enum { TERMS = 256 };


	typedef std::list< std::pair< unsigned int, boost::shared_array< t_Type > > > List;

	unsigned int m_size;
	boost::scoped_array< t_Type > m_sum;
	List m_list;
};




//============================================================================
//    ArraySum inline methods
//============================================================================


template< typename t_Type >
ArraySum< t_Type >::ArraySum() : m_size( 0 ) {
}


template< typename t_Type >
ArraySum< t_Type >::ArraySum( unsigned int const size ) :
	m_size( size ),
	m_sum( new t_Type[ m_size ] )
{
	std::fill( m_sum.get(), m_sum.get() + m_size, 0 );
}


template< typename t_Type >
void ArraySum< t_Type >::Clear() {

	m_size = 0;
	m_sum.reset( NULL );
	m_list.clear();
}


template< typename t_Type >
void ArraySum< t_Type >::Reset( unsigned int const size ) {

	m_size = size;
	m_sum.reset( new t_Type[ m_size ] );
	std::fill( m_sum.get(), m_sum.get() + m_size, 0 );
	m_list.clear();
}


template< typename t_Type >
unsigned int const ArraySum< t_Type >::Size() const {

	return m_size;
}


template< typename t_Type >
t_Type const* const ArraySum< t_Type >::Get() const {

	BOOST_ASSERT( m_size > 0 );

	const_cast< ArraySum* >( this )->Collapse();
	return m_sum.get();
}


template< typename t_Type >
t_Type* const ArraySum< t_Type >::Get() {

	BOOST_ASSERT( m_size > 0 );

	Collapse();
	return m_sum.get();
}




//============================================================================
//    ArraySum methods
//============================================================================


template< typename t_Type >
template< typename t_OtherType >
void ArraySum< t_Type >::Add( t_OtherType const* other ) {

	BOOST_ASSERT( m_size > 0 );

	if ( UNLIKELY( m_list.empty() ) ) {

		boost::shared_array< t_Type > array( new t_Type[ m_size ] );
		std::copy( other, other + m_size, array.get() );
		m_list.push_back( std::pair< unsigned int, boost::shared_array< t_Type > >( 1, array ) );
	}
	else {

		typename List::iterator ii    = m_list.begin();
		typename List::iterator iiEnd = m_list.end();
		BOOST_ASSERT( ii != iiEnd );

		{	t_OtherType const* pSource    = other;
			t_OtherType const* pSourceEnd = pSource + m_size;
			t_Type* pDestination = ii->second.get();
			for ( ; pSource != pSourceEnd; ++pSource, ++pDestination )
				*pDestination += *pSource;
			++ii->first;
		}

		while ( ii->first >= TERMS ) {

			BOOST_ASSERT( ii != iiEnd );

			typename List::iterator const iiPrevious = ii;
			++ii;
			if ( UNLIKELY( ii == iiEnd ) ) {

				boost::shared_array< t_Type > array( new t_Type[ m_size ] );
				std::fill( array.get(), array.get() + m_size, 0 );
				ii = m_list.insert( ii, std::pair< unsigned int, boost::shared_array< t_Type > >( 0, array ) );
			}
			BOOST_ASSERT( ii != iiEnd );

			t_Type* pSource    = iiPrevious->second.get();
			t_Type* pSourceEnd = pSource + m_size;
			t_Type* pDestination = ii->second.get();
			for ( ; pSource != pSourceEnd; ++pSource, ++pDestination ) {

				*pDestination += *pSource;
				*pSource = 0;
			}
			iiPrevious->first = 0;
			++ii->first;
		}
	}
}


template< typename t_Type >
void ArraySum< t_Type >::Collapse() {

	if ( ! m_list.empty() ) {

		BOOST_ASSERT( m_size > 0 );

		typename List::iterator ii    = m_list.begin();
		typename List::iterator iiEnd = m_list.end();
		for ( ; ; ) {

			typename List::iterator const iiPrevious = ii;
			++ii;
			if ( UNLIKELY( ii == iiEnd ) )
				break;

			t_Type const* pSource    = iiPrevious->second.get();
			t_Type const* pSourceEnd = pSource + m_size;
			t_Type* pDestination = ii->second.get();
			for ( ; pSource != pSourceEnd; ++pSource, ++pDestination )
				*pDestination += *pSource;
		}

		t_Type const* pSource    = m_list.back().second.get();
		t_Type const* pSourceEnd = pSource + m_size;
		t_Type* pDestination = m_sum.get();
		for ( ; pSource != pSourceEnd; ++pSource, ++pDestination )
			*pDestination += *pSource;

		m_list.clear();
	}
}




#endif    /* __cplusplus */

#endif    /* __SVM_ARRAY_SUM_HPP__ */
