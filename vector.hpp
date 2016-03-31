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
	\file vector.hpp
	\brief SparseVector, SpanVector and DenseVector implementations
*/




#ifndef __VECTOR_HPP__
#define __VECTOR_HPP__

#ifdef __cplusplus




#include "helpers.h"

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/type_traits.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>

#include <ostream>
#include <vector>
#include <algorithm>




//============================================================================
//    forward declarations
//============================================================================


template< typename t_Type >
struct SparseVector;


template< typename t_Type >
struct SpanVector;


template< typename t_Type >
struct DenseVector;




//============================================================================
//    SparseVector class
//============================================================================


template< typename t_Type >
struct SparseVector {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );


	typedef t_Type Type;


	inline SparseVector();
	inline SparseVector( SparseVector const& other );

	template< typename t_OtherType >
	inline SparseVector( SparseVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline SparseVector( SpanVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline SparseVector( DenseVector< t_OtherType > const& other );


	inline SparseVector const& operator=( SparseVector const& other );

	template< typename t_OtherType >
	inline SparseVector const& operator=( SparseVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline SparseVector const& operator=( SpanVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline SparseVector const& operator=( DenseVector< t_OtherType > const& other );


	inline unsigned int const Size() const;    // in bytes (approximate)
	inline unsigned int const Dimension() const;


	inline double const NormSquared() const;

	template< typename t_OtherType >
	inline double const InnerProduct( SparseVector< t_OtherType > const& other ) const;

	template< typename t_OtherType >
	inline double const InnerProduct( SpanVector< t_OtherType > const& other ) const;

	template< typename t_OtherType >
	inline double const InnerProduct( DenseVector< t_OtherType > const& other ) const;


	template< typename t_OtherType >
	inline SparseVector const& operator*=( t_OtherType const scale );

	template< typename t_OtherType >
	inline SparseVector const& operator/=( t_OtherType const scale );


	inline void Swap( SparseVector< t_Type > const& other );
	inline void Clear();
	inline void Append( unsigned int const index, t_Type const value );


	template< typename t_OtherType >
	friend std::ostream& operator<<( std::ostream& out, SparseVector< t_OtherType > const& vector );


private:

	template< typename t_OtherType >
	inline void Assign( SparseVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline void Assign( SpanVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline void Assign( DenseVector< t_OtherType > const& other );


	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	std::vector< std::pair< unsigned int, t_Type > > m_data;    // index, value


	template< typename t_OtherType >
	friend struct SparseVector;

	template< typename t_OtherType >
	friend struct SpanVector;

	template< typename t_OtherType >
	friend struct DenseVector;


	friend class boost::serialization::access;
};




//============================================================================
//    SpanVector class
//============================================================================


template< typename t_Type >
struct SpanVector {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );


	typedef t_Type Type;

	enum { DEFAULT_SKIP = ( sizeof( std::pair< unsigned int, unsigned int > ) + sizeof( t_Type ) - 1 ) / sizeof( t_Type ) };
	BOOST_STATIC_ASSERT( DEFAULT_SKIP > 0 );


	inline SpanVector();
	inline SpanVector( SpanVector const& other );

	template< typename t_OtherType >
	inline SpanVector( SparseVector< t_OtherType > const& other, unsigned int const skip = DEFAULT_SKIP );

	template< typename t_OtherType >
	inline SpanVector( SpanVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline SpanVector( DenseVector< t_OtherType > const& other, unsigned int const skip = DEFAULT_SKIP );


	inline SpanVector const& operator=( SpanVector const& other );

	template< typename t_OtherType >
	inline SpanVector const& operator=( SparseVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline SpanVector const& operator=( SpanVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline SpanVector const& operator=( DenseVector< t_OtherType > const& other );


	inline unsigned int const Size() const;    // in bytes (approximate)
	inline unsigned int const Dimension() const;


	inline double const NormSquared() const;

	template< typename t_OtherType >
	inline double const InnerProduct( SparseVector< t_OtherType > const& other ) const;

	template< typename t_OtherType >
	inline double const InnerProduct( SpanVector< t_OtherType > const& other ) const;

	template< typename t_OtherType >
	inline double const InnerProduct( DenseVector< t_OtherType > const& other ) const;


	template< typename t_OtherType >
	inline SpanVector const& operator*=( t_OtherType const scale );

	template< typename t_OtherType >
	inline SpanVector const& operator/=( t_OtherType const scale );


	inline void Swap( SpanVector< t_Type > const& other );
	inline void Clear();
	inline void Append( unsigned int const index, t_Type const value, unsigned int const skip = DEFAULT_SKIP );


	template< typename t_OtherType >
	friend std::ostream& operator<<( std::ostream& out, SpanVector< t_OtherType > const& vector );


private:

	template< typename t_OtherType >
	inline void Assign( SparseVector< t_OtherType > const& other, unsigned int const skip );

	template< typename t_OtherType >
	inline void Assign( SpanVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline void Assign( DenseVector< t_OtherType > const& other, unsigned int const skip );


	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	std::vector< std::pair< unsigned int, unsigned int > > m_spans;    // index, length
	std::vector< t_Type > m_data;


	template< typename t_OtherType >
	friend struct SparseVector;

	template< typename t_OtherType >
	friend struct SpanVector;

	template< typename t_OtherType >
	friend struct DenseVector;


	friend class boost::serialization::access;
};




//============================================================================
//    DenseVector class
//============================================================================


template< typename t_Type >
struct DenseVector {

	BOOST_STATIC_ASSERT( boost::is_floating_point< t_Type >::value );


	typedef t_Type Type;


	inline DenseVector( unsigned int const dimension = 0 );
	inline DenseVector( DenseVector const& other );

	template< typename t_OtherType >
	inline DenseVector( SparseVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline DenseVector( SparseVector< t_OtherType > const& other, unsigned int const dimension );

	template< typename t_OtherType >
	inline DenseVector( SpanVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline DenseVector( SpanVector< t_OtherType > const& other, unsigned int const dimension );

	template< typename t_OtherType >
	inline DenseVector( DenseVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline DenseVector( DenseVector< t_OtherType > const& other, unsigned int const dimension );


	inline DenseVector const& operator=( DenseVector const& other );

	template< typename t_OtherType >
	inline DenseVector const& operator=( SparseVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline DenseVector const& operator=( SpanVector< t_OtherType > const& other );

	template< typename t_OtherType >
	inline DenseVector const& operator=( DenseVector< t_OtherType > const& other );


	inline unsigned int const Size() const;    // in bytes (approximate)
	inline unsigned int const Dimension() const;


	inline double const NormSquared() const;

	template< typename t_OtherType >
	inline double const InnerProduct( SparseVector< t_OtherType > const& other ) const;

	template< typename t_OtherType >
	inline double const InnerProduct( SpanVector< t_OtherType > const& other ) const;

	template< typename t_OtherType >
	inline double const InnerProduct( DenseVector< t_OtherType > const& other ) const;


	template< typename t_OtherType >
	inline DenseVector const& operator*=( t_OtherType const scale );

	template< typename t_OtherType >
	inline DenseVector const& operator/=( t_OtherType const scale );


	inline void Swap( DenseVector< t_Type > const& other );
	inline void Clear();
	inline void Clear( unsigned int const dimension );
	inline void Append( t_Type const value );


	inline t_Type const& operator[]( unsigned int const index ) const;
	inline t_Type& operator[]( unsigned int const index );


	template< typename t_OtherType >
	friend std::ostream& operator<<( std::ostream& out, DenseVector< t_OtherType > const& vector );


private:

	template< typename t_OtherType >
	inline void Assign( SparseVector< t_OtherType > const& other, unsigned int const dimension );

	template< typename t_OtherType >
	inline void Assign( SpanVector< t_OtherType > const& other, unsigned int const dimension );

	template< typename t_OtherType >
	inline void Assign( DenseVector< t_OtherType > const& other, unsigned int const dimension );


	template< typename t_Archive >
	inline void serialize( t_Archive& archive, unsigned int const );


	std::vector< t_Type > m_data;


	template< typename t_OtherType >
	friend struct SparseVector;

	template< typename t_OtherType >
	friend struct SpanVector;

	template< typename t_OtherType >
	friend struct DenseVector;


	friend class boost::serialization::access;
};




//============================================================================
//    SparseVector inline methods
//============================================================================


template< typename t_Type >
SparseVector< t_Type >::SparseVector() {
}


template< typename t_Type >
SparseVector< t_Type >::SparseVector( SparseVector const& other ) : m_data( other.m_data ) {
}


template< typename t_Type >
template< typename t_OtherType >
SparseVector< t_Type >::SparseVector( SparseVector< t_OtherType > const& other ) {

	Assign( other );
}


template< typename t_Type >
template< typename t_OtherType >
SparseVector< t_Type >::SparseVector( SpanVector< t_OtherType > const& other ) {

	Assign( other );
}


template< typename t_Type >
template< typename t_OtherType >
SparseVector< t_Type >::SparseVector( DenseVector< t_OtherType > const& other ) {

	Assign( other );
}


template< typename t_Type >
SparseVector< t_Type > const& SparseVector< t_Type >::operator=( SparseVector const& other ) {

	if ( this != &other )
		m_data = other.m_data;
	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
SparseVector< t_Type > const& SparseVector< t_Type >::operator=( SparseVector< t_OtherType > const& other ) {

	m_data.clear();
	Assign( other );
	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
SparseVector< t_Type > const& SparseVector< t_Type >::operator=( SpanVector< t_OtherType > const& other ) {

	m_data.clear();
	Assign( other );
	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
SparseVector< t_Type > const& SparseVector< t_Type >::operator=( DenseVector< t_OtherType > const& other ) {

	m_data.clear();
	Assign( other );
	return *this;
}


template< typename t_Type >
unsigned int const SparseVector< t_Type >::Size() const {

	return( sizeof( SparseVector ) + m_data.size() * sizeof( std::pair< unsigned int, t_Type > ) );
}


template< typename t_Type >
unsigned int const SparseVector< t_Type >::Dimension() const {

	unsigned int result = 0;
	if ( ! m_data.empty() )
		result = m_data.back().first + 1;
	return result;
}


template< typename t_Type >
double const SparseVector< t_Type >::NormSquared() const {

	double result = 0;

	typename std::vector< std::pair< unsigned int, t_Type > >::const_iterator ii    = m_data.begin();
	typename std::vector< std::pair< unsigned int, t_Type > >::const_iterator iiEnd = m_data.end();
	for ( ; ii != iiEnd; ++ii )
		result += Square( static_cast< double >( ii->second ) );

	return result;
}


template< typename t_Type >
template< typename t_OtherType >
double const SparseVector< t_Type >::InnerProduct( SparseVector< t_OtherType > const& other ) const {

	double result = 0;

	typename std::vector< std::pair< unsigned int, t_Type > >::const_iterator ii    = m_data.begin();
	typename std::vector< std::pair< unsigned int, t_Type > >::const_iterator iiEnd = m_data.end();

	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator jj    = other.m_data.begin();
	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator jjEnd = other.m_data.end();

	while ( ( ii != iiEnd ) && ( jj != jjEnd ) ) {

		if ( ii->first < jj->first )
			++ii;
		else if ( jj->first < ii->first )
			++jj;
		else {

			result += static_cast< double >( ii->second ) * static_cast< double >( jj->second );
			++ii;
			++jj;
		}
	}

	return result;
}


template< typename t_Type >
template< typename t_OtherType >
double const SparseVector< t_Type >::InnerProduct( SpanVector< t_OtherType > const& other ) const {

	double result = 0;

	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator ii    = m_data.begin();
	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator iiEnd = m_data.end();

	std::vector< std::pair< unsigned int, unsigned int > >::const_iterator jj    = other.m_spans.begin();
	std::vector< std::pair< unsigned int, unsigned int > >::const_iterator jjEnd = other.m_spans.end();

	typename std::vector< t_Type >::const_iterator otherData = other.m_data.begin();

	for ( ; ( ii != iiEnd ) && ( jj != jjEnd ); ++jj ) {

		unsigned int const otherIndex = jj->first;
		unsigned int const otherSize  = jj->second;
		unsigned int const otherEnd   = otherIndex + otherSize;

		for ( ; ( ii != iiEnd ) && ( ii->first < otherIndex ); ++ii );
		for ( ; ( ii != iiEnd ) && ( ii->first < otherEnd   ); ++ii )
			result += static_cast< double >( ii->second ) * static_cast< double >( *( otherData + ( ii->first - otherIndex ) ) );

		otherData += otherSize;
	}

	return result;
}


template< typename t_Type >
template< typename t_OtherType >
double const SparseVector< t_Type >::InnerProduct( DenseVector< t_OtherType > const& other ) const {

	double result = 0;

	unsigned int const otherDimension = other.m_data.size();
	unsigned int const dimension = Dimension();
	if ( dimension <= otherDimension ) {

		typename std::vector< std::pair< unsigned int, t_Type > >::const_iterator ii    = m_data.begin();
		typename std::vector< std::pair< unsigned int, t_Type > >::const_iterator iiEnd = m_data.end();

		typename std::vector< t_OtherType >::const_iterator jj = other.m_data.begin();

		for ( ; ii != iiEnd; ++ii )
			result += static_cast< double >( ii->second ) * static_cast< double >( *( jj + ii->first ) );
	}
	else {

		typename std::vector< std::pair< unsigned int, t_Type > >::const_iterator ii = m_data.begin();

		typename std::vector< t_OtherType >::const_iterator jj = other.m_data.begin();

		for ( ; ii->first < otherDimension; ++ii )
			result += static_cast< double >( ii->second ) * static_cast< double >( *( jj + ii->first ) );
	}

	return result;
}


template< typename t_Type >
template< typename t_OtherType >
SparseVector< t_Type > const& SparseVector< t_Type >::operator*=( t_OtherType const scale ) {

	if ( scale == 0 )
		m_data.clear();
	else {

		typename std::vector< std::pair< unsigned int, t_Type > >::iterator ii    = m_data.begin();
		typename std::vector< std::pair< unsigned int, t_Type > >::iterator iiEnd = m_data.end();
		for ( ; ii != iiEnd; ++ii )
			ii->second *= scale;
	}

	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
SparseVector< t_Type > const& SparseVector< t_Type >::operator/=( t_OtherType const scale ) {

	BOOST_ASSERT( scale != 0 );

	typename std::vector< std::pair< unsigned int, t_Type > >::iterator ii    = m_data.begin();
	typename std::vector< std::pair< unsigned int, t_Type > >::iterator iiEnd = m_data.end();
	for ( ; ii != iiEnd; ++ii )
		ii->second /= scale;

	return *this;
}


template< typename t_Type >
void SparseVector< t_Type >::Swap( SparseVector< t_Type > const& other ) {

	m_data.swap( other.m_data );
}


template< typename t_Type >
void SparseVector< t_Type >::Clear() {

	m_data.clear();
}


template< typename t_Type >
void SparseVector< t_Type >::Append( unsigned int const index, t_Type const value ) {

	BOOST_ASSERT( index >= Dimension() );

	if ( value != 0 )
		m_data.push_back( std::pair< unsigned int, t_Type >( index, value ) );
}


template< typename t_Type >
template< typename t_OtherType >
void SparseVector< t_Type >::Assign( SparseVector< t_OtherType > const& other ) {

	BOOST_ASSERT( m_data.empty() );

	m_data.reserve( other.m_data.size() );

	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator ii    = other.m_data.begin();
	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator iiEnd = other.m_data.end();
	for ( ; ii != iiEnd; ++ii )
		m_data.push_back( std::pair< unsigned int, t_Type >( ii->first, ii->second ) );
}


template< typename t_Type >
template< typename t_OtherType >
void SparseVector< t_Type >::Assign( SpanVector< t_OtherType > const& other ) {

	BOOST_ASSERT( m_data.empty() );

	typename std::vector< std::pair< unsigned int, unsigned int > >::const_iterator ii    = other.m_spans.begin();
	typename std::vector< std::pair< unsigned int, unsigned int > >::const_iterator iiEnd = other.m_spans.end();

	typename std::vector< t_OtherType >::const_iterator otherData = other.m_data.begin();

	for ( ; ii != iiEnd; ++ii ) {

		unsigned int const otherIndex = ii->first;
		unsigned int const otherSize  = ii->second;
		unsigned int const otherEnd   = otherIndex + otherSize;

		for ( unsigned int index = otherIndex; index != otherEnd; ++index, ++otherData ) {

			t_Type const value = *otherData;
			if ( value != 0 )
				m_data.push_back( std::pair< unsigned int, t_Type >( index, value ) );
		}
	}
}


template< typename t_Type >
template< typename t_OtherType >
void SparseVector< t_Type >::Assign( DenseVector< t_OtherType > const& other ) {

	BOOST_ASSERT( m_data.empty() );

	typename std::vector< t_OtherType >::const_iterator iiBegin = other.m_data.begin();
	typename std::vector< t_OtherType >::const_iterator ii      = iiBegin;
	typename std::vector< t_OtherType >::const_iterator iiEnd   = other.m_data.end();

	for ( ; ii != iiEnd; ++ii ) {

		t_Type const value = *ii;
		if ( value != 0 )
			m_data.push_back( std::pair< unsigned int, t_Type >( ii - iiBegin, value ) );
	}
}


template< typename t_Type >
template< typename t_Archive >
inline void SparseVector< t_Type >::serialize( t_Archive& archive, unsigned int const ) {

	archive & m_data;
}




//============================================================================
//    SpanVector inline methods
//============================================================================


template< typename t_Type >
SpanVector< t_Type >::SpanVector() {
}


template< typename t_Type >
SpanVector< t_Type >::SpanVector( SpanVector const& other ) : m_spans( other.m_spans ), m_data( other.m_data ) {
}


template< typename t_Type >
template< typename t_OtherType >
SpanVector< t_Type >::SpanVector( SparseVector< t_OtherType > const& other, unsigned int const skip ) {

	Assign( other, skip );
}


template< typename t_Type >
template< typename t_OtherType >
SpanVector< t_Type >::SpanVector( SpanVector< t_OtherType > const& other ) {

	Assign( other );
}


template< typename t_Type >
template< typename t_OtherType >
SpanVector< t_Type >::SpanVector( DenseVector< t_OtherType > const& other, unsigned int const skip ) {

	Assign( other, skip );
}


template< typename t_Type >
SpanVector< t_Type > const& SpanVector< t_Type >::operator=( SpanVector const& other ) {

	if ( this != &other ) {

		m_spans = other.m_spans;
		m_data = other.m_data;
	}
	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
SpanVector< t_Type > const& SpanVector< t_Type >::operator=( SparseVector< t_OtherType > const& other ) {

	m_spans.clear();
	m_data.clear();
	Assign( other, DEFAULT_SKIP );
	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
SpanVector< t_Type > const& SpanVector< t_Type >::operator=( SpanVector< t_OtherType > const& other ) {

	m_data.clear();
	Assign( other );
	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
SpanVector< t_Type > const& SpanVector< t_Type >::operator=( DenseVector< t_OtherType > const& other ) {

	m_spans.clear();
	m_data.clear();
	Assign( other, DEFAULT_SKIP );
	return *this;
}


template< typename t_Type >
unsigned int const SpanVector< t_Type >::Size() const {

	return( sizeof( SpanVector ) + m_spans.size() * sizeof( std::pair< unsigned int, unsigned int > ) + m_data.size() * sizeof( t_Type ) );
}


template< typename t_Type >
unsigned int const SpanVector< t_Type >::Dimension() const {

	unsigned int result = 0;
	if ( ! m_spans.empty() ) {

		std::pair< unsigned int, unsigned int > back = m_spans.back();
		result = back.first + back.second;
	}
	return result;
}


template< typename t_Type >
double const SpanVector< t_Type >::NormSquared() const {

	double result = 0;

	typename std::vector< t_Type >::const_iterator ii    = m_data.begin();
	typename std::vector< t_Type >::const_iterator iiEnd = m_data.end();
	for ( ; ii != iiEnd; ++ii )
		result += Square( static_cast< double >( *ii ) );

	return result;
}


template< typename t_Type >
template< typename t_OtherType >
double const SpanVector< t_Type >::InnerProduct( SparseVector< t_OtherType > const& other ) const {

	return other.InnerProduct( *this );
}


template< typename t_Type >
template< typename t_OtherType >
double const SpanVector< t_Type >::InnerProduct( SpanVector< t_OtherType > const& other ) const {

	double result = 0;

	std::vector< std::pair< unsigned int, unsigned int > >::const_iterator ii    = m_spans.begin();
	std::vector< std::pair< unsigned int, unsigned int > >::const_iterator iiEnd = m_spans.end();

	std::vector< std::pair< unsigned int, unsigned int > >::const_iterator jj    = other.m_spans.begin();
	std::vector< std::pair< unsigned int, unsigned int > >::const_iterator jjEnd = other.m_spans.end();

	typename std::vector< t_Type      >::const_iterator data      = m_data.begin();
	typename std::vector< t_OtherType >::const_iterator otherData = other.m_data.begin();

	while ( ( ii != iiEnd ) && ( jj != jjEnd ) ) {

		unsigned int const index = ii->first;
		unsigned int const size  = ii->second;
		unsigned int const end   = index + size;

		unsigned int const otherIndex = jj->first;
		unsigned int const otherSize  = jj->second;
		unsigned int const otherEnd   = otherIndex + otherSize;

		if ( end <= otherIndex ) {

			data += size;
			++ii;
		}
		else if ( otherEnd <= index ) {

			otherData += otherSize;
			++jj;
		}
		else if ( end <= otherEnd ) {

			unsigned int const overlapIndex = std::max( index, otherIndex );
			BOOST_ASSERT( overlapIndex < end );

			data += overlapIndex - index;

			typename std::vector< t_OtherType >::const_iterator kk    = otherData + ( overlapIndex - otherIndex );
			typename std::vector< t_OtherType >::const_iterator kkEnd = otherData + ( end          - otherIndex );
			for ( ; kk != kkEnd; ++kk, ++data )
				result += static_cast< double >( *data ) * static_cast< double >( *kk );

			++ii;
		}
		else {

			unsigned int const overlapIndex = std::max( index, otherIndex );
			BOOST_ASSERT( overlapIndex < otherEnd );

			otherData += overlapIndex - otherIndex;

			typename std::vector< t_Type >::const_iterator kk    = data + ( overlapIndex - index );
			typename std::vector< t_Type >::const_iterator kkEnd = data + ( otherEnd     - index );
			for ( ; kk != kkEnd; ++kk, ++otherData )
				result += static_cast< double >( *kk ) * static_cast< double >( *otherData );

			++jj;
		}
	}

	return result;
}


template< typename t_Type >
template< typename t_OtherType >
double const SpanVector< t_Type >::InnerProduct( DenseVector< t_OtherType > const& other ) const {

	double result = 0;

	unsigned int const otherDimension = other.m_data.size();
	unsigned int const dimension = Dimension();
	if ( dimension <= otherDimension ) {

		std::vector< std::pair< unsigned int, unsigned int > >::const_iterator ii    = m_spans.begin();
		std::vector< std::pair< unsigned int, unsigned int > >::const_iterator iiEnd = m_spans.end();

		typename std::vector< t_OtherType >::const_iterator jjBegin = other.m_data.begin();

		typename std::vector< t_Type >::const_iterator data = m_data.begin();

		while ( ii != iiEnd ) {

			unsigned int const index = ii->first;
			unsigned int const size  = ii->second;

			typename std::vector< t_OtherType >::const_iterator jj    = jjBegin + index;
			typename std::vector< t_OtherType >::const_iterator jjEnd = jj + size;
			for ( ; jj != jjEnd; ++jj, ++data )
				result += static_cast< double >( *data ) * static_cast< double >( *jj );

			++ii;
		}
	}
	else {

		std::vector< std::pair< unsigned int, unsigned int > >::const_iterator ii = m_spans.begin();

		typename std::vector< t_OtherType >::const_iterator jjBegin = other.m_data.begin();

		typename std::vector< t_Type >::const_iterator data = m_data.begin();

		for ( ; ; ) {

			unsigned int const index = ii->first;
			unsigned int const size  = ii->second;
			unsigned int const end   = index + size;

			if ( end >= otherDimension ) {

				typename std::vector< t_OtherType >::const_iterator jj    = jjBegin + index;
				typename std::vector< t_OtherType >::const_iterator jjEnd = jjBegin + otherDimension;
				for ( ; jj < jjEnd; ++jj, ++data )
					result += static_cast< double >( *data ) * static_cast< double >( *jj );

				break;
			}
			else {

				typename std::vector< t_OtherType >::const_iterator jj    = jjBegin + index;
				typename std::vector< t_OtherType >::const_iterator jjEnd = jj + size;
				for ( ; jj != jjEnd; ++jj, ++data )
					result += static_cast< double >( *data ) * static_cast< double >( *jj );
			}

			++ii;
		}
	}

	return result;
}


template< typename t_Type >
template< typename t_OtherType >
SpanVector< t_Type > const& SpanVector< t_Type >::operator*=( t_OtherType const scale ) {

	if ( scale == 0 ) {

		m_spans.clear();
		m_data.clear();
	}
	else {

		typename std::vector< t_Type >::iterator ii    = m_data.begin();
		typename std::vector< t_Type >::iterator iiEnd = m_data.end();
		for ( ; ii != iiEnd; ++ii )
			*ii *= scale;
	}

	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
SpanVector< t_Type > const& SpanVector< t_Type >::operator/=( t_OtherType const scale ) {

	BOOST_ASSERT( scale != 0 );

	typename std::vector< t_Type >::iterator ii    = m_data.begin();
	typename std::vector< t_Type >::iterator iiEnd = m_data.end();
	for ( ; ii != iiEnd; ++ii )
		*ii /= scale;

	return *this;
}


template< typename t_Type >
void SpanVector< t_Type >::Swap( SpanVector< t_Type > const& other ) {

	m_spans.swap( other.m_spans );
	m_data.swap(  other.m_data  );
}


template< typename t_Type >
void SpanVector< t_Type >::Clear() {

	m_spans.clear();
	m_data.clear();
}


template< typename t_Type >
void SpanVector< t_Type >::Append( unsigned int const index, t_Type const value, unsigned int const skip ) {

	BOOST_ASSERT( index >= Dimension() );

	if ( value != 0 ) {

		if ( m_spans.empty() ) {

			m_spans.push_back( std::pair< unsigned int, unsigned int >( index, 1 ) );
			m_data.push_back( value );
		}
		else {

			std::pair< unsigned int, unsigned int >& back = m_spans.back();
			if ( index <= back.first + back.second + skip ) {

				for ( ; back.first + back.second < index; ++back.second )
					m_data.push_back( 0 );
				m_data.push_back( value );
				++back.second;
			}
			else {

				m_spans.push_back( std::pair< unsigned int, unsigned int >( index, 1 ) );
				m_data.push_back( value );
			}
		}
	}
}


template< typename t_Type >
template< typename t_OtherType >
void SpanVector< t_Type >::Assign( SparseVector< t_OtherType > const& other, unsigned int const skip ) {

	BOOST_ASSERT( m_spans.empty() );
	BOOST_ASSERT( m_data.empty() );

	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator ii    = other.m_data.begin();
	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator iiEnd = other.m_data.end();

	if ( ii != iiEnd ) {

		m_data.push_back( ii->second );
		unsigned int spanIndex = ii->first;
		unsigned int spanSize = 1;
		++ii;

		for ( ; ii != iiEnd; ++ii ) {

			unsigned int const index = ii->first;
			if ( index <= spanIndex + spanSize + skip ) {

				for ( ; spanIndex + spanSize < index; ++spanSize )
					m_data.push_back( 0 );
				m_data.push_back( ii->second );
				++spanSize;
			}
			else {

				m_spans.push_back( std::pair< unsigned int, unsigned int >( spanIndex, spanSize ) );

				m_data.push_back( ii->second );
				spanIndex = index;
				spanSize = 1;
			}
		}

		m_spans.push_back( std::pair< unsigned int, unsigned int >( spanIndex, spanSize ) );
	}
}


template< typename t_Type >
template< typename t_OtherType >
void SpanVector< t_Type >::Assign( SpanVector< t_OtherType > const& other ) {

	BOOST_ASSERT( m_data.empty() );

	m_spans = other.m_spans;
	m_data.reserve( other.m_data.size() );

	typename std::vector< t_OtherType >::const_iterator ii    = other.m_data.begin();
	typename std::vector< t_OtherType >::const_iterator iiEnd = other.m_data.end();
	for ( ; ii != iiEnd; ++ii )
		m_data.push_back( *ii );
}


template< typename t_Type >
template< typename t_OtherType >
void SpanVector< t_Type >::Assign( DenseVector< t_OtherType > const& other, unsigned int const skip ) {

	BOOST_ASSERT( m_spans.empty() );
	BOOST_ASSERT( m_data.empty() );

	typename std::vector< t_OtherType >::const_iterator ii    = other.m_data.begin();
	typename std::vector< t_OtherType >::const_iterator iiEnd = other.m_data.end();

	unsigned int index = 0;

	for ( ; ( ii != iiEnd ) && ( static_cast< t_Type >( *ii ) == 0 ); ++ii, ++index );

	if ( ii != iiEnd ) {

		m_data.push_back( *ii );
		unsigned int spanIndex = index;
		unsigned int spanSize = 1;

		for ( ++ii, ++index; ( ii != iiEnd ) && ( static_cast< t_Type >( *ii ) == 0 ); ++ii, ++index );

		while ( ii != iiEnd ) {

			if ( index <= spanIndex + spanSize + skip ) {

				for ( ; spanIndex + spanSize < index; ++spanSize )
					m_data.push_back( 0 );
				m_data.push_back( *ii );
				++spanSize;
			}
			else {

				m_spans.push_back( std::pair< unsigned int, unsigned int >( spanIndex, spanSize ) );

				m_data.push_back( *ii );
				spanIndex = index;
				spanSize = 1;
			}

			for ( ++ii, ++index; ( ii != iiEnd ) && ( static_cast< t_Type >( *ii ) == 0 ); ++ii, ++index );
		}

		m_spans.push_back( std::pair< unsigned int, unsigned int >( spanIndex, spanSize ) );
	}
}


template< typename t_Type >
template< typename t_Archive >
inline void SpanVector< t_Type >::serialize( t_Archive& archive, unsigned int const ) {

	archive & m_spans;
	archive & m_data;
}




//============================================================================
//    DenseVector inline methods
//============================================================================


template< typename t_Type >
DenseVector< t_Type >::DenseVector( unsigned int const dimension ) {

	m_data.reserve( dimension );
	for ( unsigned int ii = 0; ii < dimension; ++ii )
		m_data.push_back( 0 );
}


template< typename t_Type >
DenseVector< t_Type >::DenseVector( DenseVector const& other ) : m_data( other.m_data ) {
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type >::DenseVector( SparseVector< t_OtherType > const& other ) {

	Assign( other, other.Dimension() );
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type >::DenseVector( SparseVector< t_OtherType > const& other, unsigned int const dimension ) {

	Assign( other, dimension );
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type >::DenseVector( SpanVector< t_OtherType > const& other ) {

	Assign( other, other.Dimension() );
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type >::DenseVector( SpanVector< t_OtherType > const& other, unsigned int const dimension ) {

	Assign( other, dimension );
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type >::DenseVector( DenseVector< t_OtherType > const& other ) {

	Assign( other, other.m_data.size() );
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type >::DenseVector( DenseVector< t_OtherType > const& other, unsigned int const dimension ) {

	Assign( other, dimension );
}


template< typename t_Type >
DenseVector< t_Type > const& DenseVector< t_Type >::operator=( DenseVector const& other ) {

	if ( this != &other )
		m_data = other.m_data;
	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type > const& DenseVector< t_Type >::operator=( SparseVector< t_OtherType > const& other ) {

	m_data.clear();
	Assign( other, other.Dimension() );
	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type > const& DenseVector< t_Type >::operator=( SpanVector< t_OtherType > const& other ) {

	m_data.clear();
	Assign( other, other.Dimension() );
	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type > const& DenseVector< t_Type >::operator=( DenseVector< t_OtherType > const& other ) {

	m_data.clear();
	Assign( other, other.m_data.size() );
	return *this;
}


template< typename t_Type >
unsigned int const DenseVector< t_Type >::Size() const {

	return( sizeof( DenseVector ) + m_data.size() * sizeof( t_Type ) );
}


template< typename t_Type >
unsigned int const DenseVector< t_Type >::Dimension() const {

	return m_data.size();
}


template< typename t_Type >
double const DenseVector< t_Type >::NormSquared() const {

	double result = 0;

	typename std::vector< t_Type >::const_iterator ii    = m_data.begin();
	typename std::vector< t_Type >::const_iterator iiEnd = m_data.end();
	for ( ; ii != iiEnd; ++ii )
		result += Square( static_cast< double >( *ii ) );

	return result;
}


template< typename t_Type >
template< typename t_OtherType >
double const DenseVector< t_Type >::InnerProduct( SparseVector< t_OtherType > const& other ) const {

	return other.InnerProduct( *this );
}


template< typename t_Type >
template< typename t_OtherType >
double const DenseVector< t_Type >::InnerProduct( SpanVector< t_OtherType > const& other ) const {

	return other.InnerProduct( *this );
}


template< typename t_Type >
template< typename t_OtherType >
double const DenseVector< t_Type >::InnerProduct( DenseVector< t_OtherType > const& other ) const {

	double result = 0;

	unsigned int const dimension = std::min( m_data.size(), other.m_data.size() );

	typename std::vector< t_Type >::const_iterator ii    = m_data.begin();
	typename std::vector< t_Type >::const_iterator iiEnd = ii + dimension;

	typename std::vector< t_OtherType >::const_iterator jj = other.m_data.begin();

	for ( ; ii != iiEnd; ++ii, ++jj )
		result += static_cast< double >( *ii ) * static_cast< double >( *jj );

	return result;
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type > const& DenseVector< t_Type >::operator*=( t_OtherType const scale ) {

	typename std::vector< t_Type >::iterator ii    = m_data.begin();
	typename std::vector< t_Type >::iterator iiEnd = m_data.end();
	for ( ; ii != iiEnd; ++ii )
		*ii *= scale;

	return *this;
}


template< typename t_Type >
template< typename t_OtherType >
DenseVector< t_Type > const& DenseVector< t_Type >::operator/=( t_OtherType const scale ) {

	BOOST_ASSERT( scale != 0 );

	typename std::vector< t_Type >::iterator ii    = m_data.begin();
	typename std::vector< t_Type >::iterator iiEnd = m_data.end();
	for ( ; ii != iiEnd; ++ii )
		*ii /= scale;

	return *this;
}


template< typename t_Type >
void DenseVector< t_Type >::Swap( DenseVector< t_Type > const& other ) {

	m_data.swap( other.m_data );
}


template< typename t_Type >
void DenseVector< t_Type >::Clear() {

	m_data.clear();
}


template< typename t_Type >
void DenseVector< t_Type >::Clear( unsigned int const dimension ) {

	if ( ! m_data.empty() )
		std::fill( m_data.begin(), m_data.end(), 0 );
}


template< typename t_Type >
void DenseVector< t_Type >::Append( t_Type const value ) {

	m_data.push_back( value );
}


template< typename t_Type >
t_Type const& DenseVector< t_Type >::operator[]( unsigned int const index ) const {

	BOOST_ASSERT( index < m_data.size() );
	return m_data[ index ];
}


template< typename t_Type >
t_Type& DenseVector< t_Type >::operator[]( unsigned int const index ) {

	BOOST_ASSERT( index < m_data.size() );
	return m_data[ index ];
}


template< typename t_Type >
template< typename t_OtherType >
void DenseVector< t_Type >::Assign( SparseVector< t_OtherType > const& other, unsigned int const dimension ) {

	BOOST_ASSERT( m_data.empty() );

	BOOST_ASSERT( dimension >= other.Dimension() );
	m_data.reserve( dimension );

	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator ii    = other.m_data.begin();
	typename std::vector< std::pair< unsigned int, t_OtherType > >::const_iterator iiEnd = other.m_data.end();

	unsigned int index = 0;
	for ( ; ii != iiEnd; ++ii ) {

		for ( ; index < ii->first; ++index )
			m_data.push_back( 0 );
		m_data.push_back( ii->second );
		++index;
	}
	for ( ; index < dimension; ++index )
		m_data.push_back( 0 );
	BOOST_ASSERT( index == dimension );
}


template< typename t_Type >
template< typename t_OtherType >
void DenseVector< t_Type >::Assign( SpanVector< t_OtherType > const& other, unsigned int const dimension ) {

	BOOST_ASSERT( m_data.empty() );

	BOOST_ASSERT( dimension >= other.Dimension() );
	m_data.reserve( dimension );

	typename std::vector< std::pair< unsigned int, unsigned int > >::const_iterator ii    = other.m_spans.begin();
	typename std::vector< std::pair< unsigned int, unsigned int > >::const_iterator iiEnd = other.m_spans.end();

	typename std::vector< t_OtherType >::const_iterator otherData = other.m_data.begin();

	unsigned int index = 0;
	for ( ; ii != iiEnd; ++ii ) {

		unsigned int const otherIndex = ii->first;
		unsigned int const otherSize  = ii->second;
		unsigned int const otherEnd   = otherIndex + otherSize;

		for ( ; index < otherIndex; ++index )
			m_data.push_back( 0 );
		for ( ; index < otherEnd; ++index, ++otherData )
			m_data.push_back( *otherData );
	}
	for ( ; index < dimension; ++index )
		m_data.push_back( 0 );
	BOOST_ASSERT( index == dimension );
}


template< typename t_Type >
template< typename t_OtherType >
void DenseVector< t_Type >::Assign( DenseVector< t_OtherType > const& other, unsigned int const dimension ) {

	BOOST_ASSERT( m_data.empty() );

	BOOST_ASSERT( dimension >= other.m_data.size() );
	m_data.reserve( dimension );

	typename std::vector< t_OtherType >::const_iterator ii    = other.m_data.begin();
	typename std::vector< t_OtherType >::const_iterator iiEnd = other.m_data.end();
	for ( ; ii != iiEnd; ++ii )
		m_data.push_back( *ii );
	for ( unsigned int index = other.m_data.size(); index < dimension; ++index )
		m_data.push_back( 0 );
}


template< typename t_Type >
template< typename t_Archive >
inline void DenseVector< t_Type >::serialize( t_Archive& archive, unsigned int const ) {

	archive & m_data;
}




//============================================================================
//    printing operators
//============================================================================


template< typename t_Type >
std::ostream& operator<<( std::ostream& out, SparseVector< t_Type > const& vector ) {

	if ( vector.m_data.size() == 0 )
		out << "[]";
	else {

		typename std::vector< std::pair< unsigned int, t_Type > >::const_iterator ii    = vector.m_data.begin();
		typename std::vector< std::pair< unsigned int, t_Type > >::const_iterator iiEnd = vector.m_data.end();

		out << '[' << ii->first << ':' << ii->second;
		for ( ++ii; ii != iiEnd; ++ii )
			out << ',' << ii->first << ':' << ii->second;
		out << ']';
	}

	return out;
}


template< typename t_Type >
std::ostream& operator<<( std::ostream& out, SpanVector< t_Type > const& vector ) {

	if ( vector.m_spans.size() == 0 )
		out << "[]";
	else {

		typename std::vector< std::pair< unsigned int, unsigned int > >::const_iterator ii    = vector.m_spans.begin();
		typename std::vector< std::pair< unsigned int, unsigned int > >::const_iterator iiEnd = vector.m_spans.end();

		typename std::vector< t_Type >::const_iterator data = vector.m_data.begin();

		out << '[';
		{	unsigned int const index = ii->first;
			unsigned int const size  = ii->second;
			unsigned int const end   = index + size;

			if ( size == 1 )
				out << index << ':' << *( data++ );
			else {

				out << index << '-' << ( end - 1 ) << ":[" << *( data++ );
				for ( unsigned int jj = index + 1; jj < end; ++jj )
					out << ',' << *( data++ );
				out << ']';
			}
		}
		for ( ++ii; ii != iiEnd; ++ii ) {

			unsigned int const index = ii->first;
			unsigned int const size  = ii->second;
			unsigned int const end   = index + size;

			if ( size == 1 )
				out << ',' << index << ':' << *( data++ );
			else {

				out << ',' << index << '-' << ( end - 1 ) << ":[" << *( data++ );
				for ( unsigned int jj = index + 1; jj < end; ++jj )
					out << ',' << *( data++ );
				out << ']';
			}
		}
		out << ']';
	}

	return out;
}


template< typename t_Type >
std::ostream& operator<<( std::ostream& out, DenseVector< t_Type > const& vector ) {

	if ( vector.m_data.size() == 0 )
		out << "[]";
	else {

		typename std::vector< t_Type >::const_iterator ii    = vector.m_data.begin();
		typename std::vector< t_Type >::const_iterator iiEnd = vector.m_data.end();

		out << '[' << *ii;
		for ( ++ii; ii != iiEnd; ++ii )
			out << ',' << *ii;
		out << ']';
	}

	return out;
}




#endif    /* __cplusplus */

#endif    /* __VECTOR_HPP__ */
