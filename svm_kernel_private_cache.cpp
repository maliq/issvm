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
	\file svm_kernel_private_cache.cpp
	\brief SVM::Kernel::_Private::Cache implementation
*/




#include "svm_kernel_private_cache.hpp"




namespace SVM {


namespace Kernel {


namespace _Private {




//============================================================================
//    Kernel::_Private::Cache methods
//============================================================================


std::pair< bool, boost::shared_array< double > > const Cache::operator[]( unsigned int const index ) {

	BOOST_ASSERT( m_dimension > 0 );
	BOOST_ASSERT( index < m_dimension );

	boost::shared_array< double > result;
	bool cached = false;

	if ( m_size == 0 )
		result = boost::shared_array< double >( new double[ m_dimension ] );
	else {

		if ( UNLIKELY( m_cache.size() < m_size ) )
			result = boost::shared_array< double >( new double[ m_dimension ] );
		else {

			std::list< Element >::iterator ii = m_table[ index ];
			if ( LIKELY( ii == m_cache.end() ) ) {

				Element const element = m_cache.back();

				m_table[ element.first ] = m_cache.end();
				m_cache.pop_back();

				if ( element.second.unique() )
					result = element.second;
				else
					result = boost::shared_array< double >( new double[ m_dimension ] );
			}
			else {

				BOOST_ASSERT( ii->first == index );
				result = ii->second;

				m_table[ ii->first ] = m_cache.end();
				m_cache.erase( ii );

				cached = true;
			}
		}
		BOOST_ASSERT( m_cache.size() < m_size );

		m_table[ index ] = m_cache.insert( m_cache.begin(), Element( index, result ) );
	}

	return std::pair< bool, boost::shared_array< double > >( cached, result );
}




}    // namespace _Private


}    // namespace Kernel


}    // namespace SVM
