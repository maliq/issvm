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
	\file data.hpp
	\brief Implementation of LoadDataset function
*/




#ifndef __DATA_HPP__
#define __DATA_HPP__

#ifdef __cplusplus




#include "vector.hpp"

#include <boost/algorithm/string.hpp>
#include <boost/scoped_array.hpp>
#include <boost/regex.hpp>

#include <sstream>
#include <fstream>
#include <vector>
#include <string>




//============================================================================
//    LoadDataset function
//============================================================================


template< typename t_Type >
void LoadDataset(
	std::string const& filename,
	std::vector< SparseVector< t_Type > >& vectors,
	std::vector< double >& labels
)
{
	BOOST_ASSERT( vectors.size() == labels.size() );

	unsigned int const bufferSize = ( 1u << 20 );    // one megabyte per line should suffice
	boost::scoped_array< char > buffer( new char[ bufferSize ] );

	boost::regex const spaceRegex( "[[:space:]]*[,[:space:]][[:space:]]*" );
	boost::regex const elementRegex( "^([^:]*):([^:]*)$" );

	std::ifstream file( filename.c_str() );
	if ( file.fail() ) {

		std::stringstream stream;
		stream << "Unable to open dataset \"" << filename << "\" for reading: " << strerror( errno );
		throw std::runtime_error( stream.str() );
	}
	while ( ! file.eof() ) {

		file.getline( buffer.get(), bufferSize );
		if ( file.fail() )
			break;
		std::string lineString( buffer.get() );

		// remove comments
		unsigned int const commentIndex = lineString.find_first_of( "%#" );
		if ( commentIndex < lineString.size() )
			lineString.erase( commentIndex );

		boost::sregex_token_iterator ii(
			lineString.begin(),
			lineString.end(),
			spaceRegex,
			-1
		);
		boost::sregex_token_iterator iiEnd;

		if ( ii != iiEnd ) {    // ignore blank lines

			std::string const labelString = boost::algorithm::trim_copy( ii->str() );
			double const label = atof( labelString.c_str() );

			// if no feature index is provided, use the previous index plus 1 (starting at 0)
			unsigned int index = 0;
			SparseVector< t_Type > vector;
			for ( ++ii; ii != iiEnd; ++ii ) {

				double value = 0;
				std::string const elementString = boost::algorithm::trim_copy( ii->str() );
				boost::smatch elementMatch;
				if ( boost::regex_match( elementString, elementMatch, elementRegex, boost::match_extra ) ) {

					std::string const indexString = elementMatch[ 1 ].str();
					int const newIndex = atoi( indexString.c_str() );
					if ( newIndex < static_cast< int >( index ) ) {

						std::stringstream stream;
						stream << "Decreasing feature index on line " << ( vectors.size() + 1 ) << " of dataset \"" << filename << "\"";
						throw std::runtime_error( stream.str() );
					}
					index = newIndex;

					std::string const valueString = elementMatch[ 2 ].str();
					value = atof( valueString.c_str() );
				}
				else
					value = atof( elementString.c_str() );

				if ( value != 0 ){
					vector.Append( index, value );
					if(index < 1.0){
						std::cout << "Index:" << index << std::endl;
						std::cout << "Value:" << value << std::endl;
					}

				}
				++index;
			}

			labels.push_back(  label );
			vectors.push_back( vector );

			BOOST_ASSERT( labels.size() == vectors.size() );
		}
	}
	file.close();

	BOOST_ASSERT( vectors.size() == labels.size() );
}




//============================================================================
//    LoadVector function
//============================================================================


template< typename t_Type >
void LoadVector(
	std::string const& filename,
	std::vector< t_Type >& data
)
{
	unsigned int const bufferSize = ( 1u << 20 );    // one megabyte per line should suffice
	boost::scoped_array< char > buffer( new char[ bufferSize ] );

	std::ifstream file( filename.c_str() );
	if ( file.fail() ) {

		std::stringstream stream;
		stream << "Unable to open dataset \"" << filename << "\" for reading: " << strerror( errno );
		throw std::runtime_error( stream.str() );
	}
	while ( ! file.eof() ) {

		file.getline( buffer.get(), bufferSize );
		if ( file.fail() )
			break;
		std::string lineString( buffer.get() );

		// remove comments
		unsigned int const commentIndex = lineString.find_first_of( "%#" );
		if ( commentIndex < lineString.size() )
			lineString.erase( commentIndex );

		std::string const elementString = boost::algorithm::trim_copy( lineString );
		if ( ! elementString.empty() )    // ignore blank lines
			data.push_back( atof( elementString.c_str() ) );
	}
	file.close();
}




#endif    /* __cplusplus */

#endif    /* __DATA_HPP__ */
