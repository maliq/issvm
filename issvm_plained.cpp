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
	\file issvm_optimize.cpp
*/




// must include archives first!
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "svm.hpp"
#include "svm_serialization.hpp"
#include "random.hpp"
#include "helpers.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

#include <fstream>
#include <iostream>
#include <stdexcept>




//============================================================================
//    main function
//============================================================================


int main( int argc, char* argv[] ) {

	int resultCode = EXIT_SUCCESS;

	std::string input;
	std::string output;
	unsigned int iterations;

	boost::program_options::options_description description( "Allowed options" );
	description.add_options()
		( "help,h", "display this help" )
		( "input,i", boost::program_options::value< std::string >( &input ), "input model file" )
	;

	try {

		boost::program_options::variables_map variables;
		boost::program_options::store( boost::program_options::command_line_parser( argc, argv ).options( description ).run(), variables );
		boost::program_options::notify( variables );

		if ( variables.count( "help" ) ) {

			std::cout <<
				"Write the SupportSet on text plain file with extension -SupportSet.txt" << std::endl <<
				"to the output model file. Optimization will continue until the maximum number" << std::endl <<
				"of iterations has been exceeded." << std::endl <<
				std::endl <<
				description << std::endl;
		}
		else {

			if ( ! variables.count( "input" ) )
				throw std::runtime_error( "You must provide an input file" );

			boost::shared_ptr< SVM::Optimizer::Base > pOptimizer;
			{	std::ifstream inputFile( input.c_str(), std::ios::in | std::ios::binary );
				if ( inputFile.fail() )
					throw std::runtime_error( "Unable to open input file" );
				boost::iostreams::filtering_streambuf< boost::iostreams::input > inputStream;
				inputStream.push( boost::iostreams::gzip_decompressor() );
				inputStream.push( inputFile );
				boost::archive::binary_iarchive inputArchive( inputStream );
				inputArchive >> pOptimizer;
			}
			pOptimizer->WriteSupport(input+"-SupportSet.txt");
		}
	}
	catch( std::exception& error ) {

		std::cerr << "Error: " << error.what() << std::endl << std::endl << description << std::endl;
		resultCode = EXIT_FAILURE;
	}

	return resultCode;
}
