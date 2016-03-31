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
	\file issvm_evaluate.cpp
*/




// must include archive first!
#include <boost/archive/binary_iarchive.hpp>

#include "svm.hpp"
#include "svm_serialization.hpp"
#include "vector.hpp"
#include "data.hpp"
#include "helpers.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>


double sign(double x){
	if (x >= 0)
		return 1;
	else
		return -1;
}

double miss_loss(boost::shared_ptr< SVM::Optimizer::Base > pOptimizer, std::string testingFilename){
	std::vector< SparseVector< float > > vectors;
	std::vector< double > labels;
	LoadDataset( testingFilename.c_str(), vectors, labels );
	unsigned int const size = vectors.size();
	BOOST_ASSERT( size == labels.size() );

	Random::Generator::LinearCongruential<> seedGenerator;
	seedGenerator.Seed();

	boost::shared_array< double > classifications( new double[ size ] );
	#pragma omp parallel
	{	Random::Generator::LaggedFibonacci4<> generator;
	#pragma omp critical
		generator.Seed( seedGenerator );

	#pragma omp for schedule( static )
		for ( int ii = 0; ii < static_cast< int >( size ); ++ii )
			classifications[ ii ] = pOptimizer->Evaluate( generator, vectors[ ii ] );
	}

	//double const* tLabel    = labels;
	double mistakes = 0;
	{	//std::ofstream outputFile( output.c_str() );
		//if ( outputFile.fail() )
		//	throw std::runtime_error( "Unable to open output file" );
		for ( unsigned int ii = 0; ii < size; ++ii) {
			if(sign(classifications[ii])!=labels[ii])
				mistakes+=1;
			//outputFile << classifications[ii] << std::endl;
		}
		double miss_loss = float(mistakes)/float(size);
		return miss_loss;

		//outputFile.close();
	}
}



//============================================================================
//    main function
//============================================================================


int main( int argc, char* argv[] ) {

	int resultCode = EXIT_SUCCESS;

	std::string testingFilename;
	std::string input;
	std::string output;

	boost::program_options::options_description description( "Allowed options" );
	description.add_options()
		( "help,h", "display this help" )
		( "file,f", boost::program_options::value< std::string >( &testingFilename ), "testing dataset file (SVM-Light format)" )
		( "input,i", boost::program_options::value< std::string >( &input ), "input model file" )
	;

	try {

		boost::program_options::variables_map variables;
		boost::program_options::store( boost::program_options::command_line_parser( argc, argv ).options( description ).run(), variables );
		boost::program_options::notify( variables );

		if ( variables.count( "help" ) ) {

			std::cout <<
				"Classifies a dataset (in SVM-Light format). The result is saved to a text file" << std::endl <<
				"containing one floating-point number per line: the value of the classification" << std::endl <<
				"function applied the corresponding testing vector (the signs of these numbers" << std::endl <<
				"are the classifications)." << std::endl <<
				std::endl <<
				description << std::endl;
		}
		else {

			if (!variables.count("file"))
				throw std::runtime_error("You must provide a dataset file");
			if (!variables.count("input"))
				throw std::runtime_error("You must provide an input file");


			boost::shared_ptr <SVM::Optimizer::Base> pOptimizer;
			{
				std::ifstream inputFile(input.c_str(), std::ios::in | std::ios::binary);
				if (inputFile.fail())
					throw std::runtime_error("Unable to open input file");
				boost::iostreams::filtering_streambuf <boost::iostreams::input> inputStream;
				inputStream.push(inputFile);
				boost::archive::binary_iarchive inputArchive(inputStream);
				inputArchive >> pOptimizer;
			}

			/*if(0)
				std::cout << miss_loss(pOptimizer, testingFilename) << std::endl;
			else {
				double sum_miss_loss = 0;
				for (int i = 1; i <= 10; i++) {
					std::string testingFilenameIt = testingFilename + std::to_string(i);
					//std::cout << testingFilenameIt << std::endl;
					double missIt= miss_loss(pOptimizer, testingFilenameIt);
					sum_miss_loss += missIt;
					//std::cout << missIt << std::endl;
				}
			}*/

			double miss_loss_error = miss_loss(pOptimizer, testingFilename);
			std::cout << input << std::endl;
			std::cout << pOptimizer->Support() << " ";
			std::cout << miss_loss_error << " ";
			std::cout << 0 << " ";
			std::cout << pOptimizer->L1Norm() << std::endl;

		}
	}
	catch( std::exception& error ) {

		std::cerr << "Error: " << error.what() << std::endl << std::endl << description << std::endl;
		resultCode = EXIT_FAILURE;
	}

	return resultCode;
}
