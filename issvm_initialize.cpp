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
	\file issvm_initialize.cpp
*/




// must include archive first!
#include <boost/archive/binary_oarchive.hpp>

#include "svm.hpp"
#include "svm_serialization.hpp"
#include "vector.hpp"
#include "data.hpp"
#include "helpers.h"

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/shared_ptr.hpp>

#include <fstream>
#include <iostream>
#include <vector>
#include <stdexcept>
#include <limits>




//============================================================================
//    ConstructKernel and ConstructOptimizer helper functions
//============================================================================


template< typename t_Type >
boost::shared_ptr< SVM::Kernel::Base > ConstructKernel(
	std::string const& trainingFilename,
	std::string const& validationFilename,
	std::string const& name,
	std::vector< std::string > const& parameters,
	unsigned int const cacheSize = 0
)
{
	boost::shared_ptr< SVM::Kernel::Base > result;

	std::vector< SparseVector< t_Type > > vectors;
	std::vector< double > labels;
	LoadDataset< t_Type >( trainingFilename, vectors, labels );
	unsigned int const trainingSize = vectors.size();
	if ( ! validationFilename.empty() )
		LoadDataset< t_Type >( validationFilename, vectors, labels );

	if ( boost::iequals( name, "linear" ) )
		result = SVM::Kernel::Construct< t_Type, SVM::Kernel::Simple, SVM::Kernel::Traits::Linear >( vectors, labels, trainingSize, parameters, cacheSize );
	else if ( boost::iequals( name, "gaussian" ) )
		result = SVM::Kernel::Construct< t_Type, SVM::Kernel::VectorData, SVM::Kernel::Traits::Gaussian >( vectors, labels, trainingSize, parameters, cacheSize );
	else
		throw std::runtime_error( "kernel must be one of \"linear\" or \"gaussian\"" );

	return result;
}


boost::shared_ptr< SVM::Optimizer::Base > ConstructOptimizer(
	boost::shared_ptr< SVM::Kernel::Base > const& pKernel,
	std::string const& name,
	std::vector< std::string > const& parameters,
	bool const biased
)
{
	boost::shared_ptr< SVM::Optimizer::Base > result;

	if ( boost::iequals( name, "sparsifier" ) ) {

		if ( biased )
			result = boost::shared_ptr< SVM::Optimizer::Base >( new SVM::Optimizer::Classification::Biased::Sparsifier( pKernel, parameters ) );
		else
			result = boost::shared_ptr< SVM::Optimizer::Base >( new SVM::Optimizer::Classification::Unbiased::Sparsifier( pKernel, parameters ) );
	}
	else if ( boost::iequals( name, "perceptron" ) ) {

		if ( biased )
			result = boost::shared_ptr< SVM::Optimizer::Base >( new SVM::Optimizer::Classification::Biased::Perceptron( pKernel, parameters ) );
		else
			result = boost::shared_ptr< SVM::Optimizer::Base >( new SVM::Optimizer::Classification::Unbiased::Perceptron( pKernel, parameters ) );
	}
	else if ( boost::iequals( name, "SBP" ) ) {

		if ( biased )
			result = boost::shared_ptr< SVM::Optimizer::Base >( new SVM::Optimizer::Classification::Biased::SBP( pKernel, parameters ) );
		else
			result = boost::shared_ptr< SVM::Optimizer::Base >( new SVM::Optimizer::Classification::Unbiased::SBP( pKernel, parameters ) );
	}
	else if ( boost::iequals( name, "SMO" ) ) {

		if ( biased )
			result = boost::shared_ptr< SVM::Optimizer::Base >( new SVM::Optimizer::Classification::Biased::SMO( pKernel, parameters ) );
		else {
			result = boost::shared_ptr< SVM::Optimizer::Base >( new SVM::Optimizer::Classification::Unbiased::SMO( pKernel, parameters ) );
		}
	}
	else
		throw std::runtime_error( "classification algorithm must be one of \"sparsifier\", \"perceptron\", \"SBP\" or \"SMO\"" );

	return result;
}




//============================================================================
//    main function
//============================================================================


int main( int argc, char* argv[] ) {

	int resultCode = EXIT_SUCCESS;

	std::string trainingFilename;
	std::string validationFilename;
	std::string output;
	std::string kernelName;
	std::vector< std::string > kernelParameters;
	unsigned int cacheSize;
	std::string algorithmName;
	std::vector< std::string > algorithmParameters;
	bool biased;

	boost::program_options::options_description description( "Allowed options" );
	description.add_options()
		( "help,h", "display this help" )
		( "file,f", boost::program_options::value< std::string >( &trainingFilename ), "training dataset file (SVM-Light format)" )
		( "validation_file,v", boost::program_options::value< std::string >( &validationFilename ), "validation dataset file (SVM-Light format)" )
		( "output,o", boost::program_options::value< std::string >( &output ), "output model file" )
		( "kernel,k", boost::program_options::value< std::string >( &kernelName ), "kernel" )
		( "kernel_parameter,K", boost::program_options::value< std::vector< std::string > >( &kernelParameters ), "kernel parameter(s)" )
		( "cache,c", boost::program_options::value< unsigned int >( &cacheSize )->default_value( 0 ), "size of kernel cache" )
		( "algorithm,a", boost::program_options::value< std::string >( &algorithmName ), "optimization algorithm" )
		( "algorithm_parameter,A", boost::program_options::value< std::vector< std::string > >( &algorithmParameters ), "algorithm parameter(s)" )
		( "biased,b", boost::program_options::value< bool >( &biased )->default_value( false ), "include an unregularized bias?" )
	;

	try {

		boost::program_options::variables_map variables;
		boost::program_options::store( boost::program_options::command_line_parser( argc, argv ).options( description ).run(), variables );
		boost::program_options::notify( variables );

		if ( variables.count( "help" ) ) {

			std::cout <<
				"Creates a model file from the given training and (optionally) validation" << std::endl <<
				"datasets (in SVM-Light format), and initializes it to the zero classifier--use" << std::endl <<
				"issvm_optimize to actually find an optimal classifier. The kernel parameter" << std::endl <<
				"must be one of \"linear\" or \"gaussian\", for which the kernel functions are:" << std::endl <<
				"\tlinear   ==>  K( x, y ) = <x,y>" << std::endl <<
				"\tgaussian ==>  K( x, y ) = exp( -p1 * || x - y ||^2 )" << std::endl <<
				"Here, \"p1\" is the value given for the first kernel parameter (-K option). For" << std::endl <<
				"classification, the algorithm parameter should be one of \"SBP\", \"SMO\" or" << std::endl <<
				"\"perceptron\", which take the parameters (-A option) nu, lambda and margin," << std::endl <<
				"respectively. For sparsification, the algorithm parameter should be" << std::endl <<
				"\"sparsifier\", which takes the parameters prediction_file, norm^2, and" << std::endl <<
				"optionally eta and suboptimality. The biased parameter selects whether the" << std::endl <<
				"optimization problem should include an unregularized bias." << std::endl <<
				std::endl <<
				description << std::endl;
		}
		else {

			if ( ! variables.count( "file" ) )
				throw std::runtime_error( "You must provide a training dataset file" );
			if ( ! variables.count( "output" ) )
				throw std::runtime_error( "You must provide an output file" );
			if ( ! variables.count( "kernel" ) )
				throw std::runtime_error( "You must provide a kernel parameter" );
			if ( ! variables.count( "cache" ) )
				throw std::runtime_error( "You must provide a cache parameter" );
			if ( ! variables.count( "algorithm" ) )
				throw std::runtime_error( "You must provide an algorithm parameter" );
			if ( ! variables.count( "biased" ) )
				throw std::runtime_error( "You must provide a bias parameter" );

			boost::shared_ptr< SVM::Optimizer::Base > pOptimizer(
				ConstructOptimizer(
					ConstructKernel< float >( trainingFilename, validationFilename, kernelName, kernelParameters, cacheSize ),
					algorithmName,
					algorithmParameters,
					biased
				)
			);

			//pOptimizer->WriteSupport("mySupport.txt");

			{	std::ofstream outputFile( output.c_str(), std::ios::out | std::ios::binary | std::ios::trunc );
				if ( outputFile.fail() )
					throw std::runtime_error( "Unable to open output file" );
				boost::iostreams::filtering_streambuf< boost::iostreams::output > outputStream;
				outputStream.push( outputFile );
				boost::archive::binary_oarchive outputArchive( outputStream );
				outputArchive << pOptimizer;
			}
		}
	}
	catch( std::exception& error ) {

		std::cerr << "Error: " << error.what() << std::endl << std::endl << description << std::endl;
		resultCode = EXIT_FAILURE;

	}

	return resultCode;
}
