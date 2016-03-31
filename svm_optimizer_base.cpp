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
	\file svm_optimizer_base.cpp
	\brief SVM::Optimizer::Base implementation
*/




#include "svm_optimizer_base.hpp"




namespace SVM {


namespace Optimizer {




//============================================================================
//    Base methods
//============================================================================


Base::~Base() {
}


double const Base::Bias() const {

	return 0;
}


unsigned int const Base::Support() const {

	unsigned int const size = this->TrainingSize();
	boost::shared_array< double > alphas( new double[ size ] );
	this->GetAlphas( alphas.get(), alphas.get() + size );

	unsigned int support = 0;
	{	double const* ii    = alphas.get();
		double const* iiEnd = ii + size;
		for ( ; ii != iiEnd; ++ii )
			if ( *ii != 0 )
				++support;
	}

	return support;
}

void Base::WriteSupport(std::string filename){

		std::cout << "Calling Get Supp of Optimizer Base" << std::endl;
		std::cout << "This is not doing something useful" << std::endl;
		std::cout << "You need to create an specialized function in the specialized Optimizer Class" << std::endl;

}

double const Base::ValidationError() const {

	unsigned int const size = this->ValidationSize();
	boost::shared_array< double > responses( new double[ size ] );
	this->GetValidationResponses( responses.get(), responses.get() + size );

	double error = 0;
	{	double const* ii    = responses.get();
		double const* iiEnd = ii + size;
		for ( ; ii != iiEnd; ++ii )
			if ( *ii <= 0 )
				++error;
	}

	return( error / size );
}




}    // namespace Optimizer


}    // namespace SVM
