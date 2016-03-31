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
	\file svm_optimizer_classification_private_find_water_level.hpp
	\brief SVM::Optimizer::Classification::Biased::_Private::FindWaterLevel and SVM::Optimizer::Classification::Unbiased::_Private::FindWaterLevel declarations
*/




#ifndef __SVM_OPTIMIZER_CLASSIFICATION_PRIVATE_FIND_WATER_LEVEL_HPP__
#define __SVM_OPTIMIZER_CLASSIFICATION_PRIVATE_FIND_WATER_LEVEL_HPP__

#ifdef __cplusplus




#include <utility>




namespace SVM {


namespace Optimizer {


namespace Classification {




namespace Unbiased {


namespace _Private {


//============================================================================
//    FindWaterLevel helper function
//============================================================================


double const FindWaterLevel( double* const begin, double* const end, double const total );


}    // namespace _Private


}    // namespace Unbiased




namespace Biased {


namespace _Private {


//============================================================================
//    FindWaterLevel helper function
//============================================================================


std::pair< double, double > const FindWaterLevel( double* const positiveBegin, double* const positiveEnd, double* const negativeBegin, double* const negativeEnd, double const total );


}    // namespace _Private


}    // namespace Biased




}    // namespace Classification


}    // namespace Optimizer


}    // namespace SVM




#endif    /* __cplusplus */

#endif    /* __SVM_OPTIMIZER_CLASSIFICATION_PRIVATE_FIND_WATER_LEVEL_HPP__ */
