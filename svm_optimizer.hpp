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
	\file svm_optimizer.hpp
	\brief Includes all SVM::Optimizer headers
*/




#ifndef __SVM_OPTIMIZER_HPP__
#define __SVM_OPTIMIZER_HPP__

#ifdef __cplusplus




/**
	\namespace SVM::Optimizer
	\brief SVM::Optimizer namespace
*/

/**
	\namespace SVM::Optimizer::Classification
	\brief SVM::Optimizer::Classification namespace
*/

/**
	\namespace SVM::Optimizer::Classification::Unbiased
	\brief SVM::Optimizer::Classification::Unbiased namespace
*/

/**
	\namespace SVM::Optimizer::Classification::Unbiased::_Private
	\brief SVM::Optimizer::Classification::Unbiased::_Private namespace
*/

/**
	\namespace SVM::Optimizer::Classification::Biased
	\brief SVM::Optimizer::Classification::Biased namespace
*/

/**
	\namespace SVM::Optimizer::Classification::Biased::_Private
	\brief SVM::Optimizer::Classification::Biased::_Private namespace
*/




#include "svm_optimizer_classification_unbiased_smo.hpp"
#include "svm_optimizer_classification_unbiased_sbp.hpp"
#include "svm_optimizer_classification_unbiased_perceptron.hpp"
#include "svm_optimizer_classification_unbiased_sparsifier.hpp"

#include "svm_optimizer_classification_biased_smo.hpp"
#include "svm_optimizer_classification_biased_sbp.hpp"
#include "svm_optimizer_classification_biased_perceptron.hpp"
#include "svm_optimizer_classification_biased_sparsifier.hpp"




#endif    /* __cplusplus */

#endif    /* __SVM_OPTIMIZER_HPP__ */
