#ifndef LIBSVM_INTERFACE_H
#define LIBSVM_INTERFACE_H

#include "svm.h"
#include "Matrix.h"

// initialize a svm problem in syntax of libSVM using a
// positive and negative data set produced by Reader class
//
// @param: a positive and negative data set in Matrix format
//
// @return: svm_problem in libSVM struct
//
struct svm_problem* construct_svm_problem(Matrix<float>& positive_data, Matrix<float>& negative_data);


// initialize parameters for SVM in syntax of libSVM
//
// @return: svm_param in libSVM struct
struct svm_parameter* construct_svm_param();


#endif /* LIBSVM_INTERFACE_H */
