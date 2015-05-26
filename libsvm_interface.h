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
struct svm_problem* construct_svm_problem(Matrix<double>& positive_data, Matrix<double>& negative_data);


// initialize parameters for SVM in syntax of libSVM
//
// @param: rbf_gamma - kernel parameter (we use a RBF kernel)
//         svm_c     - slack variable
//
// @return: svm_param in libSVM struct
//
struct svm_parameter* construct_svm_param(double rbf_gamma, double svm_c);


// train best set of parameters by using k-fold cross validation
//
// @param: prob - problem to solve
//
// @return: the best found parameter set
//
struct svm_parameter* train_params(const struct svm_problem* prob);


// get distance between a given problem and (by cross validation) predicted labels
// if predicted label differ from real label distance is increased by 1
//
// @param: prob - libSVM problem on which cross validation was executed
//         predicted - cross validation output
//
// @return: computed distance
//
int labelset_distance(const struct svm_problem* prob, const double* predicted);


// split given libsvm problem into 2 different sets (10% will be drawn randomly for the evaluation set
// and the rest will be assigned to the training set)
//
// @param:  prob - the libsvm problem (data set) that should be splitted
//          training_set - empty svm_problem that will hold the training set
//          eval_set - empty svm_problem that will hold the evaluation set
//
void split_training_set(const struct svm_problem* prob, struct svm_problem* training_set, struct svm_problem* eval_set);


// evaluate a trained model on a given sample set
//
// @param:  eval_prob - the sample set on which the model should be evaluated on
//          model - the trained model that should be tested
//
// @return: the number of falsely classified samples
//
pair<int, int> evaluate_model(const struct svm_problem* train_prob, const struct svm_problem* eval_prob, const svm_model* model);


#endif /* LIBSVM_INTERFACE_H */
