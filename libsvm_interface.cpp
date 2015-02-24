#include <stdlib.h>

#include "libsvm_interface.h"




struct svm_problem* construct_svm_problem(Matrix<float>& positive_data, Matrix<float>& negative_data) {

    struct svm_problem* svm_prob = (svm_problem*)malloc(sizeof(*svm_prob));

    // initialize number of training data samples
    svm_prob->l = positive_data.get_number_of_lines() + negative_data.get_number_of_lines();

    // initialize an array containing the labels for all training data samples
    double* data_labels = (double*)malloc(sizeof(*data_labels) * svm_prob->l);
    svm_prob->y = data_labels;

    int i_positive;
    for (i_positive = 0; i_positive < positive_data.get_number_of_lines(); ++i_positive) {

        data_labels[i_positive] = 1;
    }
    for (int i_negative = i_positive + 1; i_negative < svm_prob->l; ++i_negative) {

        data_labels[i_negative] = -1;
    }

    // convert the training data in libSVM format
    struct svm_node** training_set = (struct svm_node**)malloc(sizeof(*training_set) * svm_prob->l);

    //      positive samples
    for (int i_sample = 0; i_sample < positive_data.get_number_of_lines(); ++i_sample) {

        struct svm_node* training_sample = (struct svm_node*)malloc(sizeof(*training_sample) * positive_data.get_number_of_columns() - 2);

        int sample_index = 0;
        // the first 3 columns must be skipped, as they only contain chromosome region information
        for (int i_feature = 3; i_feature < positive_data.get_number_of_columns(); ++i_feature) {

            if (positive_data(i_sample, i_feature) != 0) {

                struct svm_node* sample_feature = (svm_node*)malloc(sizeof(*sample_feature));

                // note: libsvm start indexing at one and not zero
                sample_feature->index = i_feature - 2;
                sample_feature->value = positive_data(i_sample, i_feature);
                training_sample[sample_index++] = *sample_feature;
            }
        }
        struct svm_node* fin_flag = (svm_node*)malloc(sizeof(*fin_flag));
        fin_flag->index = -1;
        training_sample[sample_index] = *fin_flag;

        training_set[i_sample] = training_sample;
    }

    //      negative samples
    for (int i_sample = 0; i_sample < negative_data.get_number_of_lines(); ++i_sample) {

        struct svm_node* training_sample = (struct svm_node*)malloc(sizeof(*training_sample) * positive_data.get_number_of_columns() - 2);

        int sample_index = 0;
        // the first 3 columns must be skipped, as they only contain chromosome region information
        for (int i_feature = 3; i_feature < negative_data.get_number_of_columns(); ++i_feature) {

            if (negative_data(i_sample, i_feature) != 0) {

                struct svm_node* sample_feature = (svm_node*)malloc(sizeof(*sample_feature));

                // note: libsvm start indexing at one and not zero
                sample_feature->index = i_feature - 2;
                sample_feature->value = negative_data(i_sample, i_feature);
                training_sample[sample_index++] = *sample_feature;
            }
        }
        struct svm_node* fin_flag = (svm_node*)malloc(sizeof(*fin_flag));
        fin_flag->index = -1;
        training_sample[sample_index] = *fin_flag;

        training_set[i_sample + positive_data.get_number_of_lines()] = training_sample;
    }



    svm_prob->x = &training_set[0];
    for (int sample = 0; sample < svm_prob->l; ++sample) {

        int i = 0;
        while (svm_prob->x[sample][i].index != -1) {

            cerr << "(" << svm_prob->x[sample][i].index << "/" << svm_prob->x[sample][i].value << ")\n\n";
            ++i;
        }

    }
    return svm_prob;
}





struct svm_parameter* construct_svm_param(double rbf_gamma, double svm_c) {

    struct svm_parameter* params = (svm_parameter*)malloc(sizeof(*params));

    params->svm_type = C_SVC;
    params->C = svm_c;

    params->kernel_type = RBF;
    params->gamma = rbf_gamma;

    // cache size in MB
    params->cache_size = 50000;

    // recommended stopping criterion by libSVM
    params->eps = 0.01;

    // disable weighted data
    params->nr_weight = 0;

    // disable shrinking heuristics and probabilty estimates
    params->shrinking = 0;
    params->probability = 0;

    return params;
}





// number of folds for cross validation in train_params
// defined in libsvm_interface.cpp
constexpr int k_fold = 10;

struct svm_parameter* train_params(const struct svm_problem* prob) {

    // get initial parameter set
    // OMP: private
    struct svm_parameter* params;
    // OMP: shared
    struct svm_parameter* best_params = construct_svm_param(0, 0);

    // OMP: private
    // this variable will contain the predicted labels for prob
    double* predicted;

    // init best_result worse than worst case (i.e. in every run is everything wrong + 1)
    int best_result = (k_fold * prob->l) + 1;


    // train params in parallel
    // outer loop: slack variable c 0% (hard margin) to 5% of training data
    for (int c = 1; c <= prob->l /*/ 20*/; c++) {

        // inner loop: rbf kernel parameter from 1 to 10^(-6)
        for (double gamma = 1; gamma >= 0.000001; gamma /= 10) {

            // omp init
            params = construct_svm_param(0, 0);
            predicted = (double*) malloc(sizeof(*predicted) * prob->l * k_fold);

            params->C = c;
            params->gamma = gamma;
            svm_cross_validation(prob, params, k_fold, predicted);

            // omp private
            int new_best;
            // if this is the best found set of parameters update the current best
// #ifdef _OPNEMP
// #pragma omp atomic
// #endif
            if ((new_best = labelset_distance(prob, predicted)) < best_result) {

                best_result = new_best;
                best_params->C = c;
                best_params->gamma = gamma;
            }

            free(predicted);

        }
    }

    return best_params;
}





int labelset_distance(const struct svm_problem* prob, const double* predicted) {

    int distance = 0;

    // for (int fold = 0; fold  < k_fold; ++fold) {

        for (int index = 0; index < prob->l; ++index) {

            // libSVM seems to fill predicted with probabilitys for c_svm - change comparison?
            if (prob->y[index] == predicted[/*fold * prob->l +*/ index]) {

                ++distance;
            }
        }
    // }
    return distance;
}
