
#include <stdlib.h>
#include <vector>

#include "libsvm_interface.h"

struct svm_problem* construct_svm_problem(Matrix<float>& positive_data, Matrix<float>& negative_data) {

    struct svm_problem* svm_prob = (svm_problem*)malloc(sizeof(*svm_prob));

    // initialize number of training data samples
    svm_prob->l = positive_data.get_number_of_lines() + negative_data.get_number_of_lines();

    // initialize an array containing the labels for all training data samples
    vector<double> data_labels = vector<double>(svm_prob->l);
    svm_prob->y = &data_labels[0];

    int i_positive;
    for (i_positive = 0; i_positive < positive_data.get_number_of_lines(); ++i_positive) {

        data_labels[i_positive] = 1;
    }
    for (int i_negative = i_positive + 1; i_negative < svm_prob->l; ++i_negative) {

        data_labels[i_negative] = -1;
    }

    // convert the training data in libSVM format
    vector<struct svm_node*> training_set = vector<struct svm_node*>(svm_prob->l);

    //      positive samples
    for (int i_sample = 0; i_sample < positive_data.get_number_of_lines(); ++i_sample) {

        vector<struct svm_node> training_sample;
        training_sample.reserve(positive_data.get_number_of_columns() - 2);

        // the first 3 columns must be skipped, as they only contain chromosome region information
        for (int i_feature = 3; i_feature < positive_data.get_number_of_columns(); ++i_feature) {

            if (positive_data(i_sample, i_feature) != 0) {

                struct svm_node* sample_feature = (svm_node*)malloc(sizeof(*sample_feature));
                sample_feature->index = i_feature - 2;
                sample_feature->value = positive_data(i_sample, i_feature);
                training_sample.push_back(*sample_feature);
            }
        }
        struct svm_node* fin_flag = (svm_node*)malloc(sizeof(*fin_flag));
        fin_flag->index = -1;
        training_sample.push_back(*fin_flag);

        training_set[i_sample] = &training_sample[0];
    }

    //      negative samples
    for (int i_sample = 0; i_sample < negative_data.get_number_of_lines(); ++i_sample) {

        vector<struct svm_node> training_sample;
        training_sample.reserve(negative_data.get_number_of_columns() - 2);

        // the first 3 columns must be skipped, as they only contain chromosome region information
        for (int i_feature = 3; i_feature < negative_data.get_number_of_columns(); ++i_feature) {

            if (negative_data(i_sample, i_feature) != 0) {

                struct svm_node* sample_feature = (svm_node*)malloc(sizeof(*sample_feature));
                sample_feature->index = i_feature - 2;
                sample_feature->value = negative_data(i_sample, i_feature);
                training_sample.push_back(*sample_feature);
            }
        }
        struct svm_node* fin_flag = (svm_node*)malloc(sizeof(*fin_flag));
        fin_flag->index = -1;
        training_sample.push_back(*fin_flag);

        training_set[i_sample + positive_data.get_number_of_lines()] = &training_sample[0];
    }



    svm_prob->x = &training_set[0];
    return svm_prob;
}




struct svm_parameter* construct_svm_param() {

    struct svm_parameter* params = (svm_parameter*)malloc(sizeof(*params));

    params->svm_type = C_SVC;
    params->kernel_type = RBF;

    params->cache_size = 0 /* TODO */;

    // recommended stopping criterion by libSVM
    params->eps = 0.01;

    params->C = 0 /* TODO */;

    // Do we need this
    params->nr_weight = 0;


    return params;
}
