#include <stdlib.h>
#include <random>
#include <functional> // bind
#include <list>

#include <fstream>

#include "libsvm_interface.h"


// part divided by 10 that should be used for the evaluation set
const float EVAL_SET_SIZE = 8;


struct svm_problem* construct_svm_problem(Matrix<double>& positive_data, Matrix<double>& negative_data) {


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

    // enable shrinking heuristics and disable probabilty estimates
    params->shrinking = 1;
    params->probability = 0;

    return params;
}





// number of folds for cross validation in train_params
// defined in libsvm_interface.cpp
constexpr int k_fold = 10;

struct svm_parameter* train_params(const struct svm_problem* prob) {

    // get initial parameter set
    struct svm_parameter* params;
    struct svm_parameter* best_params = construct_svm_param(0, 0);

    // this variable will contain the predicted labels for prob
    double* predicted;

    // init best_result worse than worst case (i.e. in every run is everything wrong + 1)
    int best_result = (k_fold * prob->l) + 1;

    ofstream os("log");
    os << "C_log10\tgamma\tdistance" << endl;

    // SEARCH GRID
    // train params in parallel
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) private(params, predicted) num_threads(16)
#endif
    for (int c = -6; c <= 5; ++c) {

        for (double gamma = pow(10,6); gamma >= pow(10,-6); gamma /= 10) {

            // omp init
            params = construct_svm_param(0, 0);
            predicted = (double*) malloc(sizeof(*predicted) * prob->l);

            params->C = pow(10, c);
            params->gamma = gamma;
            svm_cross_validation(prob, params, k_fold, predicted);

            // omp private
            int new_best = labelset_distance(prob, predicted);
            cerr << new_best << "\n";
            // if this is the best found set of parameters update the current best
#ifdef _OPENMP
#pragma omp critical
#endif
            {
            os << c << "\t" << gamma << "\t" << new_best << endl;
            if (new_best < best_result) {

                best_result = new_best;
                best_params->C = pow(10, c);
                best_params->gamma = gamma;
            }
            }

            free(predicted);

        }
    }


    // REFINED GRID FOR GAMMA

    double g = best_params->gamma;
    double fine_grid[] = {0.125, 0.25, 0.5, 2, 4, 8};
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) private(params, predicted) num_threads(16)
#endif
    for (int grid_index = 0; grid_index <= 5; ++grid_index) {

        // omp init
        params = construct_svm_param(0, 0);
        predicted = (double*) malloc(sizeof(*predicted) * prob->l);

        params->C = best_params->C;
        params->gamma = g * fine_grid[grid_index];
        svm_cross_validation(prob, params, k_fold, predicted);

        // omp private
        int new_best = labelset_distance(prob, predicted);
        cerr << new_best << "\n";
        // if this is the best found set of parameters update the current best
#ifdef _OPENMP
#pragma omp critical
#endif
        {
        os << best_params->C << "\t" << params->gamma << "\t" << new_best << endl;
        if (new_best < best_result) {

            best_result = new_best;
            best_params->gamma = params->gamma;
        }
        }

        free(predicted);

    }

    cout << "\n\nBest C: " << best_params->C << "\nBest gamma: " << best_params->gamma << "\nDistance: " << best_result << "\n\n";

    return best_params;
}





void split_training_set(const struct svm_problem* prob, struct svm_problem* training_set, struct svm_problem* eval_set) {


    // split the size to get 10% of the set as evaluation set
    eval_set->l = prob->l/10 * EVAL_SET_SIZE;
    training_set->l = prob->l - eval_set->l;
    cout  << "\nWhole set: " << prob->l << "\nEval set: " << eval_set->l << "\nTraining set: " << training_set->l << endl;
    // init sample flag memory
    eval_set->y = (double*)(malloc(eval_set->l * sizeof(*eval_set->y)));
    training_set->y = (double*)(malloc(training_set->l * sizeof(*training_set->y)));
    // init sample features memory
    eval_set->x = (struct svm_node **)(malloc(eval_set->l * sizeof(*eval_set->x)));
    training_set->x = (struct svm_node **)(malloc(training_set->l * sizeof(*training_set->x)));

    // save which samples has been drawn so far
    list<int> used_samples;

    // draw the evaluation set
    // draw positive samples
    // init random seed for mersenne twister and make drawn values uniformally distributed in problem size range
    mt19937::result_type seed = 100;
    auto dice_rand_pos = std::bind(std::uniform_int_distribution<int>(0, (prob->l/2) - 1), mt19937(seed));
    for (int i = 0; i < eval_set->l/2; ++i) {

        // draw a sample and test if it has been drawn before
        int new_sample;
        bool drawn_before;
        do {

            drawn_before = false;
            new_sample = dice_rand_pos();
            for (int used_sample : used_samples) {

                if (new_sample == used_sample) {

                    drawn_before = true;
                    break;
                }
            }

        } while(drawn_before);

        used_samples.push_back(new_sample);
        eval_set->y[i] = prob->y[new_sample];
        eval_set->x[i] = prob->x[new_sample];
    }

    auto dice_rand_neg = std::bind(std::uniform_int_distribution<int>(prob->l/2, prob->l - 1), mt19937(seed));
    for (int i = eval_set->l/2; i < eval_set->l; ++i) {

        // draw a sample and test if it has been drawn before
        int new_sample;
        bool drawn_before;
        do {

            drawn_before = false;
            new_sample = dice_rand_neg();
            for (int used_sample : used_samples) {

                if (new_sample == used_sample) {

                    drawn_before = true;
                    break;
                }
            }

        } while(drawn_before);

        used_samples.push_back(new_sample);
        eval_set->y[i] = prob->y[new_sample];
        eval_set->x[i] = prob->x[new_sample];
    }

    // assign the training set

    // index for initial set which we try to split
    int prob_index = 0;
    for (int i = 0; i < training_set->l; ++i) {

        bool used_in_eval;
        do {

            used_in_eval = false;
            for (int used_sample : used_samples) {

                if (prob_index == used_sample) {

                    used_in_eval = true;
                    ++prob_index;
                    break;
                }
            }

        } while(used_in_eval);

        training_set->y[i] = prob->y[prob_index];
        training_set->x[i] = prob->x[prob_index];
        ++prob_index;
    }
}





int labelset_distance(const struct svm_problem* prob, const double* predicted) {

    int distance = 0;

    // for (int fold = 0; fold < k_fold; ++fold) {
    //     for (int index = 0; index < prob->l; ++index) {
    //
    //         // compare flags assigned by the classifier
    //         if (prob->y[index] != predicted[index + fold * prob->l]) {
    //
    //             ++distance;
    //         }
    //     }
    // }
    for (int index = 0; index < prob->l; ++index) {

        if (prob->y[index] != predicted[index]) {

            ++distance;
        }

    }

    return distance;
}





pair<Matrix<int>,Matrix<int> > evaluate_model(const struct svm_problem* train_prob, const struct svm_problem* eval_prob, const svm_model* model) {


    // confusion matrix
    //                          model outcome
    //                          P    N
    //     real outcome    P    TP   FN
    //                     N    FP   TN
    Matrix<int> confusion_test (2,2,0);
    Matrix<int> confusion_train (2,2,0);

    for (int i = 0; i < eval_prob->l; ++i) {

        double pred = svm_predict(model, eval_prob->x[i]);
        double outcome = eval_prob->y[i];
        if (pred == 1) {

            pred == outcome ? ++confusion_test(0,0) : ++confusion_test(1,0);
        } else {

            pred == outcome ? ++confusion_test(1,1) : ++confusion_test(0,1);
        }
    }


    for (int i = 0; i < train_prob->l; ++i) {

        double pred = svm_predict(model, train_prob->x[i]);
        double outcome = train_prob->y[i];
        if (pred == 1) {

            pred == outcome ? ++confusion_train(0,0) : ++confusion_train(1,0);
        } else {

            pred == outcome ? ++confusion_train(1,1) : ++confusion_train(0,1);
        }
    }
    return pair<Matrix<int>, Matrix<int> > (confusion_test, confusion_train);
}
