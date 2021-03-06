#include <stdio.h>     // fprintf
#include <string.h>    // strcmp
#include <unistd.h>    // getopt
#include <stdlib.h>    // exit
#include <fstream>     // ofstream

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Controller.h"
#include "libsvm_interface.h"

using namespace std;


Controller::Controller() :
        number_of_data_types_(0)
    ,   files_()
    ,   data_type_names_()
    ,   reader_class_positive_set_()
    ,   reader_class_negative_set_()
    ,   reader_class_apply_()
    ,   model_output_file_("")
{
}





Controller::~Controller() {

    if (training_ || to_apply_) {

        svm_free_model_content(svm_);
        svm_free_and_destroy_model(&svm_);
    }
}





void Controller::parse_arguments(int argc, char* argv[]) {

    // return value of getopt
    int c;

    // argument booleans
    // file with log ratio values
    bool file_flag_wlog = false;
    // file with not-log ratio values
    bool file_flag_wolog = false;

    to_apply_ = false;
    // should training be made
    training_ = false;

    // flag for a file with mean values and single bin resolution
    // value gives the position in argument string under all files that have to be read
    int mean_file_flag = -1;

    bool model_given = false;

    bool directory_flag = false;
    // region file flag
    bool peak_flag = false;
    // negative training data region file flag
    bool negative_data_flag = false;

    bool model_output_flag = false;

    // matrix file flags
    bool positive_matrix_flag = false;
    bool negative_matrix_flag = false;


    // count the number of files or directorys
    // set flags and exit program if there are multiple files for the same operation
    for (int i = 0; i < argc; ++i) {

        // note: strcmp returns 0 if strings are equal
        if (!strcmp(argv[i], "-f") || !strcmp(argv[i], "-l")){

            ++number_of_data_types_;
        }

        if (!strcmp(argv[i], "-s")){

            mean_file_flag = ++number_of_data_types_;
        }

        if (!strcmp(argv[i], "-m")) {

            if (positive_matrix_flag || peak_flag) {

                fprintf(stderr, "Error: Multiple positive sets (-m -m  OR  -m -p  OR -m -b)\n");
                exit(EXIT_FAILURE);
            }
            positive_matrix_flag = true;
        }

        if (!strcmp(argv[i], "-t")) {

            if (negative_matrix_flag || negative_data_flag) {

                fprintf(stderr, "Error: Multiple negative sets (-t -t  OR  -t -n)\n");
                exit(EXIT_FAILURE);
            }
            negative_matrix_flag = true;
        }

        if (!strcmp(argv[i], "-n")) {

            if (negative_data_flag) {

                fprintf(stderr, "Error: Multiple negative region files given (-n)\n");
                exit(EXIT_FAILURE);
            }
            negative_data_flag = true;
        }
        if (!strcmp(argv[i], "-b") || (!strcmp(argv[i], "-p"))) {

            if (peak_flag) {

                fprintf(stderr, "Only one peak file (-p or -b) allowed.");
            } else {

                peak_flag = true;
            }
        }
    }

    training_ = peak_flag && (negative_matrix_flag || negative_data_flag);

    reader_class_apply_ = Reader(number_of_data_types_);

    reader_class_positive_set_ = Reader(number_of_data_types_);
    reader_class_negative_set_ = Reader(number_of_data_types_);


    // scan files and directory arguments
    while ((c = getopt(argc, argv, "x:a:m:t:n:p:b:l:f:d:o:s:")) != - 1) {

        switch (c) {

            // apply on this model
            case 'x':

                svm_ = svm_load_model(optarg);
                model_given = true;
                break;

            // apply on this set
            case 'a':

                reader_class_apply_.read_simplebed_file(optarg);
                to_apply_ = true;
                break;

            case 'n':

                reader_class_negative_set_.read_simplebed_file(optarg);
                break;

            // output filepath to save trained model
            case 'o':

                if (model_output_flag) {

                    fprintf(stderr, "Error: Multiple output pathes given (-o)\n");
                    exit(EXIT_FAILURE);
                }
                model_output_file_ = optarg;
                model_output_flag = true;
                break;

            case 'm':

                reader_class_positive_set_.read_matrix_file(optarg);
                break;

            case 't':

                reader_class_negative_set_.read_matrix_file(optarg);
                break;

            case 'p':

                reader_class_positive_set_.read_broadpeak_file(optarg);
                break;

            case 'b':

                reader_class_positive_set_.read_simplebed_file(optarg);
                break;

            case 's':

                file_flag_wolog = true;
                break;

            case 'l':

                file_flag_wlog = true;
                break;

            case 'f':

                file_flag_wolog = true;
                break;

            case 'd':

                directory_flag = true;
                break;

            case '?':

                if (optopt == 'x' || optopt == 'a' || optopt == 'm' || optopt == 't' || optopt == 'n' || optopt == 'b' || optopt == 'p' || optopt == 'f' || optopt == 'd' || optopt == 'o' || optopt == 'l') {

                    fprintf(stderr, "Missing argument for -%c option.\n", optopt);
                } else {

                    fprintf(stderr, "Unknown option: -%c\n", optopt);
                }
                exit(EXIT_FAILURE);

            default:

                fprintf(stderr, "An error occured while parsing the arguments.\n");
                exit(EXIT_FAILURE);
        }

        if (file_flag_wlog) {

            file_flag_wlog = false;
            files_.push_back(optarg);
            log_ratio_.push_back(true);

        } else if (file_flag_wolog) {

            file_flag_wolog = false;
            files_.push_back(optarg);
            log_ratio_.push_back(false);

        } else if (directory_flag) {

            // TODO: implement this with boost
            // directory_flag = false;
        }
    }

    // scan if names for the type of the files are provided
    data_type_names_ = vector<string>(number_of_data_types_, "<>");

    int i = optind;
    for (int counter = 0; i < argc; ++i, ++counter) {

        if (counter >= files_.size()) {

            fprintf(stderr, "Too many column descriptions provided.\n");
            exit(EXIT_FAILURE);
        }
        data_type_names_[counter] = argv[i];
    }

    int counter = 0;
    for (string file_name : files_) {

        // if no name is provided set column description to given file path
        if (data_type_names_[counter] == "<>") {

            data_type_names_[counter] = file_name;
        }
        ++counter;
    }


    vector<bool> log_buffer(begin(log_ratio_), end(log_ratio_));
    vector<string> files_buffer(begin(files_), end(files_));

    vector<bool> pos_avg_set;
    vector<bool> neg_avg_set;
    vector<bool> apply_avg_set;

    // do only if svm should be trained
    if (training_) {

        // read files
        // init matrix
        reader_class_positive_set_.init_matrix(number_of_data_types_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) num_threads(16)
#endif
        for (int counter = 0; counter < number_of_data_types_; ++counter) {

            if (mean_file_flag -1 != counter) {

                reader_class_positive_set_.read_file(files_buffer[counter], counter, log_buffer[counter]);

            } else {

                pos_avg_set = reader_class_positive_set_.read_mean_file(files_buffer[counter], counter);
            }
        }

        // read files again to build negative training samples
        if (negative_data_flag) {

            reader_class_negative_set_.init_matrix(number_of_data_types_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) num_threads(16)
#endif
            for (int counter = 0; counter < number_of_data_types_; ++counter) {

                if (mean_file_flag -1 != counter) {

                    reader_class_negative_set_.read_file(files_buffer[counter], counter, log_buffer[counter]);

                } else {

                    neg_avg_set = reader_class_negative_set_.read_mean_file(files_buffer[counter], counter);
                }
            }
        }

        // read files again to build negative training samples
        if (negative_data_flag) {

            reader_class_negative_set_.rescale_data();
            reader_class_negative_set_.collateBinnedRegions();
        }

        reader_class_positive_set_.rescale_data();
        reader_class_positive_set_.collateBinnedRegions();
    }


    // if the model should be applied read the data for the new regions
    if (to_apply_) {

        if (!(training_ || model_given)) {

            fprintf(stderr, "No model or training data provided!.\n");
            exit(EXIT_FAILURE);
        }

        reader_class_apply_.init_matrix(number_of_data_types_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) num_threads(16)
#endif
        for (int counter = 0; counter < number_of_data_types_; ++counter) {

            if (mean_file_flag -1 != counter) {

                reader_class_apply_.read_file(files_buffer[counter], counter, log_buffer[counter]);

            } else {

                apply_avg_set = reader_class_apply_.read_mean_file(files_buffer[counter], counter);
            }
        }
        reader_class_apply_.rescale_data();
        reader_class_apply_.collateBinnedRegions();
    }

    // Throw out lines with empty bins
    if (mean_file_flag != -1) {


        if (training_) {

            remove_zero_bins(pos_avg_set, reader_class_positive_set_.get_prev_read_data());

            if (negative_data_flag) {

                remove_zero_bins(neg_avg_set, reader_class_negative_set_.get_prev_read_data());

            }
        }

        if (to_apply_) {

            remove_zero_bins(apply_avg_set, reader_class_apply_.get_prev_read_data());

        }

    }

}






void Controller::print_prev_read_data() {

    ofstream os;
    os.open("positive_samples.matrix");
    os << "chrom\t" << "gen_start\t" << "gen_end\t";

    for (int i = 0; i < number_of_data_types_; ++i) {

        os << data_type_names_[i] << "\t";
    }
    os << endl;
    reader_class_positive_set_.print_prev_read_data(os);
    os.close();

    os.open("negative_samples.matrix");
    // wiggle file specific
    os << "gen\t" << "gen_start\t" << "gen_end\t";

    for (int i = 0; i < number_of_data_types_; ++i) {

        os << data_type_names_[i] << "\t";
    }
    os << endl;
    reader_class_negative_set_.print_prev_read_data(os);
    os.close();
}




void Controller::build_svm_model() {

    if (training_) {

        // translate matrix into libsvm problem
        struct svm_problem* prob = construct_svm_problem(reader_class_positive_set_.get_prev_read_data(), reader_class_negative_set_.get_prev_read_data());
        fprintf(stdout, "Successfully constructed a SVM problem out of the given data.\n");

        struct svm_problem* train_prob = (struct svm_problem*)(malloc(sizeof(*train_prob)));
        struct svm_problem* eval_prob = (struct svm_problem*)(malloc(sizeof(*eval_prob)));
        split_training_set(prob, train_prob, eval_prob);
        fprintf(stdout, "Splitted data into evaluation and training set.\n");

        struct svm_parameter* params = train_params(train_prob);
        fprintf(stdout, "Trained SVM parameters.\n");


        // train actual model
        svm_ = svm_train(train_prob, params);
        fprintf(stdout, "Trained final model.\n");

        FILE * result_out;
        result_out = fopen("result_summary.txt", "w+");

        pair<Matrix<int>,Matrix<int> > confusion = evaluate_model(train_prob, eval_prob, svm_);
        fprintf(result_out, "\n\nConfusion matrix test set:\n\n\
\t\t\t\t model outcome\n\
\t\t\t\tP\tN\n\
real outcome\tP\t%d\t%d\n\
\t\t\tN\t%d\t%d\n\n\
\n\nConfusion matrix train set:\n\
\t\t\t\t model outcome\n\
\t\t\t\tP\tN\n\
real outcome\tP\t%d\t%d\n\
\t\t\tN\t%d\t%d\n\n\
Number of Support Vectors: %d\n"
                , confusion.first(0,0), confusion.first(0,1), confusion.first(1,0), confusion.first(1,1)
                ,confusion.second(0,0), confusion.second(0,1), confusion.second(1,0), confusion.second(1,1)
                ,svm_get_nr_sv(svm_));
        fprintf(result_out, "\n\nBest parameters are:\n\n\
                C: %f\n\
                gamma: %f\n", params->C, params->gamma);
        fclose(result_out);
    }
}




void Controller::print_svm_model() {

    if (training_) {

        if (model_output_file_ != "") {

            if (svm_save_model(model_output_file_.c_str(), svm_) != 0) {

                fprintf(stderr, ("An error occured while saving the SVM model to " + model_output_file_ + "\n").c_str());
            }
        }
    }
}




void Controller::apply_svm_model() {

    // if there is no data to apply skip
    if (!to_apply_) {

        return;
    }


    fprintf(stdout, "Apply model.\n");

    Matrix<double> dummy = Matrix<double>();
    struct svm_problem* feature_matrix = construct_svm_problem(reader_class_apply_.get_prev_read_data(), dummy);

    vector<int> predicted_flags(feature_matrix->l, -1);
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) num_threads(20)
#endif
    for (int i = 0; i < feature_matrix->l; ++i) {

        predicted_flags[i] = svm_predict(svm_, feature_matrix->x[i]);
    }

    ofstream result_file("openChrom_call_flags");

    for (int i = 0; i < predicted_flags.size(); ++i) {

        if (predicted_flags[i] == 1) {

            result_file << "1\n";
        } else {

            result_file << "0\n";
        }
    }
    result_file.close();
}
