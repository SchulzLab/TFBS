#include <stdio.h>     // fprintf
#include <string.h>    // strcmp
#include <unistd.h>    // getopt
#include <stdlib.h>    // exit

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Controller.h"
#include "libsvm_interface.h"

using namespace std;


// ###################### MAIN #########################

int main(int argc, char* argv[]) {

    Controller controller;

    controller.parse_arguments(argc, argv);
    controller.print_prev_read_data(cout);

    // controller.build_svm_model();

    return 0;
}

// #####################################################





Controller::Controller() :
        number_of_data_types_(0)
    ,   files_()
    ,   data_type_names_()
    ,   reader_class_()
    ,   model_output_file_("")
{
}





Controller::~Controller() {

    // svm_free_model_content(svm_);
    // svm_free_and_destroy_model(&svm_);
}





void Controller::parse_arguments(int argc, char* argv[]) {

    // return value of getopt
    int c;

    // argument booleans
    // file with log ratio values
    bool file_flag_wlog = false;
    // file with not-log ratio values
    bool file_flag_wolog = false;
    bool directory_flag = false;
    bool peak_flag = false;

    // count the number of files or directorys
    for (int i = 0; i < argc; ++i) {

        // note: strcmp returns 0 if strings are equal
        !strcmp(argv[i], "-f") || !strcmp(argv[i], "-l") ? ++number_of_data_types_ : 0;

        if (!strcmp(argv[i], "-b") || (!strcmp(argv[i], "-p"))) {

            if (peak_flag) {

                fprintf(stderr, "Only one peak file (-p or -b) allowed.");
            } else {

                peak_flag = true;
            }
        }
    }

    reader_class_ = Reader(number_of_data_types_);

    if (!peak_flag) {

        fprintf(stderr, "No peak file (-p or -b) given.");
        exit(EXIT_FAILURE);
    }

    bool model_output_flag = false;
    // scan files and directory arguments
    while ((c = getopt(argc, argv, "p:b:l:f:d:o:")) != - 1) {

        switch (c) {

            case 'o':

                if (model_output_flag) {

                    fprintf(stderr, "Error: Multiple output pathes given (-o)\n");
                    exit(EXIT_FAILURE);
                }
                model_output_file_ = optarg;
                model_output_flag = true;
                break;

            case 'p':

                reader_class_.read_broadpeak_file(optarg);
                break;

            case 'b':

                reader_class_.read_simplebed_file(optarg);
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

                if (optopt == 'b' || optopt == 'p' || optopt == 'f' || optopt == 'd' || optopt == 'o' || optopt == 'l') {

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

    // read files
    // init matrix
    reader_class_.init_matrix(number_of_data_types_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) num_threads(16)
#endif
    for (int counter = 0; counter < number_of_data_types_; ++counter) {

        reader_class_.read_file(files_buffer[counter], counter, log_buffer[counter]);
    }
}





void Controller::print_prev_read_data(ostream& os) {

    os << "Data that has been read:\n";
    // wiggle file specific
    os << "gen\t" << "gen_start\t" << "gen_end\t";

    for (int i = 0; i < number_of_data_types_; ++i) {

        os << data_type_names_[i] << "\t";
    }
    os << endl << endl;
    reader_class_.print_prev_read_data(os);
}




void Controller::build_svm_model() {

    Matrix<float> m = Matrix<float>();
    struct svm_problem* prob = construct_svm_problem(reader_class_.get_prev_read_data(), m);

    struct svm_parameter* params = train_params(prob);

    // check parameters doesn't make sense since you must init all values of params which is cluelsess
    // const char* errormessage = svm_check_parameter(prob, params);
    //
    // if (errormessage != NULL) {
    //
    //     fprintf(stderr, errormessage);
    // }

    // train actual model
    svm_ = svm_train(prob, params);
}




void Controller::print_svm_model() {

    if (model_output_file_ != "") {

        if (svm_save_model(model_output_file_.c_str(), svm_) != 0) {

            fprintf(stderr, ("An error occured while saving the SVM model to " + model_output_file_ + "\n").c_str());
        }
    }
}
