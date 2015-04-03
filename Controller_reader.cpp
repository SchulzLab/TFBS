#include <stdio.h>     // fprintf
#include <string.h>    // strcmp
#include <unistd.h>    // getopt
#include <stdlib.h>    // exit
#include <fstream>     // ofstream

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Controller.h"

using namespace std;


Controller::Controller() :
        number_of_data_types_(0)
    ,   files_()
    ,   data_type_names_()
    ,   reader_class_positive_set_()
    ,   reader_class_negative_set_()
    ,   model_output_file_("")
{
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
    // region file flag
    bool peak_flag = false;
    // negative training data region file flag
    bool negative_data_flag = false;

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

    reader_class_positive_set_ = Reader(number_of_data_types_);
    reader_class_negative_set_ = Reader(number_of_data_types_);

    if (!peak_flag) {

        fprintf(stderr, "No peak file (-p or -b) given.");
        exit(EXIT_FAILURE);
    }

    // scan files and directory arguments
    while ((c = getopt(argc, argv, "m:t:n:p:b:l:f:d:")) != - 1) {

        switch (c) {



            case 'n':

                reader_class_negative_set_.read_simplebed_file(optarg);
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

                if (optopt == 'm' || optopt == 't' || optopt == 'n' || optopt == 'b' || optopt == 'p' || optopt == 'f' || optopt == 'd' || optopt == 'l') {

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
    reader_class_positive_set_.init_matrix(number_of_data_types_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) num_threads(16)
#endif
    for (int counter = 0; counter < number_of_data_types_; ++counter) {

        reader_class_positive_set_.read_file(files_buffer[counter], counter, log_buffer[counter]);
    }

    // read files again to build negative training samples
    if (negative_data_flag) {

        reader_class_negative_set_.init_matrix(number_of_data_types_);
#ifdef _OPENMP
#pragma omp parallel for schedule(static,1) num_threads(16)
#endif
        for (int counter = 0; counter < number_of_data_types_; ++counter) {

            reader_class_negative_set_.read_file(files_buffer[counter], counter, log_buffer[counter]);
        }
    }
}





void Controller::print_prev_read_data() {

    ofstream os;
    os.open("positive_samples.matrix");
    os << "Data that has been read:\n";
    // wiggle file specific
    os << "gen\t" << "gen_start\t" << "gen_end\t";

    for (int i = 0; i < number_of_data_types_; ++i) {

        os << data_type_names_[i] << "\t";
    }
    os << endl << endl;
    reader_class_positive_set_.print_prev_read_data(os);
    os.close();

    os.open("negative_samples.matrix");
    os << "Data that has been read:\n";
    // wiggle file specific
    os << "gen\t" << "gen_start\t" << "gen_end\t";

    for (int i = 0; i < number_of_data_types_; ++i) {

        os << data_type_names_[i] << "\t";
    }
    os << endl << endl;
    reader_class_negative_set_.print_prev_read_data(os);
    os.close();
}
