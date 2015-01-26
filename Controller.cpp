#include <stdio.h>     // fprintf
#include <unistd.h>    // getopt
#include <stdlib.h>    // exit

#include "Controller.h"

// ###################### MAIN #########################

int main(int argc, char* argv[]) {

    Controller controller;

    controller.parse_arguments(argc, argv);
    return 0;
}

// #####################################################





Controller::Controller() :
        number_of_data_types_(0)
    ,   data_type_names_()
    ,   reader_class_()
{
}





void Controller::parse_arguments(int argc, char* argv[]) {

    // return value of getopt
    int c;

    // argument booleans
    bool file_flag = false;
    bool directory_flag = false;

    // count the number of files or directorys
    while ((c = getopt(argc, argv, "f:d:")) != - 1) {


        switch (c) {

            case 'f':

                ++number_of_data_types_;
                break;

            case 'd':

                ++number_of_data_types_;
                break;

            case '?':

                if (optopt == 'f' || optopt == 'd') {

                    fprintf(stderr, "Missing argument for -%c option.\n", optopt);
                } else {

                    fprintf(stderr, "Unknown option: -%c\n", optopt);
                }
                exit(EXIT_FAILURE);

            default:

                fprintf(stderr, "An error occured while parsing the arguments.\n");
                exit(EXIT_FAILURE);
        }
    }

    reader_class_ = Reader(number_of_data_types_);

    // counter required for data type assignment for names
    int file_counter = 0;

    // scan files and direcorys arguments
    while ((c = getopt(argc, argv, "f:d:")) != - 1) {

        switch (c) {

            case 'f':

                file_flag = true;
                break;

            case 'd':

                directory_flag = true;
                break;

            case '?':

                if (optopt == 'f' || optopt == 'd') {

                    fprintf(stderr, "Missing argument for -%c option.\n", optopt);
                } else {

                    fprintf(stderr, "Unknown option: -%c\n", optopt);
                }
                exit(EXIT_FAILURE);

            default:

                fprintf(stderr, "An error occured while parsing the arguments.\n");
                exit(EXIT_FAILURE);
        }

        if (file_flag) {

            file_flag = false;
            reader_class_.read_file(optarg, file_counter);

        } else if (directory_flag) {

            // TODO: fix boost in Reader class
            // directory_flag = false;
            // reader_class_read_files_in_directory(optarg, file_counter);
        }

        ++file_counter;
    }

    // scan if names for the type of the files are provided
    data_type_names_ = vector<string>(number_of_data_types_, "<>");

    int i = optind;
    for (int counter = 0; i < argc; ++i, ++counter) {

        data_type_names_[c] = argv[i];
    }
}





void Controller::print_prev_read_data(ostream& os) {

    os << "Data that hast been read:\n";
    // wiggle file specific
    os << "gen_start" << "\t" << "gen_end" << "\t";

    for (int i = 0; i < number_of_data_types_; ++i) {

        os << data_type_names_[i] << "\t";
    }
    os << endl << endl;
    os << reader_class_.get_prev_read_data();
}
