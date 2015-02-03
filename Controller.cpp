#include <stdio.h>     // fprintf
#include <string.h>    // strcmp
#include <unistd.h>    // getopt
#include <stdlib.h>    // exit

#include "Controller.h"

// ###################### MAIN #########################

int main(int argc, char* argv[]) {

    Controller controller;

    controller.parse_arguments(argc, argv);
    controller.print_prev_read_data(cout);
    return 0;
}

// #####################################################





Controller::Controller() :
        number_of_data_types_(0)
    ,   files_()
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
    bool peak_flag = false;

    // count the number of files or directorys
    for (int i = 0; i < argc; ++i) {

        // note: strcmp returns 0 if strings are equal
        !strcmp(argv[i], "-f") || !strcmp(argv[i], "-d") ? ++number_of_data_types_ : 0;

        if (!strcmp(argv[i], "-p")) {

            if (peak_flag) {

                fprintf(stderr, "Only one peak file (-p) allowed.");
            } else {

                peak_flag = true;
            }
        }
    }

    reader_class_ = Reader(number_of_data_types_);

    // scan files and directory arguments
    while ((c = getopt(argc, argv, "p:f:d:")) != - 1) {

        switch (c) {

            case 'p':

                reader_class_.read_peak_file(optarg);
                break;

            case 'f':

                file_flag = true;
                break;

            case 'd':

                directory_flag = true;
                break;

            case '?':

                if (optopt == 'p' || optopt == 'f' || optopt == 'd') {

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
            files_.push_back(optarg);

        } else if (directory_flag) {

            // TODO: fix boost in Reader class
            // directory_flag = false;
        }
    }

    // scan if names for the type of the files are provided
    data_type_names_ = vector<string>(number_of_data_types_, "<>");

    int i = optind;
    for (int counter = 0; i < argc; ++i, ++counter) {

        data_type_names_[counter] = argv[i];
    }
    for (int counter = 0; counter < number_of_data_types_; ++counter) {

        reader_class_.read_file(files_.front(), counter);
        files_.pop_front();
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
    os << reader_class_.get_prev_read_data();
}
