#ifndef CONTROLLER_READER_H
#define CONTROLLER_READER_H

#include <list>
#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "Reader.h"

using namespace std;

class Controller {

    public:

        Controller();

        ~Controller() = default;

        // parse the given arguments
        //
        // @param: main method will delegate user arguments to this method
        //
        // WARNING: This method exits the whole program on wrong input
        void parse_arguments(int argc, char* argv[]);

        // print out the data that has been read
        //
        // positive set: positive_samples.matrix
        // negative set: negative_samples.matrix
        void print_prev_read_data();


    private:

        // number of data types that will be read
        int number_of_data_types_;

        // files to read
        list<string> files_;
        // says if file (corresponding entry in files) contains log ratio (TRUE) or normal (FALSE) values
        list<bool> log_ratio_;

        // corresponding names for the data type (e.g. "DNAse" "Histone_mod")
        vector<string> data_type_names_;

        // Reader for the positive and negative training data set
        Reader reader_class_positive_set_;
        Reader reader_class_negative_set_;

};

#endif /* CONTROLLER_READER_H */
