#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "Reader.h"
#include "svm.h"

using namespace std;

class Controller {

    public:

        Controller();

        ~Controller();

        // parse the given arguments
        //
        // @param: main method will delegate user arguments to this method
        //
        // WARNING: This method exits the whole program on wrong input
        void parse_arguments(int argc, char* argv[]);

        // print out the data that has been read
        void print_prev_read_data(ostream& os);


        // build a svm model
        // see also libsvm_interface.h
        void build_svm_model();

        // outputs the svm_model to model_output_file
        void print_svm_model();

    private:

        // number of data types that will be read
        int number_of_data_types_;

        // files to read
        list<string> files_;

        // corresponding names for the data type (e.g. "DNAse" "Histone_mod")
        vector<string> data_type_names_;

        // reads all files
        Reader reader_class_;

        // trained svm model
        struct svm_model* svm_;

        // contains a filepath for saving the svm model
        // see print_svm_model
        string model_output_file_;

};

#endif /* CONTROLLER_H */
