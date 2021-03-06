#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <list>
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
        //
        // positive set: positive_samples.matrix
        // negative set: negative_samples.matrix
        void print_prev_read_data();


        // build a svm model
        // see also libsvm_interface.h
        void build_svm_model();

        // outputs the svm_model to model_output_file
        void print_svm_model();

        // apply the model to data saved in reader_class_apply_
        //
        // called open Chromatin regions will be saved as flags congruent to the given file in openChrom_call_flags
        void apply_svm_model();

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
        // Reader for the set that the model should be applied to
        // and the flag if there was a file provided for the appliance
        Reader reader_class_apply_;
        bool to_apply_;
        // flag if training should be done
        bool training_;

        // trained svm model
        struct svm_model* svm_;

        // contains a filepath for saving the svm model
        // see print_svm_model
        string model_output_file_;

};

#endif /* CONTROLLER_H */
