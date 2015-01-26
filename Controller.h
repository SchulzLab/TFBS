#ifndef CONTROLLER_H
#define CONTROLLER_H

#include <stdio.h>
#include <iostream>
#include <iomanip>

#include "Reader.h"

using namespace std;

class Controller {

    public:

        Controller();

        // parse the given arguments
        //
        // main method will delegate user arguments to this method
        //
        // WARNING: This method exits the whole program on wrong input
        void parse_arguments(int argc, char* argv[]);

        // print out the data that has been read
        void print_prev_read_data(ostream& os);

    private:

        // number of data types that will be read
        int number_of_data_types_;

        // corresponding names for the data type (e.g. "DNAse" "Histone_mod")
        vector<string> data_type_names_;

        // reads all files
        Reader reader_class_;

};

#endif /* CONTROLLER_H */
