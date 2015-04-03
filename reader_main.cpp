#include "Controller_reader.h"


// ###################### MAIN #########################




int main(int argc, char* argv[]) {

    Controller controller;

    controller.parse_arguments(argc, argv);
    controller.print_prev_read_data();

    return 0;
}




// #####################################################
