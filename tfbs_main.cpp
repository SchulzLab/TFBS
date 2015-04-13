#include "Controller.h"


// ###################### MAIN #########################




int main(int argc, char* argv[]) {

    Controller controller;

    controller.parse_arguments(argc, argv);
    controller.print_prev_read_data();

    controller.build_svm_model();
    controller.print_svm_model();

    return 0;
}




// #####################################################
