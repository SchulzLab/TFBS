// standard lib
#include <iostream>
#include <fstream>
#include <math.h>  // abs() exp()

// external lib
// #include "boost/filesystem.hpp"

// project
#include "Matrix.h"
#include "Reader.h"

using namespace std;

// namespace bf = boost::filesystem;

Reader::Reader() :
        matrix_()
    ,   number_of_data_types_(0)
    ,   line_counter_(0)
{
}





Reader::Reader(int number_of_data_types) :
        matrix_()
    ,   number_of_data_types_(number_of_data_types)
    ,   line_counter_(0)
    ,   chrom_numerical_(0)
{
}





Reader::Reader(Reader&& other_reader) :
        matrix_(move(other_reader.matrix_))
    ,   number_of_data_types_(other_reader.number_of_data_types_)
    ,   map_str_to_chr_(other_reader.map_str_to_chr_)
    ,   map_chr_to_str_(other_reader.map_chr_to_str_)
    ,   chrom_numerical_(other_reader.chrom_numerical_)
{
    other_reader.matrix_ = Matrix<float>();
    other_reader.map_str_to_chr_ = unordered_map<string, int>();
    other_reader.map_chr_to_str_ = unordered_map<int, string>();
}





Reader& Reader::operator=(Reader&& other_reader){

    matrix_ = move(other_reader.matrix_);
    number_of_data_types_ = other_reader.number_of_data_types_;
    map_str_to_chr_ = move(other_reader.map_str_to_chr_);
    map_chr_to_str_ = move(other_reader.map_chr_to_str_);
    chrom_numerical_ = other_reader.chrom_numerical_;

    return *this;
}





// void Reader::read_files_in_directory(const string& directory_path) {
//
//     bf::path directory(directory_path);
//     bf::directory_iterator iterator_end;
//
//     // tests if path ist valid
//     if (bf::exists(directory) && bf::is_directory(directory)) {
//
//         // iterate over all files contained in the directory
//         for (bf::directory_iterator dir_iterator(directory);
//                 dir_iterator != iterator_end;
//                 ++dir_iterator) {
//
//             // perform file reading on all valid files
//             if (bf::is_regular_file(dir_iterator->status())) {
//
//                 read_file(dir_iterator->path().string());
//             }
//         }
//     } else {
//
//         cerr << "Not able to read files - no such directory: "
//              << directory_path << endl;
//     }
// }





void Reader::read_file(const string& file_path, int data_type, bool is_log) {

    FILE* actual_file = fopen (file_path.c_str(), "r");


    // variables holding the content of one line of the file temporarily
    int chrom_begin, chrom_end, chromosome;
    char chrom[50];
    float peak;

    vector<float> new_feature(line_counter_, .0);

    matrix_.append_new_column(new_feature);

    // read one line in the file per loop
    while (fscanf(actual_file, "%s %d %d %f", chrom, &chrom_begin, &chrom_end, &peak) == 4) {

        if (map_str_to_chr_.find(chrom) != map_str_to_chr_.end()) {

            chromosome = map_str_to_chr_[chrom];

            // store just read data in matrix
            binary_search(chromosome, chrom_begin, chrom_end, 0, line_counter_ - 1, data_type, peak, is_log);

        }
    }

    if (!feof(actual_file)) {

        fprintf(stderr, ("\nA reading error occured while reading \"" + file_path  + "\"\n").c_str());
    }

    fclose(actual_file);
}





void Reader::read_peak_file(const string& file_path) {


    // adapt this for faster peak file reading
    // it should be bigger then the actual number of lines in the file
    constexpr int INITIAL_V_CAP = 70000;
    constexpr int V_CAP_STEPSIZE = 10000;

    // actual capabilty of one column
    int actually_reserved = INITIAL_V_CAP;

    FILE* peak_file = fopen (file_path.c_str(), "r");

    char chrom[50];
    int chromosome, chrom_begin, chrom_end;

    // initialize and reserve memory for peak file values
    vector<float> chromosome_v;
    chromosome_v.reserve(INITIAL_V_CAP);
    vector<float> chrom_begin_v;
    chrom_begin_v.reserve(INITIAL_V_CAP);
    vector<float> chrom_end_v;
    chrom_end_v.reserve(INITIAL_V_CAP);

    matrix_.append_new_column(chromosome_v);
    matrix_.append_new_column(chrom_begin_v);
    matrix_.append_new_column(chrom_end_v);

    // read one line in the file per loop
    // skip the additional information of broadpeak format
    while (fscanf(peak_file, "%s %d %d %*s %*s %*s %*s %*s %*s", chrom, &chrom_begin, &chrom_end) == 3) {


        ++line_counter_;
        if (map_str_to_chr_.find(chrom) != map_str_to_chr_.end()) {

            chromosome = map_str_to_chr_[chrom];

        } else {

            map_str_to_chr_[chrom] = chrom_numerical_;
            map_chr_to_str_[chrom_numerical_] = chrom;
            chromosome = chrom_numerical_++;
        }

        // if necessary allocate more memory for matrix columns
        if (line_counter_ > actually_reserved) {

            matrix_.get_column(0).reserve(actually_reserved + V_CAP_STEPSIZE);
            matrix_.get_column(1).reserve(actually_reserved + V_CAP_STEPSIZE);
            matrix_.get_column(2).reserve(actually_reserved + V_CAP_STEPSIZE);
            actually_reserved += V_CAP_STEPSIZE;
        }

        matrix_.get_column(0).push_back(chromosome);
        matrix_.get_column(1).push_back(chrom_begin);
        matrix_.get_column(2).push_back(chrom_end);

    }

    if (!feof(peak_file)) {

        fprintf(stderr, ("\nA reading error occured while reading peak file \"" + file_path  + "\"\n").c_str());
    }

    fclose(peak_file);
}





void Reader::print_prev_read_data(ostream& os) {


    os << endl << endl;
    if (matrix_.get_number_of_lines() == 0) {

        os << "<Empty>" << endl;
    } else {

        if (matrix_.get_number_of_columns() == 0) {

            os << "<Empty>" << endl;
        }
    }

    // iterate over lines
    for (int lines = 0; lines < matrix_.get_number_of_lines(); ++lines) {


        os << setw(10) << right << map_chr_to_str_[matrix_(lines, 0)] << "\t";
        // iterate over columns
        for (int columns = 1; columns < matrix_.get_number_of_columns(); ++columns) {

            os << setw(5) << right << matrix_(lines, columns) << "\t";
        }
        os << "\n";
    }
}





Matrix<float>& Reader::get_prev_read_data() {

    return matrix_;
}





void Reader::binary_search(const int chrom, const int chrom_begin, const int chrom_end, const int start_point, const int end_point, const int data_type, const float peak, bool is_log) {


    // termination
    if (start_point > end_point) {

        return;
    }

    // center of one vector
    const int central_point = start_point + ((end_point - start_point) / 2);

    // tests if queue is on the same chromosome as new sample
    if (matrix_(central_point, 0) == chrom) {

        // region overlap tests between queue and sample
        if (matrix_(central_point, 1) >= chrom_begin) {

            if (matrix_(central_point, 1) == chrom_begin) {

                // sample:  |---|
                // queue:   |---|
                //  or
                // sample:  |--|
                // queue:   |---|
                if (matrix_(central_point, 2) >= chrom_end) {

                    matrix_(central_point, data_type + 3) += peak;

                // sample: |----|
                // queue:  |---|
                } else {

                    // remaining peak value
                    float partial_peak_buffer;
                    float partial_peak;

                    if (is_log) {

                        // peak value that overlaps with the queue
                        partial_peak = (float)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * pow(2, peak);
                        matrix_(central_point, data_type + 3) += partial_peak;
                        partial_peak_buffer = pow(2, peak) - partial_peak;

                    } else {

                        partial_peak = (float)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * peak;
                        matrix_(central_point, data_type + 3) += partial_peak;
                        partial_peak_buffer = peak - partial_peak;

                    }


                    // increase position in each step by one
                    int i = 1;
                    int mod_start = matrix_(central_point, 2);

                    while (matrix_(central_point + i, 0) == chrom && matrix_(central_point + i, 1) < chrom_end) {

                        if (matrix_(central_point + i, 2) >= chrom_end) {

                            partial_peak = (float)(chrom_end - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
                            matrix_(central_point, data_type + 3) += partial_peak;
                            break;

                        } else {

                            partial_peak = (float)(matrix_(central_point + i, 2) - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
                            matrix_(central_point + i, data_type + 3) += partial_peak;
                            partial_peak_buffer -= partial_peak;
                            mod_start = matrix_(central_point + i, 2);
                            ++i;
                        }
                    }
                }

            // matrix(central_point, 1) > chrom_begin
            } else {

                if (matrix_(central_point, 1) < chrom_end) {

                    // sample:  |---|
                    // queue:    |-|
                    if (matrix_(central_point, 2) < chrom_end) {

                        // remaining peak value
                        float partial_peak_buffer;

                        if (is_log) {

                            // peak value that overlaps with the queue
                            const float partial_peak = (float)(matrix_(central_point, 2) - matrix_(central_point, 1))/(chrom_end - chrom_begin) * pow(2, peak);
                            matrix_(central_point, data_type + 3) += partial_peak;

                            partial_peak_buffer = (float)(matrix_(central_point, 1) - chrom_begin)/(chrom_end - chrom_begin) * pow(2, peak);

                        } else {

                            // peak value that overlaps with the queue
                            const float partial_peak = (float)(matrix_(central_point, 2) - matrix_(central_point, 1))/(chrom_end - chrom_begin) * peak;
                            matrix_(central_point, data_type + 3) += partial_peak;

                            partial_peak_buffer = (float)(matrix_(central_point, 1) - chrom_begin)/(chrom_end - chrom_begin) * peak;
                        }

                        int i = 1;
                        int mod_end = matrix_(central_point, 1);

                        while (matrix_(central_point - i, 0) == chrom && matrix_(central_point - i, 2) > chrom_begin) {

                            if (matrix_(central_point - i, 1) <= chrom_begin) {

                                const float partial_peak_left = (float)(matrix_(central_point - i, 2) - chrom_begin)/(mod_end - chrom_begin) * partial_peak_buffer;
                                matrix_(central_point - i, data_type + 3) += partial_peak_left;
                                break;

                            } else {

                                const float partial_peak_left = (float)(matrix_(central_point - i, 2) - matrix_(central_point - i, 1))/(mod_end - chrom_begin) * partial_peak_buffer;
                                matrix_(central_point - i, data_type + 3) += partial_peak_left;
                                partial_peak_buffer -= partial_peak_left;
                                mod_end = matrix_(central_point - i, 1);
                                ++i;
                            }
                        }


                        if (is_log) {

                            partial_peak_buffer = (float)(chrom_end - matrix_(central_point, 2))/(chrom_end - chrom_begin) * pow(2, peak);

                        } else {

                            partial_peak_buffer = (float)(chrom_end - matrix_(central_point, 2))/(chrom_end - chrom_begin) * peak;

                        }

                        i = 1;
                        int mod_start = matrix_(central_point, 2);

                        while (matrix_(central_point + i, 0) == chrom && matrix_(central_point + i, 1) < chrom_end) {

                            if (matrix_(central_point + i, 2) >= chrom_end) {

                                const float partial_peak_right = (float)(chrom_end - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
                                matrix_(central_point + i, data_type + 3) += partial_peak_right;
                                break;

                            } else {

                                const float partial_peak_right = (float)(matrix_(central_point + i, 2) - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
                                matrix_(central_point + i, data_type + 3) += partial_peak_right;
                                partial_peak_buffer -= partial_peak_right;
                                mod_start = matrix_(central_point + i, 2);
                                ++i;
                            }
                        }


                    // sample:  |----|
                    // queue:     |---|
                    //  or
                    // sample:  |---|
                    // queue:     |-|
                    //  or
                    // sample:  |---|
                    // queue:       |--|
                    } else {

                        float partial_peak_buffer;
                        float partial_peak;

                        if (is_log) {

                            // peak value that overlaps with the queue
                            partial_peak = (float)(chrom_end - matrix_(central_point, 1))/(chrom_end - chrom_begin) * pow(2, peak);
                            matrix_(central_point, data_type + 3) += partial_peak;
                            partial_peak_buffer = pow(2, peak) - partial_peak;

                        } else {

                            // peak value that overlaps with the queue
                            partial_peak = (float)(chrom_end - matrix_(central_point, 1))/(chrom_end - chrom_begin) * peak;
                            matrix_(central_point, data_type + 3) += partial_peak;
                            partial_peak_buffer = peak - partial_peak;

                        }

                        int i = 1;
                        int mod_end = matrix_(central_point, 1);

                        while (matrix_(central_point - i, 0) == chrom && matrix_(central_point - i, 2) > chrom_begin) {

                            if (matrix_(central_point - i, 1) <= chrom_begin) {

                                partial_peak = (float)(matrix_(central_point - i, 2) - chrom_begin)/(mod_end - chrom_begin) * partial_peak_buffer;
                                matrix_(central_point - i, data_type + 3) += partial_peak;
                                break;

                            } else {

                                partial_peak = (float)(matrix_(central_point - i, 2) - matrix_(central_point + i, 1))/(mod_end - chrom_begin) * partial_peak_buffer;
                                matrix_(central_point - i, data_type + 3) += partial_peak;
                                partial_peak_buffer -= partial_peak;
                                mod_end = matrix_(central_point - i, 1);
                                ++i;
                            }
                        }
                    }

                // sample: |---|
                // queue:        |---|
                } else {

                    binary_search(chrom, chrom_begin, chrom_end, start_point, central_point - 1, data_type, peak, is_log);
                }

            }

        // matrix_(central_point, 1) < chrom_begin
        } else {

            // sample:        |---|
            // queue:   |---|
            if (matrix_(central_point, 2) <= chrom_begin) {

                binary_search(chrom, chrom_begin, chrom_end, central_point + 1, end_point, data_type, peak, is_log);

            // matrix_(central_point, 2) > chrom_begin
            } else {

                // sample:   |--|
                // queue:  |----|
                //  or
                // sample:   |--|
                // queue:  |------|
                if (matrix_(central_point, 2) >= chrom_end) {

                    if (is_log) {

                        matrix_(central_point, data_type + 3) += pow(2, peak);

                    } else {

                        matrix_(central_point, data_type + 3) += peak;

                    }

                // sample:   |---|
                // queue:  |---|
                } else {

                    float partial_peak;
                    float partial_peak_buffer;

                    if (is_log) {

                        // peak value that overlaps with the queue
                        partial_peak = (float)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * pow(2, peak);
                        partial_peak_buffer = pow(2, peak) - partial_peak;
                        matrix_(central_point, data_type + 3) += partial_peak;

                    } else {

                        // peak value that overlaps with the queue
                        partial_peak = (float)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * peak;
                        partial_peak_buffer = peak - partial_peak;
                        matrix_(central_point, data_type + 3) += partial_peak;

                    }

                    int i = 1;
                    int mod_start = matrix_(central_point, 2);

                    while (matrix_(central_point + i, 0) == chrom && matrix_(central_point + i, 1) < chrom_end) {

                        if (matrix_(central_point + i, 2) >= chrom_end) {

                            partial_peak = (float)(chrom_end - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
                            matrix_(central_point + i, data_type + 3) += partial_peak;
                            break;

                        } else {

                            partial_peak = (float)(matrix_(central_point + i, 2) - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
                            matrix_(central_point + i, data_type + 3) += partial_peak;
                            partial_peak_buffer -= partial_peak;
                            mod_start = matrix_(central_point + i, 2);
                            ++i;
                        }
                    }
                }
            }
        }

    } else {

        // sample:        |---|
        // queue:   |---|
        if (matrix_(central_point, 0) < chrom) {

            binary_search(chrom, chrom_begin, chrom_end, central_point + 1, end_point, data_type, peak, is_log);

        // sample:  |---|
        // queue:         |---|
        } else {

            binary_search(chrom, chrom_begin, chrom_end, start_point, central_point - 1, data_type, peak, is_log);
        }
    }
}
