// standard lib
#include <iostream>
#include <fstream>

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
    ,   has_peak_file_(false)
    ,   last_found_pos_(-1)
    ,   actual_file_(NULL)
{
}





Reader::Reader(int number_of_data_types) :
        matrix_()
    ,   number_of_data_types_(number_of_data_types)
    ,   has_peak_file_(false)
    ,   last_found_pos_(-1)
    ,   actual_file_(NULL)
{
}





Reader::Reader(Reader&& other_reader) :
        matrix_(move(other_reader.matrix_))
    ,   number_of_data_types_(other_reader.number_of_data_types_)
    ,   last_found_pos_(other_reader.last_found_pos_)
    ,   last_found_it_(other_reader.last_found_it_)
    ,   actual_file_(other_reader.actual_file_)
{
    other_reader.matrix_ = Matrix<float>();
}





Reader& Reader::operator=(Reader&& other_reader){

    matrix_ = move(other_reader.matrix_);
    number_of_data_types_ = other_reader.number_of_data_types_;
    last_found_pos_ = other_reader.last_found_pos_;
    last_found_it_ = other_reader.last_found_it_;
    actual_file_ = other_reader.actual_file_;
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





void Reader::read_file(const string& file_path, int data_type) {

    actual_file_ = fopen (file_path.c_str(), "r");

    // set last found pos for initial binary search
    last_found_pos_ = matrix_.get_number_of_lines() / 2;

    // set last found it for initial binary search
    last_found_it_ = matrix_.first_line();
    for (int i = 0; i < last_found_pos_; ++i) {

        ++last_found_it_;
    }

    // variables holding the content of one line of the file temporarily
    int chrom_begin, chrom_end, chrom;
    float peak;

    // read one line in the file per loop
    while (fscanf(actual_file_, "%d %d %d %f", &chrom, &chrom_begin, &chrom_end, &peak) != EOF) {

        // store just read data in matrix
        store_data(chrom, chrom_begin, chrom_end, peak, data_type);

    }

    fclose(actual_file_);
}





void Reader::read_peak_file(const string& file_path) {

    has_peak_file_ = true;
    FILE* peak_file = fopen (file_path.c_str(), "r");

    int chrom, chrom_begin, chrom_end;

    // read one line in the file per loop
    while (fscanf(peak_file, "%d %d %d", &chrom, &chrom_begin, &chrom_end) != EOF) {

        vector<float> new_data(number_of_data_types_ + 3, 0);
        new_data[0] = chrom;
        new_data[1] = chrom_begin;
        new_data[2] = chrom_end;
        matrix_.append_new_line(new_data);

    }

    fclose(peak_file);
}





Matrix<float>& Reader::get_prev_read_data() {

    return matrix_;

}





void Reader::store_data(const int chrom, const int chrom_begin, const int chrom_end, const float peak, const int data_type) {

    // checks if chromosome region has been read before
    //   if not -> insert new matrix line holding the data
    if (matrix_.get_number_of_lines() == 0) {

        vector<float> new_data(number_of_data_types_ + 3, 0);
        new_data[0] = chrom;
        new_data[1] = chrom_begin;
        new_data[2] = chrom_end;
        new_data[data_type + 3] = peak;
        matrix_.append_new_line(new_data);

    } else if (chrom_begin > (*(matrix_.last_line()))[2]) {

        if (!has_peak_file_) {

            vector<float> new_data(number_of_data_types_ + 3, 0);
            new_data[0] = chrom;
            new_data[1] = chrom_begin;
            new_data[2] = chrom_end;
            new_data[data_type + 3] = peak;
            matrix_.append_new_line(new_data);
        }

    // else do ovelarp computation
    } else {

        binary_search(chrom, chrom_begin, chrom_end, last_found_pos_, data_type, peak);
    }
}





void Reader::binary_search(const int chrom, const int chrom_begin, const int chrom_end, const int starting_point, const int data_type, const float peak) {


    cout << "new line:" << endl
        << matrix_ << "\n\n\n";
    cout << starting_point << endl;

    cout << "params:   " << chrom_begin << "/" << chrom_end << "/" << peak << "\n\n";

    // if starting point is below zero -> prepend new line
    // or
    // if starting point is greater then the actual number of lines -> append new line
    if (starting_point < 0 || starting_point >= matrix_.get_number_of_lines()) {

        vector<float> new_line (number_of_data_types_ + 3, 0);
        new_line[0] = chrom;
        new_line[1] = chrom_begin;
        new_line[2] = chrom_end;
        new_line[data_type + 3] = peak;

        last_found_it_ = matrix_.insert_new_line(starting_point, new_line);
        last_found_pos_ = 0;
        // reset binary search rememberer
        gone_left = false;
        gone_right = false;
        return;
    }

    // if insertion is needed use last found line as origin for iterator finding
    // for better perfomance
    auto start_it = last_found_it_;

    // find iterator for matrix[starting_point]
    if (starting_point >= last_found_pos_) {

        for (int i = 0; starting_point > last_found_pos_ + i; ++i) {

            ++start_it;
        }
    } else {

        for (int i = 0; starting_point < last_found_pos_ - i; ++i) {

            --start_it;
        }
    }


    if (matrix_(start_it, 0) == chrom) {

        // Info: the +1 offset in partial_peaks is due to the fact that begin as well as end
        //       counts as bin for the chromosome region
        //
        // region overlap tests
        if (matrix_(start_it, 1) >= chrom_begin) {

            if (matrix_(start_it, 1) == chrom_begin) {

                // sample:  |---|
                // queue:   |---|
                //  or
                // sample:  |--|
                // queue:   |---|
                if (matrix_(start_it, 2) >= chrom_end) {

                    matrix_(start_it, data_type + 3) += peak;
                    last_found_pos_ = starting_point;
                    last_found_it_ = start_it;

                // sample: |----|
                // queue:  |---|
                } else {

                    // peak value that overlaps with the queue
                    const float partial_peak = (float)(matrix_(start_it, 1) - chrom_begin + 1)/(chrom_end - chrom_begin + 1) * peak;
                    matrix_(start_it, data_type + 3) += partial_peak;
                    last_found_pos_ = starting_point;
                    last_found_it_ = start_it;
                    // start binary search for the rest of the sample
                    binary_search(chrom, matrix_(start_it, 2) + 1, chrom_end, starting_point + 1, data_type, peak - partial_peak);
                }

            // matrix(start_it, 1) > chrom_begin
            } else {

                if (matrix_(start_it, 1) <= chrom_end) {

                    // sample:  |---|
                    // queue:    |-|
                    if (matrix_(start_it, 2) < chrom_end) {

                        // peak value that overlaps with the queue
                        const float partial_peak = (float)(matrix_(start_it, 2) - matrix_(start_it, 1) + 1)/(chrom_end - chrom_begin + 1) * peak;
                        matrix_(start_it, data_type + 3) += partial_peak;
                        last_found_pos_ = starting_point;
                        last_found_it_ = start_it;

                        // start binary seach for remaining left and right part of the sample
                        const float partial_peak_left = (float)(matrix_(start_it, 1) - chrom_begin)/(chrom_end - chrom_begin + 1) * peak;

                        binary_search(chrom, chrom_begin, matrix_(start_it, 1) - 1, starting_point - 1, data_type, partial_peak_left);

                        const float partial_peak_right = (float)(chrom_end - matrix_(start_it, 2))/(chrom_end - chrom_begin + 1) * peak;
                        binary_search(chrom, matrix_(start_it, 2) + 1, chrom_end, starting_point + 1, data_type, partial_peak_right);

                    // sample:  |----|
                    // queue:     |---|
                    //  or
                    // sample:  |---|
                    // queue:     |-|
                    //  or
                    // sample:  |---|
                    // queue:       |--|
                    } else {

                        // peak value that overlaps with the queue
                        const float partial_peak = (float)(chrom_end - matrix_(start_it, 1) + 1)/(chrom_end - chrom_begin + 1) * peak;
                        matrix_(start_it, data_type + 3) += partial_peak;
                        last_found_pos_ = starting_point;
                        last_found_it_ = start_it;
                        // binary search for remaining peak of the sample
                        binary_search(chrom, chrom_begin, matrix_(start_it, 1) - 1, starting_point - 1, data_type, peak - partial_peak);

                    }

                // sample: |---|
                // queue:        |---|
                } else {

                    // real binary search
                    const int starting_point_shift = (starting_point + 1) / 2;
                    // if binary search not finished
                    if (starting_point_shift > 0 && !gone_right) {

                        gone_left = true;
                        binary_search(chrom, chrom_begin, chrom_end, starting_point - starting_point_shift, data_type, peak);

                    // else insert new element at actual position
                    } else {

                        if (!has_peak_file_) {

                            vector<float> new_line (number_of_data_types_ + 3, 0);
                            new_line[0] = chrom;
                            new_line[1] = chrom_begin;
                            new_line[2] = chrom_end;
                            new_line[data_type + 3] = peak;

                            last_found_it_ = matrix_.insert_new_line(start_it, new_line);
                            last_found_pos_ = starting_point;
                        }
                    }
                }

            }

        // matrix_(start_it, 1) < chrom_begin
        } else {

            if (matrix_(start_it, 2) <= chrom_begin) {

                // sample:        |---|
                // queue:   |---|
                if (matrix_(start_it, 2) < chrom_begin) {

                    const int starting_point_shift = (matrix_.get_number_of_lines() - starting_point) / 2;
                    // if binary search not finished
                    if (starting_point_shift > 0 && !gone_left) {

                        gone_right = true;
                        binary_search(chrom, chrom_begin, chrom_end, starting_point + starting_point_shift, data_type, peak);

                    // else insert new element at actual position
                    } else {

                        if (!has_peak_file_) {

                            vector<float> new_line (number_of_data_types_ + 3, 0);
                            new_line[0] = chrom;
                            new_line[1] = chrom_begin;
                            new_line[2] = chrom_end;
                            new_line[data_type + 3] = peak;

                            last_found_it_ = matrix_.insert_new_line(++start_it, new_line);
                            last_found_pos_ = starting_point + 1;
                        }
                    }

                // sample:      |--|
                // queue:  |----|
                } else {

                    // peak value that overlaps with the queue
                    const float partial_peak = (float)(matrix_(start_it, 2) - chrom_begin + 1)/(chrom_end - chrom_begin + 1) * peak;
                    matrix_(start_it, data_type + 3) += partial_peak;
                    last_found_pos_ = starting_point;
                    last_found_it_ = start_it;
                    // binary search for all bits except the one that overlaps
                    binary_search(chrom, chrom_begin + 1, chrom_end, starting_point + 1, data_type, peak - partial_peak);

                }

            // matrix_(start_it, 2) > chrom_begin
            } else {

                // sample:   |--|
                // queue:  |----|
                //  or
                // sample:   |--|
                // queue:  |------|
                if (matrix_(start_it, 2) >= chrom_end) {

                    matrix_(start_it, data_type + 3) += peak;
                    last_found_pos_ = starting_point;
                    last_found_it_ = start_it;

                // sample:   |---|
                // queue:  |---|
                } else {

                    // peak value that overlaps with the queue
                    const float partial_peak = (float)(matrix_(start_it, 2) - chrom_begin + 1)/(chrom_end - chrom_begin + 1) * peak;
                    matrix_(start_it, data_type + 3) += partial_peak;
                    last_found_pos_ = starting_point;
                    last_found_it_ = start_it;
                    // binary search for all bits except the one that overlaps
                    binary_search(chrom, matrix_(start_it, 2) + 1, chrom_end, starting_point + 1, data_type, peak - partial_peak);

                }
            }
        }

    } else {

        // sample:        |---|
        // queue:   |---|
        if (matrix_(start_it, 0) < chrom) {

            const int starting_point_shift = (matrix_.get_number_of_lines() - starting_point) / 2;
            // if binary search not finished
            if (starting_point_shift > 0 && !gone_left) {

                gone_right = true;
                binary_search(chrom, chrom_begin, chrom_end, starting_point + starting_point_shift, data_type, peak);

            // else insert new element at actual position
            } else {

                if (!has_peak_file_) {

                    vector<float> new_line (number_of_data_types_ + 3, 0);
                    new_line[0] = chrom;
                    new_line[1] = chrom_begin;
                    new_line[2] = chrom_end;
                    new_line[data_type + 3] = peak;

                    last_found_it_ = matrix_.insert_new_line(++start_it, new_line);
                    last_found_pos_ = starting_point + 1;
                }
            }

        // sample:  |---|
        // queue:         |---|
        } else {

            // real binary search
            const int starting_point_shift = (starting_point + 1) / 2;
            // if binary search not finished
            if (starting_point_shift > 0 && !gone_right) {

                gone_left = true;
                binary_search(chrom, chrom_begin, chrom_end, starting_point - starting_point_shift, data_type, peak);

            // else insert new element at actual position
            } else {

                if (!has_peak_file_) {

                    vector<float> new_line (number_of_data_types_ + 3, 0);
                    new_line[0] = chrom;
                    new_line[1] = chrom_begin;
                    new_line[2] = chrom_end;
                    new_line[data_type + 3] = peak;

                    last_found_it_ = matrix_.insert_new_line(start_it, new_line);
                    last_found_pos_ = starting_point;
                }
            }
        }
    }

    // reset binary search rememberer
    gone_left = false;
    gone_right = false;

}
