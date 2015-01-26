// standard lib
#include <iostream>
#include <fstream>

#include <math.h>  // round

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
    ,   last_found_pos_(-1)
    ,   actual_file_(NULL)
{
}





Reader::Reader(int number_of_data_types) :
        matrix_()
    ,   number_of_data_types_(number_of_data_types)
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
    other_reader.matrix_ = Matrix<int>();
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

    // set last foun pos for initial binary search
    last_found_pos_ = matrix_.get_number_of_lines() / 2;

    // variables holding the content of one line of the file temporarily
    int chrom_begin, chrom_end, peak;

    // read one line in the file
    while (fscanf(actual_file_, "%d %d %d", &chrom_begin, &chrom_end, &peak) != EOF) {

        // store just read data in matrix
        store_data(chrom_begin, chrom_end, peak, data_type);

    }

    fclose(actual_file_);
}





Matrix<int>& Reader::get_prev_read_data() {

    return matrix_;

}





void Reader::store_data(const int chrom_begin, const int chrom_end, const int peak, const int data_type) {

    // checks if chromosome region has been read before
    //   if not -> insert new matrix line holding the data
    if (matrix_.get_number_of_lines() == 0) {

        vector<int> new_data(number_of_data_types_ + 2);
        new_data[0] = chrom_begin;
        new_data[1] = chrom_end;
        new_data[data_type + 2] = peak;
        matrix_.append_new_line(new_data);

    } else if (chrom_begin > matrix_.last_line()[1]) {

        vector<int> new_data(number_of_data_types_ + 2);
        new_data[0] = chrom_begin;
        new_data[1] = chrom_end;
        new_data[data_type + 2] = peak;
        matrix_.append_new_line(new_data);

    // else do ovelarp computation
    } else {

        binary_search(chrom_begin, chrom_end, last_found_pos_, data_type, peak);
    }
}





void Reader::binary_search(const int chrom_begin, const int chrom_end, const int starting_point, const int data_type, const int peak) {

            // TODO: Shouldn't be float more appropriate?

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

    // Info: the +1 offset in partial_peaks is due to the fact that begin as well as end
    //       counts as bin for the chromosome region
    //
    // region overlap tests
    if (matrix_(start_it, 0) >= chrom_begin) {

        if (matrix_(start_it, 0) == chrom_begin) {

            // sample:  |---|
            // queue:   |---|
            //  or
            // sample:  |--|
            // queue:   |---|
            if (matrix_(start_it, 1) <= chrom_end) {

                matrix_(start_it, data_type + 2) += peak;
                last_found_pos_ = starting_point;
                last_found_it_ = start_it;

            // sample: |----|
            // queue:  |---|
            } else {

                // peak value that overlaps with the queue
                const int partial_peak = round(
                        (matrix_(start_it, 1) - chrom_begin + 1)/(chrom_end - chrom_begin)
                        ) * peak;
                matrix_(start_it, data_type + 2) += partial_peak;
                last_found_pos_ = starting_point;
                last_found_it_ = start_it;
                binary_search(chrom_begin, chrom_end, starting_point + 1, data_type, peak);
            }

        // matrix(starting_point, 0) > chrom_begin
        } else {

            if (matrix_(starting_point, 0) <= chrom_end) {

                // sample:  |---|
                // queue:    |-|
                if (matrix_(starting_point, 1) < chrom_end) {

                    // peak value that overlaps with the queue
                    const int partial_peak = round(
                            (matrix_(start_it, 1) - matrix_(start_it, 0) + 1)/(chrom_end - chrom_begin)
                            ) * peak;
                    matrix_(start_it, data_type + 2) += partial_peak;
                    last_found_pos_ = starting_point;
                    last_found_it_ = start_it;
                    binary_search(chrom_begin, chrom_end, starting_point - 1, data_type, peak);
                    binary_search(chrom_begin, chrom_end, starting_point + 1, data_type, peak);

                // sample:  |----|
                // queue:     |--|
                } else {

                    // peak value that overlaps with the queue
                    const int partial_peak = round(
                            (matrix_(start_it, 1) - matrix_(start_it, 0) + 1)/(chrom_end - chrom_begin)
                            ) * peak;
                    matrix_(start_it, data_type + 2) += partial_peak;
                    last_found_pos_ = starting_point;
                    last_found_it_ = start_it;
                    binary_search(chrom_begin, chrom_end, starting_point - 1, data_type, peak);

                }

            // sample: |---|
            // queue:        |---|
            } else {

                // real binary search
                const int starting_point_shift = (starting_point + 1) / 2;
                // if binary search not finished
                if (starting_point_shift > 0) {

                    binary_search(chrom_begin, chrom_end, starting_point - starting_point_shift, data_type, peak);

                // else insert new element at actual position
                } else {

                    vector<int> new_line (number_of_data_types_ + 2, - 1);
                    new_line[0] = chrom_begin;
                    new_line[1] = chrom_end;
                    new_line[data_type] = peak;

                    last_found_it_ = matrix_.insert_new_line(start_it, new_line);
                    last_found_pos_ = starting_point;
                }
            }

        }

    // matrix_(starting_point, 0) < chrom_begin
    } else {

        if (matrix_(start_it, 1) <= chrom_end) {

            // sample:        |---|
            // queue:   |---|
            if (matrix_(start_it, 1) < chrom_end) {

                const int starting_point_shift = (matrix_.get_number_of_lines() - starting_point + 1) / 2;
                // if binary search not finished
                if (starting_point_shift > 0) {

                    binary_search(chrom_begin, chrom_end, starting_point - starting_point_shift, data_type, peak);

                // else insert new element at actual position
                } else {

                    vector<int> new_line (number_of_data_types_ + 2, - 1);
                    new_line[0] = chrom_begin;
                    new_line[1] = chrom_end;
                    new_line[data_type] = peak;

                    last_found_it_ = matrix_.insert_new_line(start_it, new_line);
                    last_found_pos_ = starting_point;
                }

            // sample:      |--|
            // queue:  |----|
            } else {

                // peak value that overlaps with the queue
                const int partial_peak = round(
                        (matrix_(start_it, 1) - chrom_begin + 1)/(chrom_end - chrom_begin)
                        ) * peak;
                matrix_(start_it, data_type + 2) += partial_peak;
                last_found_pos_ = starting_point;
                last_found_it_ = start_it;

            }

        // sample:   |--|
        // queue:  |-----|
        } else {

            matrix_(start_it, data_type + 2) += peak;
            last_found_pos_ = starting_point;
            last_found_it_ = start_it;
        }
    }
}
