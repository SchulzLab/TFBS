// standard lib
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>  // abs() exp()
#include <algorithm> // sort()

#ifdef _OPENMP
#include <omp.h>
#endif
// external lib
// #include "boost/filesystem.hpp"

// project
#include "Matrix.h"
#include "Reader.h"

using namespace std;


// adapt this for faster peak file reading
// it should be bigger then the actual number of lines in the file
constexpr int INITIAL_V_CAP = 700000;
constexpr int V_CAP_STEPSIZE = 100000;

// Binning information for regions
constexpr int NUM_BINS = 6;
constexpr int BIN_SIZE = 100;

bool bin_comp (pair<double, int> a, pair<double, int> b) { return a.first < b.first; }




Reader::Reader() :
        matrix_()
    ,   number_of_data_types_(0)
    ,   line_counter_(0)
    ,   chrom_numerical_(0)
{
    init_chr_mapping();
}





Reader::Reader(int number_of_data_types) :
        matrix_()
    ,   number_of_data_types_(number_of_data_types)
    ,   line_counter_(0)
    ,   chrom_numerical_(0)
{
    init_chr_mapping();
}





Reader::Reader(Reader&& other_reader) :
        matrix_(move(other_reader.matrix_))
    ,   number_of_data_types_(other_reader.number_of_data_types_)
    ,   map_str_to_chr_(other_reader.map_str_to_chr_)
    ,   map_chr_to_str_(other_reader.map_chr_to_str_)
    ,   chrom_numerical_(other_reader.chrom_numerical_)
{
    other_reader.matrix_ = Matrix<double>();
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





void Reader::init_matrix(int num) {

    for (int i = 0; i < num; ++i) {

        vector<double> new_feature(line_counter_, .0);
        matrix_.append_new_column(new_feature);
    }
}





void Reader::read_broadpeak_file(const string& file_path) {


    // actual capabilty of one column
    int actually_reserved = INITIAL_V_CAP;

    FILE* peak_file = fopen (file_path.c_str(), "r");

    char chrom[50];
    int chromosome, chrom_begin, chrom_end;

    // initialize and reserve memory for peak file values
    vector<double> chromosome_v;
    chromosome_v.reserve(INITIAL_V_CAP);
    vector<double> chrom_begin_v;
    chrom_begin_v.reserve(INITIAL_V_CAP);
    vector<double> chrom_end_v;
    chrom_end_v.reserve(INITIAL_V_CAP);

    matrix_.append_new_column(chromosome_v);
    matrix_.append_new_column(chrom_begin_v);
    matrix_.append_new_column(chrom_end_v);

    // TODO test  buffer for overlapping regions
    int last_chrom_reg = -(NUM_BINS * BIN_SIZE);

    // read one line in the file per loop
    // skip the additional information of broadpeak format
    while (fscanf(peak_file, "%s %d %d %*s %*s %*s %*s %*s %*s", chrom, &chrom_begin, &chrom_end) == 3) {


        line_counter_ += NUM_BINS;
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

        const int region_center = chrom_begin + (chrom_end - chrom_begin) / 2;
        // TODO overlap tests
        if (region_center - ((NUM_BINS/2)*BIN_SIZE) < last_chrom_reg) {

            cerr << "Overlap of region " << line_counter_/NUM_BINS << " with " << line_counter_/NUM_BINS - 1 << endl;
        }

        for (int bin = -(NUM_BINS/2); bin < (NUM_BINS/2); ++bin) {

            matrix_.get_column(0).push_back(chromosome);
            matrix_.get_column(1).push_back(region_center + bin * BIN_SIZE);
            matrix_.get_column(2).push_back(region_center + (bin + 1) * BIN_SIZE);
        }
        // TODO act. overlapping test buffer
        last_chrom_reg = region_center + (NUM_BINS/2 * BIN_SIZE);

    }

    if (!feof(peak_file)) {

        fprintf(stderr, ("\nA reading error occured while reading peak file \"" + file_path  + "\"\n").c_str());
    }

    fclose(peak_file);
}





void Reader::read_simplebed_file(const string& file_path) {

    // actual capabilty of one column
    int actually_reserved = INITIAL_V_CAP;

    FILE* peak_file = fopen (file_path.c_str(), "r");

    char chrom[50];
    int chromosome, chrom_begin, chrom_end;

    // initialize and reserve memory for peak file values
    vector<double> chromosome_v;
    chromosome_v.reserve(INITIAL_V_CAP);
    vector<double> chrom_begin_v;
    chrom_begin_v.reserve(INITIAL_V_CAP);
    vector<double> chrom_end_v;
    chrom_end_v.reserve(INITIAL_V_CAP);

    matrix_.append_new_column(chromosome_v);
    matrix_.append_new_column(chrom_begin_v);
    matrix_.append_new_column(chrom_end_v);

    // TODO test  buffer for overlapping regions
    int last_chrom_reg = -(NUM_BINS * BIN_SIZE);
    int last_chrom = -1;

    // read one line in the file per loop
    // skip the additional information of broadpeak format
    while (fscanf(peak_file, "%s %d %d", chrom, &chrom_begin, &chrom_end) == 3) {


        line_counter_ += NUM_BINS;
        if (map_str_to_chr_.find(chrom) != map_str_to_chr_.end()) {

            chromosome = map_str_to_chr_[chrom];

        } else {

            map_str_to_chr_[chrom] = chrom_numerical_;
            map_chr_to_str_[chrom_numerical_] = chrom;
            chromosome = chrom_numerical_++;
        }

        const int region_center = chrom_begin + (chrom_end - chrom_begin)/2;
        // TODO overlap tests
        if (chromosome == last_chrom && region_center - ((NUM_BINS/2)*BIN_SIZE) < last_chrom_reg) {

            cerr << "Overlap of region " << line_counter_/NUM_BINS << " with " << line_counter_/NUM_BINS - 1 << endl;
            // skip this region
            line_counter_ -= NUM_BINS;
            continue;
        }

        // if necessary allocate more memory for matrix columns
        if (line_counter_ > actually_reserved) {

            matrix_.get_column(0).reserve(actually_reserved + V_CAP_STEPSIZE);
            matrix_.get_column(1).reserve(actually_reserved + V_CAP_STEPSIZE);
            matrix_.get_column(2).reserve(actually_reserved + V_CAP_STEPSIZE);
            actually_reserved += V_CAP_STEPSIZE;
        }

        for (int bin = -(NUM_BINS/2); bin < (NUM_BINS/2); ++bin) {

            matrix_.get_column(0).push_back(chromosome);
            matrix_.get_column(1).push_back(region_center + bin * BIN_SIZE);
            matrix_.get_column(2).push_back(region_center + (bin + 1) * BIN_SIZE);
        }
        // TODO act. overlapping test buffer
        last_chrom_reg = region_center + (NUM_BINS/2 * BIN_SIZE);
        last_chrom = chromosome;

    }

    if (!feof(peak_file)) {

        fprintf(stderr, ("\nA reading error occured while reading peak file \"" + file_path  + "\"\n").c_str());
    }

    fclose(peak_file);

}





void Reader::read_matrix_file(const string& matrix_file_path) {

    fstream fs (matrix_file_path, fstream::in);

    // extract first line of the file
    string first_line;
    getline(fs, first_line);
    istringstream extractor_first_line(first_line);

    string column_desc;
    int number_of_columns = 0;
    while (extractor_first_line >> column_desc) {

        ++number_of_columns;
    }

    fs.close();

    // Capacity of vectors in Matrix class - reservation due to performance reasons
    constexpr int INITIAL_V_CAP = 100000;
    // constexpr int V_CAP_STEPSIZE = 10000;

    // actual capacity of one column
    // int actually_reserved = INITIAL_V_CAP;

    FILE* matrix_file = fopen (matrix_file_path.c_str(), "r");
    // TODO skip first line

    // initialize and reserve memory for information
    for (int i = 0; i < number_of_columns; ++i) {

        vector<double> chromosome_v;
        chromosome_v.reserve(INITIAL_V_CAP);

        matrix_.append_new_column(chromosome_v);

    }

    // TODO read one single entry per loop and write it into matrix[entry_number/number_of_columns][entry_number%number_of_columns]
    //
    // read one line in the file per outer loop
    //
    // number of entrys totally read
    int entry_number = 0;

    float data;
    // we can assume that even the name of the chromosome is a number since we require
    // this as a format convention for our matrix and our produced output will always
    // match chrNUMBER to NUMBER when printing the just read data
    while (fscanf(matrix_file, "%f", &data) == 1) {


        // TODO: solve outprint issue with missing map entrys
        //          i.e. string reading issue every number_of_columns loopcycle
        //
        // if (map_str_to_chr_.find(chrom) != map_str_to_chr_.end()) {
        //
        //     chromosome = map_str_to_chr_[chrom];
        //
        // } else {
        //
        //     map_str_to_chr_[chrom] = chrom_numerical_;
        //     map_chr_to_str_[chrom_numerical_] = chrom;
        //     chromosome = chrom_numerical_++;
        // }
        //
        // // if necessary allocate more memory for matrix columns
        // if (line_counter_ > actually_reserved) {
        //
        //     matrix_.get_column(0).reserve(actually_reserved + V_CAP_STEPSIZE);
        //     matrix_.get_column(1).reserve(actually_reserved + V_CAP_STEPSIZE);
        //     matrix_.get_column(2).reserve(actually_reserved + V_CAP_STEPSIZE);
        //     actually_reserved += V_CAP_STEPSIZE;
        // }
        //
        // matrix_.get_column(0).push_back(chromosome);
        // matrix_.get_column(1).push_back(chrom_begin);
        // matrix_.get_column(2).push_back(chrom_end);

        ++entry_number;
    }

    line_counter_ = entry_number / number_of_columns;

    if (!feof(matrix_file) || entry_number % number_of_columns) {

        fprintf(stderr, ("\nA reading error occured while reading peak file \"" + matrix_file_path  + "\"\n").c_str());
    }

    fclose(matrix_file);
}





void Reader::print_prev_read_data(ostream& os) {


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

            os << setw(10) << right << fixed << matrix_(lines, columns) << "\t";
        }
        os << "\n";
    }
}





void Reader::collateBinnedRegions() {

    if (NUM_BINS != 1) {

        Matrix<double> collated_matrix (line_counter_ / NUM_BINS, number_of_data_types_ * NUM_BINS + 3, 0);

        // iterate over all binned regions
        for (int region_num = 0; region_num < line_counter_/ NUM_BINS; ++region_num) {

            collated_matrix(region_num, 0) = matrix_(region_num * NUM_BINS, 0);
            collated_matrix(region_num, 1) = matrix_(region_num * NUM_BINS, 1);
            collated_matrix(region_num, 2) = matrix_(region_num * NUM_BINS + NUM_BINS - 1, 2);

            // convert row entries for all bins of this region to column entries
            // make this for all features
            for (int feature = 0; feature < number_of_data_types_; ++feature) {

                for (int bin = 0; bin < NUM_BINS; ++bin) {

                    collated_matrix(region_num, (feature * NUM_BINS) + bin + 3) = matrix_(region_num * NUM_BINS + bin, feature + 3);
                }
            }
        }

        // // collate all binned regions
        // for (int region_num = 0; region_num < line_counter_ / NUM_BINS; ++region_num) {
        //
        //     // init new matrix with top values
        //     collated_matrix(region_num, 0) = matrix_(region_num * NUM_BINS, 0);
        //     collated_matrix(region_num, 1) = matrix_(region_num * NUM_BINS, 1);
        //     collated_matrix(region_num, 2) = matrix_(region_num * NUM_BINS + NUM_BINS - 1, 2);
        //     // for each feature search the best bins
        //     for (int feature = 3; feature < matrix_.get_number_of_columns(); ++feature) {
        //
        //         // get all bins of one region
        //         vector<pair<double, int>> bins(NUM_BINS);
        //         for (int bin_num = 0; bin_num < NUM_BINS; ++bin_num) {
        //
        //             bins[bin_num] = pair<double, int>(matrix_(region_num * NUM_BINS + bin_num, feature), bin_num);
        //         }
        //
        //         // sort into ascending order
        //         sort(bins.begin(), bins.end(), bin_comp);
        //
        //         // pick best bins
        //         auto region_it = --bins.end();
        //         for (int top_bin = 0; top_bin < NUM_CHOSEN_BINS; ++top_bin, --region_it) {
        //
        //             const int bin_region = region_num * NUM_BINS + region_it->second;
        //             // collate and normalize on length
        //             // collated_matrix(region_num, feature) += matrix_(bin_region, feature) / (matrix_(bin_region, 2) - matrix_(bin_region, 1));
        //             collated_matrix(region_num, feature) += matrix_(bin_region, feature);
        //         }
        //     }
        // }
        matrix_ = move(collated_matrix);
        line_counter_ /= NUM_BINS;

    }
}





Matrix<double>& Reader::get_prev_read_data() {

    return matrix_;
}





void Reader::rescale_data() {

    // search maxima and minima
    vector<double> max_val(number_of_data_types_);
    vector<double> min_val(number_of_data_types_);
    for (int sample_ind = 0; sample_ind < matrix_.get_number_of_lines(); ++sample_ind) {

        for (int feature_ind = 0; feature_ind < matrix_.get_number_of_columns() - 3; ++feature_ind) {

            max_val[feature_ind] < matrix_(sample_ind, feature_ind + 3) ? max_val[feature_ind] = matrix_(sample_ind, feature_ind + 3) : 0;
            min_val[feature_ind] > matrix_(sample_ind, feature_ind + 3) ? min_val[feature_ind] = matrix_(sample_ind, feature_ind + 3) : 0;
        }
    }

    // rescale every feature
    for (int sample_ind = 0; sample_ind < matrix_.get_number_of_lines(); ++sample_ind) {

        for (int feature_ind = 0; feature_ind < matrix_.get_number_of_columns() - 3; ++feature_ind) {

            matrix_(sample_ind, feature_ind + 3) = (matrix_(sample_ind, feature_ind + 3) - min_val[feature_ind])/(max_val[feature_ind] - min_val[feature_ind]);
        }
    }
}





// TODO TODO TODO fix everywhere
//
// TODO: implement raw read count based on peak and remove partial peak -> safer and better performance
// void Reader::binary_search(const int chrom, const int chrom_begin, const int chrom_end, const int start_point, const int end_point, const int data_type, const double peak, bool is_log) {
//
//
//     // termination
//     if (start_point > end_point) {
//
//         return;
//     }
//
//     // center of one vector
//     const int central_point = start_point + ((end_point - start_point) / 2);
//
//     // tests if query is on the same chromosome as new sample
//     if (matrix_(central_point, 0) == chrom) {
//
//         // region overlap tests between query and sample
//         if (matrix_(central_point, 1) >= chrom_begin) {
//
//             if (matrix_(central_point, 1) == chrom_begin) {
//
//                 // sample:  |---|
//                 // query:   |---|
//                 //  or
//                 // sample:  |--|
//                 // query:   |---|
//                 if (matrix_(central_point, 2) >= chrom_end) {
//
//                     if (is_log) {
//
//                         matrix_(central_point, data_type + 3) += pow(2, peak);
//
//                     } else {
//
//                         matrix_(central_point, data_type + 3) += peak;
//                     }
//
//                 // sample: |----|
//                 // query:  |---|
//                 } else {
//
//                     // remaining peak value
//                     double partial_peak_buffer;
//                     double partial_peak;
//
//                     if (is_log) {
//
//                         // peak value that overlaps with the query
//                         partial_peak = (double)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * pow(2, peak);
//                         matrix_(central_point, data_type + 3) += partial_peak;
//                         partial_peak_buffer = pow(2, peak) - partial_peak;
//
//                     } else {
//
//                         partial_peak = (double)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * peak;
//                         matrix_(central_point, data_type + 3) += partial_peak;
//                         partial_peak_buffer = peak - partial_peak;
//
//                     }
//
//
//                     // increase position in each step by one
//                     int i = 1;
//                     int mod_start = matrix_(central_point, 2);
//
//                     while (central_point + i <= line_counter_ - 1 && matrix_(central_point + i, 0) == chrom && matrix_(central_point + i, 1) < chrom_end) {
//
//                         if (matrix_(central_point + i, 2) >= chrom_end) {
//
//                             partial_peak = (double)(chrom_end - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
//                             matrix_(central_point, data_type + 3) += partial_peak;
//                             break;
//
//                         } else {
//
//                             partial_peak = (double)(matrix_(central_point + i, 2) - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
//                             matrix_(central_point + i, data_type + 3) += partial_peak;
//                             partial_peak_buffer -= partial_peak;
//                             mod_start = matrix_(central_point + i, 2);
//                             ++i;
//                         }
//                     }
//                 }
//
//             // matrix(central_point, 1) > chrom_begin
//             } else {
//
//                 if (matrix_(central_point, 1) < chrom_end) {
//
//                     // sample:  |---|
//                     // query:    |-|
//                     if (matrix_(central_point, 2) < chrom_end) {
//
//                         // remaining peak value
//                         double partial_peak_buffer;
//
//                         if (is_log) {
//
//                             // peak value that overlaps with the query
//                             const double partial_peak = (double)(matrix_(central_point, 2) - matrix_(central_point, 1))/(chrom_end - chrom_begin) * pow(2, peak);
//                             matrix_(central_point, data_type + 3) += partial_peak;
//
//                             partial_peak_buffer = (double)(matrix_(central_point, 1) - chrom_begin)/(chrom_end - chrom_begin) * pow(2, peak);
//
//                         } else {
//
//                             // peak value that overlaps with the query
//                             const double partial_peak = (double)(matrix_(central_point, 2) - matrix_(central_point, 1))/(chrom_end - chrom_begin) * peak;
//                             matrix_(central_point, data_type + 3) += partial_peak;
//
//                             partial_peak_buffer = (double)(matrix_(central_point, 1) - chrom_begin)/(chrom_end - chrom_begin) * peak;
//                         }
//
//
//                         // overlap of left side
//                         int i = 1;
//                         int mod_end = matrix_(central_point, 1);
//
//                         while (central_point - i >= 0 && matrix_(central_point - i, 0) == chrom && matrix_(central_point - i, 2) > chrom_begin) {
//
//                             if (matrix_(central_point - i, 1) <= chrom_begin) {
//
//                                 const double partial_peak_left = (double)(matrix_(central_point - i, 2) - chrom_begin)/(mod_end - chrom_begin) * partial_peak_buffer;
//                                 matrix_(central_point - i, data_type + 3) += partial_peak_left;
//                                 break;
//
//                             } else {
//
//                                 const double partial_peak_left = (double)(matrix_(central_point - i, 2) - matrix_(central_point - i, 1))/(mod_end - chrom_begin) * partial_peak_buffer;
//                                 matrix_(central_point - i, data_type + 3) += partial_peak_left;
//                                 partial_peak_buffer -= partial_peak_left;
//                                 mod_end = matrix_(central_point - i, 1);
//                                 ++i;
//                             }
//                         }
//
//
//                         if (is_log) {
//
//                             partial_peak_buffer = (double)(chrom_end - matrix_(central_point, 2))/(chrom_end - chrom_begin) * pow(2, peak);
//
//                         } else {
//
//                             partial_peak_buffer = (double)(chrom_end - matrix_(central_point, 2))/(chrom_end - chrom_begin) * peak;
//
//                         }
//
//
//                         // overlap of right side
//                         i = 1;
//                         int mod_start = matrix_(central_point, 2);
//
//                         while (central_point + i <= line_counter_ - 1 && matrix_(central_point + i, 0) == chrom && matrix_(central_point + i, 1) < chrom_end) {
//
//                             if (matrix_(central_point + i, 2) >= chrom_end) {
//
//                                 const double partial_peak_right = (double)(chrom_end - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
//                                 matrix_(central_point + i, data_type + 3) += partial_peak_right;
//                                 break;
//
//                             } else {
//
//                                 const double partial_peak_right = (double)(matrix_(central_point + i, 2) - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
//                                 matrix_(central_point + i, data_type + 3) += partial_peak_right;
//                                 partial_peak_buffer -= partial_peak_right;
//                                 mod_start = matrix_(central_point + i, 2);
//                                 ++i;
//                             }
//                         }
//
//
//                     // sample:  |----|
//                     // query:     |---|
//                     //  or
//                     // sample:  |---|
//                     // query:     |-|
//                     //  or
//                     // sample:  |---|
//                     // query:       |--|
//                     } else {
//
//                         double partial_peak_buffer;
//                         double partial_peak;
//
//                         if (is_log) {
//
//                             // peak value that overlaps with the query
//                             partial_peak = (double)(chrom_end - matrix_(central_point, 1))/(chrom_end - chrom_begin) * pow(2, peak);
//                             matrix_(central_point, data_type + 3) += partial_peak;
//                             partial_peak_buffer = pow(2, peak) - partial_peak;
//
//                         } else {
//
//                             // peak value that overlaps with the query
//                             partial_peak = (double)(chrom_end - matrix_(central_point, 1))/(chrom_end - chrom_begin) * peak;
//                             matrix_(central_point, data_type + 3) += partial_peak;
//                             partial_peak_buffer = peak - partial_peak;
//
//                         }
//
//                         int i = 1;
//                         int mod_end = matrix_(central_point, 1);
//
//                         while (central_point - i >= 0 && matrix_(central_point - i, 0) == chrom && matrix_(central_point - i, 2) > chrom_begin) {
//
//                             if (matrix_(central_point - i, 1) <= chrom_begin) {
//
//                                 partial_peak = (double)(matrix_(central_point - i, 2) - chrom_begin)/(mod_end - chrom_begin) * partial_peak_buffer;
//                                 matrix_(central_point - i, data_type + 3) += partial_peak;
//                                 break;
//
//                             } else {
//
//                                 partial_peak = (double)(matrix_(central_point - i, 2) - matrix_(central_point - i, 1))/(mod_end - chrom_begin) * partial_peak_buffer;
//                                 matrix_(central_point - i, data_type + 3) += partial_peak;
//                                 partial_peak_buffer -= partial_peak;
//                                 mod_end = matrix_(central_point - i, 1);
//                                 ++i;
//                             }
//                         }
//                     }
//
//                 // sample: |---|
//                 // query:        |---|
//                 } else {
//
//                     binary_search(chrom, chrom_begin, chrom_end, start_point, central_point - 1, data_type, peak, is_log);
//                 }
//
//             }
//
//         // matrix_(central_point, 1) < chrom_begin
//         } else {
//
//             // sample:        |---|
//             // query:   |---|
//             if (matrix_(central_point, 2) <= chrom_begin) {
//
//                 binary_search(chrom, chrom_begin, chrom_end, central_point + 1, end_point, data_type, peak, is_log);
//
//             // matrix_(central_point, 2) > chrom_begin
//             } else {
//
//                 // sample:   |--|
//                 // query:  |----|
//                 //  or
//                 // sample:   |--|
//                 // query:  |------|
//                 if (matrix_(central_point, 2) >= chrom_end) {
//
//                     if (is_log) {
//
//                         matrix_(central_point, data_type + 3) += pow(2, peak);
//
//                     } else {
//
//                         matrix_(central_point, data_type + 3) += peak;
//
//                     }
//
//                 // sample:   |---|
//                 // query:  |---|
//                 } else {
//
//                     double partial_peak;
//                     double partial_peak_buffer;
//
//                     if (is_log) {
//
//                         // peak value that overlaps with the query
//                         partial_peak = (double)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * pow(2, peak);
//                         partial_peak_buffer = pow(2, peak) - partial_peak;
//                         matrix_(central_point, data_type + 3) += partial_peak;
//
//                     } else {
//
//                         // peak value that overlaps with the query
//                         partial_peak = (double)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * peak;
//                         partial_peak_buffer = peak - partial_peak;
//                         matrix_(central_point, data_type + 3) += partial_peak;
//
//                     }
//
//                     int i = 1;
//                     int mod_start = matrix_(central_point, 2);
//
//                     while (central_point + i<= line_counter_ - 1 && matrix_(central_point + i, 0) == chrom && matrix_(central_point + i, 1) < chrom_end) {
//
//                         if (matrix_(central_point + i, 2) >= chrom_end) {
//
//                             partial_peak = (double)(chrom_end - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
//                             matrix_(central_point + i, data_type + 3) += partial_peak;
//                             break;
//
//                         } else {
//
//                             partial_peak = (double)(matrix_(central_point + i, 2) - matrix_(central_point + i, 1))/(chrom_end - mod_start) * partial_peak_buffer;
//                             matrix_(central_point + i, data_type + 3) += partial_peak;
//                             partial_peak_buffer -= partial_peak;
//                             mod_start = matrix_(central_point + i, 2);
//                             ++i;
//                         }
//                     }
//                 }
//             }
//         }
//
//     } else {
//
//         // sample:        |---|
//         // query:   |---|
//         if (matrix_(central_point, 0) < chrom) {
//
//             binary_search(chrom, chrom_begin, chrom_end, central_point + 1, end_point, data_type, peak, is_log);
//
//         // sample:  |---|
//         // query:         |---|
//         } else {
//
//             binary_search(chrom, chrom_begin, chrom_end, start_point, central_point - 1, data_type, peak, is_log);
//         }
//     }
// }
void Reader::binary_search(const int chrom, const int chrom_begin, const int chrom_end, const int start_point, const int end_point, const int data_type, const double peak, bool is_log) {


    // termination
    if (start_point > end_point) {

        return;
    }

    // center for binary search jump
    const int central_point = start_point + ((end_point - start_point) / 2);

    // tests if query is on the same chromosome as new sample
    if (matrix_(central_point, 0) == chrom) {

        // region overlap tests between query and sample
        if (matrix_(central_point, 1) >= chrom_begin) {

            if (matrix_(central_point, 1) == chrom_begin) {

                // sample:  |---|
                // query:   |---|
                //  or
                // sample:  |--|
                // query:   |---|
                if (matrix_(central_point, 2) >= chrom_end) {

                    matrix_(central_point, data_type + 3) += (is_log ? pow(2, peak) : peak);

                // sample: |----|
                // query:  |---|
                } else {


                    // peak value that overlaps with the query
                    const double partial_peak = (double)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                    matrix_(central_point, data_type + 3) += partial_peak;


                    // increase position in each step by one
                    int i = 1;

                    while (central_point + i <= line_counter_ - 1 && matrix_(central_point + i, 0) == chrom && matrix_(central_point + i, 1) < chrom_end) {

                        if (matrix_(central_point + i, 2) >= chrom_end) {

                            const double partial_peak2 = (double)(chrom_end - matrix_(central_point + i, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                            matrix_(central_point, data_type + 3) += partial_peak2;
                            break;

                        } else {

                            const double partial_peak2 = (double)(matrix_(central_point + i, 2) - matrix_(central_point + i, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                            matrix_(central_point + i, data_type + 3) += partial_peak2;
                            ++i;
                        }
                    }
                }

            // matrix(central_point, 1) > chrom_begin
            } else {

                if (matrix_(central_point, 1) < chrom_end) {

                    // sample:  |---|
                    // query:    |-|
                    if (matrix_(central_point, 2) < chrom_end) {

                        // peak value that overlaps with the query
                        const double partial_peak = (double)(matrix_(central_point, 2) - matrix_(central_point, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                        matrix_(central_point, data_type + 3) += partial_peak;


                        // overlap of left side
                        int i = 1;

                        while (central_point - i >= 0 && matrix_(central_point - i, 0) == chrom && matrix_(central_point - i, 2) > chrom_begin) {

                            if (matrix_(central_point - i, 1) <= chrom_begin) {

                                const double partial_peak_left = (double)(matrix_(central_point - i, 2) - chrom_begin)/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                                matrix_(central_point - i, data_type + 3) += partial_peak_left;
                                break;

                            } else {

                                const double partial_peak_left = (double)(matrix_(central_point - i, 2) - matrix_(central_point - i, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                                matrix_(central_point - i, data_type + 3) += partial_peak_left;
                                ++i;
                            }
                        }

                        // overlap of right side
                        i = 1;

                        while (central_point + i <= line_counter_ - 1 && matrix_(central_point + i, 0) == chrom && matrix_(central_point + i, 1) < chrom_end) {

                            if (matrix_(central_point + i, 2) >= chrom_end) {

                                const double partial_peak_right = (double)(chrom_end - matrix_(central_point + i, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                                matrix_(central_point + i, data_type + 3) += partial_peak_right;
                                break;

                            } else {

                                const double partial_peak_right = (double)(matrix_(central_point + i, 2) - matrix_(central_point + i, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                                matrix_(central_point + i, data_type + 3) += partial_peak_right;
                                ++i;
                            }
                        }


                    // sample:  |----|
                    // query:     |---|
                    //  or
                    // sample:  |---|
                    // query:     |-|
                    //  or
                    // sample:  |---|
                    // query:       |--|
                    } else {

                        // peak value that overlaps with the query
                        const double partial_peak = (double)(chrom_end - matrix_(central_point, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                        matrix_(central_point, data_type + 3) += partial_peak;

                        int i = 1;

                        while (central_point - i >= 0 && matrix_(central_point - i, 0) == chrom && matrix_(central_point - i, 2) > chrom_begin) {

                            if (matrix_(central_point - i, 1) <= chrom_begin) {

                                const double partial_peak2 = (double)(matrix_(central_point - i, 2) - chrom_begin)/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                                matrix_(central_point - i, data_type + 3) += partial_peak2;
                                break;

                            } else {

                                const double partial_peak2 = (double)(matrix_(central_point - i, 2) - matrix_(central_point - i, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                                matrix_(central_point - i, data_type + 3) += partial_peak2;
                                ++i;
                            }
                        }
                    }

                // sample: |---|
                // query:        |---|
                } else {

                    binary_search(chrom, chrom_begin, chrom_end, start_point, central_point - 1, data_type, peak, is_log);
                }

            }

        // matrix_(central_point, 1) < chrom_begin
        } else {

            // sample:        |---|
            // query:   |---|
            if (matrix_(central_point, 2) <= chrom_begin) {

                binary_search(chrom, chrom_begin, chrom_end, central_point + 1, end_point, data_type, peak, is_log);

            // matrix_(central_point, 2) > chrom_begin
            } else {

                // sample:   |--|
                // query:  |----|
                //  or
                // sample:   |--|
                // query:  |------|
                if (matrix_(central_point, 2) >= chrom_end) {

                    matrix_(central_point, data_type + 3) += (is_log ? pow(2, peak) : peak);


                // sample:   |---|
                // query:  |---|
                } else {

                    // peak value that overlaps with the query
                    const double partial_peak = (double)(matrix_(central_point, 2) - chrom_begin)/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                    matrix_(central_point, data_type + 3) += partial_peak;


                    int i = 1;

                    while (central_point + i<= line_counter_ - 1 && matrix_(central_point + i, 0) == chrom && matrix_(central_point + i, 1) < chrom_end) {

                        if (matrix_(central_point + i, 2) >= chrom_end) {

                            const double partial_peak2 = (double)(chrom_end - matrix_(central_point + i, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                            matrix_(central_point + i, data_type + 3) += partial_peak2;
                            break;

                        } else {

                            const double partial_peak2 = (double)(matrix_(central_point + i, 2) - matrix_(central_point + i, 1))/(chrom_end - chrom_begin) * (is_log ? pow(2, peak) : peak);
                            matrix_(central_point + i, data_type + 3) += partial_peak2;
                            ++i;
                        }
                    }
                }
            }
        }

    } else {

        // sample:        |---|
        // query:   |---|
        if (matrix_(central_point, 0) < chrom) {

            binary_search(chrom, chrom_begin, chrom_end, central_point + 1, end_point, data_type, peak, is_log);

        // sample:  |---|
        // query:         |---|
        } else {

            binary_search(chrom, chrom_begin, chrom_end, start_point, central_point - 1, data_type, peak, is_log);
        }
    }
}





// TODO TODO TODO
void Reader::binary_search_mean(const int chrom, const int chrom_begin, const int chrom_end, const int start_point, const int end_point, const int data_type, const double value) {

    // termination
    if (start_point > end_point) {

        return;
    }
}





void Reader::init_chr_mapping() {

    map_str_to_chr_["chr1"] = chrom_numerical_;
    map_str_to_chr_["1"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "1";

    map_str_to_chr_["chr2"] = chrom_numerical_;
    map_str_to_chr_["2"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "2";

    map_str_to_chr_["chr3"] = chrom_numerical_;
    map_str_to_chr_["3"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "3";

    map_str_to_chr_["chr4"] = chrom_numerical_;
    map_str_to_chr_["4"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "4";

    map_str_to_chr_["chr5"] = chrom_numerical_;
    map_str_to_chr_["5"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "5";

    map_str_to_chr_["chr6"] = chrom_numerical_;
    map_str_to_chr_["6"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "6";

    map_str_to_chr_["chr7"] = chrom_numerical_;
    map_str_to_chr_["7"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "7";

    map_str_to_chr_["chr8"] = chrom_numerical_;
    map_str_to_chr_["8"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "8";

    map_str_to_chr_["chr9"] = chrom_numerical_;
    map_str_to_chr_["9"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "9";

    map_str_to_chr_["chr10"] = chrom_numerical_;
    map_str_to_chr_["10"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "10";

    map_str_to_chr_["chr11"] = chrom_numerical_;
    map_str_to_chr_["11"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "11";

    map_str_to_chr_["chr12"] = chrom_numerical_;
    map_str_to_chr_["12"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "12";

    map_str_to_chr_["chr13"] = chrom_numerical_;
    map_str_to_chr_["13"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "13";

    map_str_to_chr_["chr14"] = chrom_numerical_;
    map_str_to_chr_["14"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "14";

    map_str_to_chr_["chr15"] = chrom_numerical_;
    map_str_to_chr_["15"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "15";

    map_str_to_chr_["chr16"] = chrom_numerical_;
    map_str_to_chr_["16"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "16";

    map_str_to_chr_["chr17"] = chrom_numerical_;
    map_str_to_chr_["17"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "17";

    map_str_to_chr_["chr18"] = chrom_numerical_;
    map_str_to_chr_["18"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "18";

    map_str_to_chr_["chr19"] = chrom_numerical_;
    map_str_to_chr_["19"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "19";

    map_str_to_chr_["chr20"] = chrom_numerical_;
    map_str_to_chr_["20"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "20";

    map_str_to_chr_["chr21"] = chrom_numerical_;
    map_str_to_chr_["21"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "21";

    map_str_to_chr_["chr22"] = chrom_numerical_;
    map_str_to_chr_["22"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "22";

    map_str_to_chr_["chrX"] = chrom_numerical_;
    map_str_to_chr_["X"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "X";

    map_str_to_chr_["chrY"] = chrom_numerical_;
    map_str_to_chr_["Y"] = chrom_numerical_;
    map_chr_to_str_[chrom_numerical_++] = "Y";
}
