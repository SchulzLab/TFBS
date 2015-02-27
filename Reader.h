#ifndef READER_H
#define READER_H

#include <string>
#include <stdio.h> // fscanf
#include <utility> // pair
#include <unordered_map>

#include "svm.h"   // struct svm_node
#include "Matrix.h"

// general file reader class
class Reader
{

    public:

        // -- Constructors --

        Reader();

        // @param: number of different data types (e.g. Histone modification measurements
        // and DNAse I data would mean 2 different data types)
        //      i.e. number of additional columns in the matrix
        explicit Reader(int number_of_data_types);

        // move constructor
        Reader(Reader&& other_reader);

        // ---


        // move assignment operator
        Reader& operator=(Reader&& other_reader);

        // reads all files contained in the given directory
        //
        // @param: path to the directory containing files to read
        // void read_files_in_directory(const string& directory_path);

        // reads a single file given by the specified path
        //
        // @param: path to a file to read
        //         corresponding column in matrix_ for the type of data in the file
        //         boolean if peaks in file are log ratio or not
        void read_file(const string& file_path, int data_type, bool is_log);

        // reads a peak file of ENCODE broadpeak format
        //
        // @param: path to a peak file
        void read_peak_file(const string& peak_file_path);

        // get matrix - i.e. get previously read data
        Matrix<float>& get_prev_read_data();
        void print_prev_read_data(ostream& os);

        // get data in libSVM format
        //
        // @param: (will later hold the result)
        //          number_of_data_points -> _number of data points hence size of the following vector
        //          data_points -> vector of data points, each holding all features. see libSVM for more info
        void get_as_libSVM_data(int* number_of_data_points, struct svm_node** data_points);


    private:


        // search in matrix for given genome region, compute overlap and save peak
        //
        // called by   read_file
        //
        // @param: chromosome region start and end point
        //         starting point and end point for binary search
        //         number of data_type (in wiggle files matrix columns = data type + 2)
        //         peak of the chromosome region
        //         boolean if peak should be treated as log ratio
        void binary_search(const int chrom, const int chrom_begin, const int chromend, const int starting_point, const int end_point, const int data_type, const float peak, bool is_log);



        // since we have a strongly modified version of binary search we need
        // to remember if we've just did a right or left jump for the binary search
        // to avoid endless recursive descent
        // bool gone_right = false;
        // bool gone_left = false;


        // matrix that holds data that have been read
        // for wiggle files this should be:
        //      for each genome region a different line containing
        //          genome_region_start   genome_region_end   peak_test_1  peak_test_2 ...
        Matrix<float> matrix_;

        // number of different data types
        // i.e. number of columns in matrix_
        int number_of_data_types_;

        // counts the number of lines given by a peak file
        // required for efficient allocation of memory for new features
        long line_counter_;

        // maps for accessing intern numerical values for chromsomes in matrix
        unordered_map<string, int> map_str_to_chr_;
        unordered_map<int, string> map_chr_to_str_;

        // internal counter for numerical values representing chromosomes in matrix
        int chrom_numerical_;
};

#endif /* READER_H */
