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
        //         and DNAse I data would mean 2 different data types)
        //          i.e. number of additional columns in the matrix
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
        void read_broadpeak_file(const string& peak_file_path);

        // reads a peak file of BED format, with only the first three required fields (chromosome, start, end)
        //
        // @param: path to a peak file
        void read_simplebed_file(const string& peak_file_path);

        // reads a matrix out of a given file
        //
        // @param: a file to read from
        void read_matrix_file(const string& matrix_file_path);

        // initialize matrix to avoid unordered modification of wrong columns in parallel region
        //
        // @param: number of files (i.e. number of columns to append to matrix)
        void init_matrix(int num);

        // get matrix - i.e. get previously read data
        Matrix<double>& get_prev_read_data();
        void print_prev_read_data(ostream& os);

        // get data in libSVM format
        //
        // @param: (will later hold the result)
        //          number_of_data_points - _number of data points, i.e. size of the following vector
        //          data_points - vector of data points, each holding all features. see libSVM for more info
        void get_as_libSVM_data(int* number_of_data_points, struct svm_node** data_points);


        // rescale data such that all features have values between 0 and 1
        // min - max rescaling
        void rescale_data();

        // collates the binned regions, adds up only some of the bins to get value for former region
        void collateBinnedRegions();


    private:


        // search in matrix for given genome region, compute overlap and save
        // binary_search computes overlap for raw read counts (or log transformed raw counts)
        // binary_search_mean computes overlap for values that reflect a mean of a region
        //
        // called by   read_file
        //
        // @param: chromosome region start and end point
        //         starting point and end point for binary search
        //         number of data_type (in wiggle files matrix columns = data type + 2)
        //         read count of the chromosome region
        //         boolean if peak should be treated as log ratio
        void binary_search(const int chrom, const int chrom_begin, const int chromend, const int starting_point, const int end_point, const int data_type, const double peak, bool is_log);
        void binary_search_mean(const int chrom, const int chrom_begin, const int chrom_end, const int start_point, const int end_point, const int data_type, const double value);

        // initialize the mapping of external names for chromosomes used in files to internal
        // this is required to avoid problems matching different notations for the same chromosome in files of a single run
        void init_chr_mapping();


        // matrix that holds data that have been read
        // for wiggle files this should be:
        //      for each genome region a different line containing
        //          genome_region_start   genome_region_end   overlap_assay_1  overlap_assay_2 ...
        Matrix<double> matrix_;

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
