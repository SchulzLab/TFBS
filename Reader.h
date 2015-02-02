#ifndef READER_H
#define READER_H

#include <string>
#include <stdio.h> // fscanf
#include <utility> // pair

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
        void read_file(const string& file_path, int data_type);

        // reads a peak file in wiggle format
        //
        // @param: path to a peak file
        void read_peak_file(const string& peak_file_path);

        // get matrix - i.e. get previously read data
        Matrix<float>& get_prev_read_data();


    private:


        // store data in matrix
        // do overlap computation if genome region corresponds to previously read data
        //
        // called by   read_file
        //
        // @param:  chromosome region start and end point
        //          measured peak
        //          number of data_type (in wiggle files matrix columns = data type + 3)
        //          if peak_file is provided, -1 is data type for peak file
        void store_data(const int chrom, const int chrom_begin, const int chrom_end, const float peak, const int data_type);



        // search in matrix for given genome region, compute overlap and save peak
        //
        // called by   store_data
        //
        // @param: chromosome region start and end point
        //         starting point for binary search
        //         number of data_type (in wiggle files matrix columns = data type + 2)
        //         peak of the chromosome region
        void binary_search(const int chrom, const int chrom_begin, const int chromend, const int starting_point, const int data_type, const float peak);



        // since we have a strongly modified version of binary search we need
        // to remember if we've just did a right or left jump for the binary search
        // to avoid endless recursive descent
        bool gone_right = false;
        bool gone_left = false;

        // matrix that holds data that have been read
        // for wiggle files this should be:
        //      for each genome region a different line containing
        //          genome_region_start   genome_region_end   peak_test_1  peak_test_2 ...
        Matrix<float> matrix_;

        // number of different data types
        // i.e. number of columns in matrix_
        int number_of_data_types_;

        // specifies if peak file was provided
        bool has_peak_file_;

        // last found position of binary search
        // this speeds up the next binary search because the genome regions
        // are ordered within the files
        int last_found_pos_;

        // last used iterator for the lines of the matrix used by binary
        // search if not yet mapped regions are found in new data types
        list<vector<float>>::iterator last_found_it_;

        // file on which the reader actually works on
        FILE* actual_file_;
};

#endif /* READER_H */
