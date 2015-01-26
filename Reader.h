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

        // get matrix - i.e. get previously read data
        Matrix<int>& get_prev_read_data();


    private:


        // store data in matrix
        // do overlap computation if genome region corresponds to previously read data
        //
        // called by   read_file
        //
        // @param:  chromosome region start and end point
        //          measured peak
        //          number of data_type (in wiggle files matrix columns = data type + 2)
        void store_data(const int chrom_begin, const int chrom_end, const int peak, const int data_type);



        // search in matrix for given genome region, compute overlap and save peak
        //
        // called by   store_data
        //
        // @param: chromosome region start and end point
        //         starting point for binary search
        //         number of data_type (in wiggle files matrix columns = data type + 2)
        //         peak of the chromosome region
        void binary_search(const int chrom_begin, const int chromend, const int starting_point, const int data_type, const int peak);



        // matrix that holds data that have been read
        // for wiggle files this should be:
        //      for each genome region a different line containing
        //          genome_region_start   genome_region_end   peak_test_1  peak_test_2 ...
        Matrix<int> matrix_;

        // number of different data types
        // i.e. number of columns in matrix_
        int number_of_data_types_;

        // last found position of binary search
        // this speeds up the next binary search because the genome regions
        // are ordered within the files
        int last_found_pos_;

        // last used iterator for the lines of the matrix used by binary
        // search if not yet mapped regions are found in new data types
        list<vector<int>>::iterator last_found_it_;

        // file on which the reader actually works on
        FILE* actual_file_;
};

#endif /* READER_H */
