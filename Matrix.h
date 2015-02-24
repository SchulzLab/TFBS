#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

template <typename T>
class Matrix {

    public:


        // -- Constructors --

        Matrix();

        Matrix(int lines, int columns, const T& initial_value);

        Matrix(const Matrix& M) = default;

        // move constructor
        Matrix(Matrix&& M);

        // ---



        // move assignment operator
        Matrix& operator=(Matrix&& M);

        // get the potential number of lines
        // i.e. space reserved to fill without the need of reallocation
        int get_line_capacity();

        int get_number_of_lines();
        int get_number_of_columns();


        // element access
        T& operator() (int line, int colnum);

        vector<T>& get_column(int colnum);




        // appends a new column to the matrix
        //
        // @param: the new column that should be appended
        void append_new_column(const vector<T>& new_line);

        // returns a proper representation of the matrix
        friend ostream& operator<<(ostream& os, const Matrix& M) {

            os << endl << endl;
            if (M.matrix_.size() == 0) {
                os << "<Empty>" << endl;
                return os;
            } else {
                if (M.matrix_.front().size() == 0) {
                    os << "<Empty>" << endl;
                    return os;
                }
            }

            // iterate over lines
            for (int lines = 0; lines < M.matrix_[0].size(); ++lines) {
                // iterate over columns
                for (int columns = 0; columns < M.matrix_.size(); ++columns) {
                    os << setw(5) << right << M.matrix_[columns][lines] << "\t";
                }
                os << "\n\n\n";
            }
            return os;
        }

    private:

        // matrix columns<lines>
        vector<vector<T> > matrix_;
};



// --- template definitions ---

template <typename T> Matrix<T>::Matrix() :
        matrix_()
{
}





template <typename T> Matrix<T>::Matrix(int lines, int columns, const T& initial_value) :
        matrix_(columns, vector<T>(lines, initial_value))
{
}





template <typename T> Matrix<T>::Matrix(Matrix&& M) :
        matrix_(move(M.matrix_))
{
    M.matrix_ = vector<vector<T>>();
}





template <typename T> Matrix<T>& Matrix<T>::operator=(Matrix&& M) {

    matrix_ = move(M.matrix_);
    return *this;
}





template <typename T> int Matrix<T>::get_number_of_lines() {

    return matrix_.size() > 0 ? matrix_[0].size() : 0;
}





template <typename T> int Matrix<T>::get_number_of_columns() {

    return matrix_.size();
}





template <typename T> T& Matrix<T>::operator()(int linenum, int colnum) {

    return matrix_[colnum][linenum];
}





template <typename T> vector<T>& Matrix<T>::get_column(int colnum) {

    return matrix_[colnum];
}





template <typename T> void Matrix<T>::append_new_column(const vector<T>& new_column) {

    matrix_.push_back(new_column);
}

#endif /* MATRIX_H */
