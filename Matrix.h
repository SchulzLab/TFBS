#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <list>
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
        // element access with iterator for lines - better performance
        T& operator() (const typename list<vector<T>>::iterator line_it, int colnum);

        // appends a new line to the matrix
        //
        // @param: the new line that should be appended
        void append_new_line(const vector<T>& new_line);

        // inserts new line at specified position
        //
        // @ param: position (or iterator to position) where the new line should be inserted
        //          content for the new line
        typename list<vector<T>>::iterator insert_new_line(const typename list<vector<T>>::iterator pos, const vector<T>& new_line);
        typename list<vector<T>>::iterator insert_new_line(int pos, const vector<T>& new_line);

        // access for the last line
        //
        // @return: iterator pointing to the last line of matrix_
        typename list<vector<T>>::iterator last_line();

        // access for the first line
        //
        // @return: iterator pointing to the first line of matrix_
        typename list<vector<T>>::iterator first_line();

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
            for ( vector<T> line : M.matrix_ ) {
                // iterate over columns
                for (int j = 0; j < line.size(); ++j) {
                    os << setw(5) << right << line[j] << "\t";
                }
                os << "\n\n\n";
            }
            return os;
        }

    private:

        // matrix line -> column access
        list<vector<T> > matrix_;
};



// --- template definitions ---

template <typename T> Matrix<T>::Matrix() :
        matrix_()
{
}





template <typename T> Matrix<T>::Matrix(int lines, int columns, const T& initial_value)
{
    for (int i = 0; i < lines; ++i) {

        matrix_.push_back(vector<T>(columns, initial_value));
    }
}





template <typename T> Matrix<T>::Matrix(Matrix&& M) :
        matrix_(move(M.matrix_))
{
    M.matrix_ = list<vector<T>>();
}





template <typename T> Matrix<T>& Matrix<T>::operator=(Matrix&& M) {

    matrix_ = move(M.matrix_);
    return *this;
}





template <typename T> int Matrix<T>::get_number_of_lines() {

    return matrix_.size();
}





template <typename T> int Matrix<T>::get_number_of_columns() {

    return matrix_.front().size();
}





template <typename T> T& Matrix<T>::operator()(int linenum, int colnum) {

    // visit backwards
    if (linenum > matrix_.size() / 2) {

        // we start counting from one. Sadly size() doesn't.
        int counter = matrix_.size() - linenum - 1;
        auto line_it = matrix_.end();
        while (counter >= 0) {

            --line_it;
            --counter;
        }
        return (*line_it)[colnum];

    // visit forwards
    } else {

        int counter = linenum;
        auto line_it = matrix_.begin();
        while (counter > 0) {

            ++line_it;
            --counter;
        }
        return (*line_it)[colnum];
    }
}





template <typename T> T& Matrix<T>::operator() (const typename list<vector<T>>::iterator line_it, int colnum) {

    return (*line_it)[colnum];
}





template <typename T> void Matrix<T>::append_new_line(const vector<T>& new_line) {

    matrix_.push_back(new_line);
}





template <typename T> typename list<vector<T>>::iterator Matrix<T>::insert_new_line(const typename list<vector<T>>::iterator pos, const vector<T>& new_line) {

    return matrix_.insert(pos, new_line);
}





template <typename T> typename list<vector<T>>::iterator Matrix<T>::insert_new_line(int pos, const vector<T>& new_line) {

    typename list<vector<T>>::iterator line_it;

    // visit backwards
    if (pos > matrix_.size() / 2) {

        // we start counting from one. Sadly size() doesn't.
        int counter = matrix_.size() - pos - 1;
        line_it = matrix_.end();
        while (counter >= 0) {

            --line_it;
            --counter;
        }
        matrix_.insert(line_it, new_line);

    // visit forwards
    } else {

        int counter = pos;
        line_it = matrix_.begin();
        while (counter > 0) {

            ++line_it;
            --counter;
        }
        matrix_.insert(line_it, new_line);
    }
    return line_it;
}





template <typename T> typename list<vector<T>>::iterator Matrix<T>::first_line() {

    return matrix_.begin();
}





template <typename T> typename list<vector<T>>::iterator Matrix<T>::last_line() {

    return --matrix_.end();
}

#endif /* MATRIX_H */
