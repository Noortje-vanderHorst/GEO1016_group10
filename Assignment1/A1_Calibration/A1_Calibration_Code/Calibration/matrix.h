
#ifndef _MATH_MATRIX_H_
#define _MATH_MATRIX_H_

#include <cassert>

#include <vector>
#include <complex>
#include <algorithm> // for std::min, std::max
#include <cmath>	 // for std::sqrt


/**
 * A templated matrix class which is designed for basic linear algebra operations such as:
 *              A + x    x + A    A += x    A1 + A2   A1 += A2
 *              A - x    x - A    A -= x    A1 - A2   A1 -= A2
 *              A * x    x * A    A *= x    A1 * A2
 *              A / x    x / A    A /= x
 *              sum,	min, max, mean, swap, norm, identity, diagonal,
 *			    transpose,  conjugate_transpose
 *              transpose_mult,   mult_transpose
 *			    elem_mult, elem_mult_Eq, elem_divd, elem_divd_Eq
 * These operators and functions can be applied to both real matrix and complex matrices.
 *
 * @tparam FT The data type, which can be float, double, or std::complex.
 */

namespace easy3d {

    template<typename FT>
    class Matrix {
    public:
        // constructors and destructor
        Matrix();

        /**
         * Copy constructor (i.e., construct a matrix from another matrix). The new matrix will be the same as the input
         * matrix.
         * @param A Input matrix.
         */
        Matrix(const Matrix<FT> &A);

        /**
         * Construct a matrix by specifying its dimension (i.e., number of rows and number of columns).
         * @param rows Number of rows
         * @param cols Number of columns.
         * @param x Default value of each element of the matrix.
         */
        Matrix(int rows, int cols, const FT &x = FT(0));

        /**
         * Construct a matrix by specifying its dimension (i.e., number of rows and number of columns) and initialize the
         * matrix from a 1D array.
         * @param rows Number of rows
         * @param cols Number of columns.
         * @param v Default values the matrix elements.
         * Attention: Matrices are stored internally as row-major.
         */
        Matrix(int rows, int cols, const FT *v);

        ~Matrix();

        // assignments
        Matrix<FT> &operator=(const Matrix<FT> &A);    // overload evaluate operator = from matrix to matrix
        Matrix<FT> &operator=(const FT &x);            // overload evaluate operator = from scalar to matrix

        // Return the row_th row as a 1D array.
        FT *operator[](int row);

        // Return the row_th row as a 1D array.
        const FT *operator[](int row) const;

        // Return the element at (row, col).
        FT &operator()(int row, int col);

        // Return the element at (row, col).
        const FT &operator()(int row, int col) const;

        // Set the value of element at (row, col).
        void set(int row, int col, FT v);

        // Get the value of element at (row, col).
        const FT &get(int row, int col) const;

        /**
         * Return the data of the matrix as a 1D array.
         * Attention: Matrices are stored internally as row-major.
         */
        operator FT *();
        operator const FT *() const;

        // The number of rows
        int rows() const;

        // The number of cols
        int cols() const;

        // Reallocate matrix's size
        Matrix<FT> &resize(int rows, int cols);

        // get the matrix's row vector as a 1D array
        std::vector<FT> get_row(int row) const;

        // get the matrix's column vector
        std::vector<FT> get_column(int col) const;

        // Set the matrix's row vector
        void set_row(const std::vector<FT> &v, int row);

        // Set the matrix's column vector
        void set_column(const std::vector<FT> &v, int col);

        // Set all elements to zero.
        void load_zero();

        // Make this matrix identity (i.e., all elements on the diagonal have a value of 1). This function also allow to set
        // the elements on the diagonal to have values other than 1.
        void load_identity(FT v = FT(1));

        // Return the transposed matrix.
        Matrix<FT> transpose() const;

        // Return= the trace of this matrix, i.e. the sum of the coefficients on the main diagonal.
        // NOTE: the matrix can be any matrix, not necessarily square.
        FT trace() const;

        // computed assignment
        Matrix<FT> &operator+=(const FT &);    //compound assignment operators +=
        Matrix<FT> &operator-=(const FT &);    //compound assignment operators -=
        template<typename FT2>
        Matrix<FT> &operator*=(const FT2 &);    //compound assignment operators *=
        template<typename FT2>
        Matrix<FT> &operator/=(const FT2 &);    //compound assignment operators /=

        Matrix<FT> &operator+=(const Matrix<FT> &);

        Matrix<FT> &operator-=(const Matrix<FT> &);

        Matrix<FT> &operator*=(const Matrix<FT> &);    // WARNING: element-wise multiplication
        Matrix<FT> &operator/=(const Matrix<FT> &);    // WARNING: element-wist division

    private:

        // 0-based data pointer
        FT *data_;

        // 0-based row pointer's pointer
        FT **prow_;

        // row number, column number and total number
        int nRow;
        int nColumn;
        long nTotal;

        void init(int rows, int cols);

        void copy_from_array(const FT *v); // copy matrix from normal array
        void set_by_scalar(const FT &x); // set matrix by a scalar
        void destroy();    // destroy the matrix

    }; // class Matrix



//----------------------------- Utilities for Matrix ------------------------------//

// input and output
    template<typename FT>
    std::ostream &operator<<(std::ostream &, const Matrix<FT> &);    // Overload of the output stream.
    template<typename FT>
    std::istream &operator>>(std::istream &, Matrix<FT> &);        // Overload of the input stream.

// arithmetic operators
    template<typename FT>
    Matrix<FT> operator-(const Matrix<FT> &);  // get negative matrix
    template<typename FT>
    Matrix<FT> operator+(const Matrix<FT> &, const FT &); // matrix-scalar addition
    template<typename FT>
    Matrix<FT> operator+(const FT &, const Matrix<FT> &);    // scalar-matrix addition
    template<typename FT>
    Matrix<FT> operator+(const Matrix<FT> &, const Matrix<FT> &); // matrix-matrix addition
    template<typename FT>
    Matrix<FT> operator-(const Matrix<FT> &, const FT &); // matrix-scalar subtraction
    template<typename FT>
    Matrix<FT> operator-(const FT &, const Matrix<FT> &);    // scalar-matrix subtraction
    template<typename FT>
    Matrix<FT> operator-(const Matrix<FT> &, const Matrix<FT> &); // matrix-matrix subtraction
    template<typename FT>
    Matrix<FT> operator*(const Matrix<FT> &, const Matrix<FT> &);        // matrix-matrix multiplication
    template<typename FT>
    std::vector<FT> operator*(const Matrix<FT> &, const std::vector<FT> &);    // matrix-vector multiplication
    template<typename FT1, typename FT2>
    Matrix<FT1> operator*(const Matrix<FT1> &, const FT2 &);    // matrix-scalar multiplication
    template<typename FT1, typename FT2>
    Matrix<FT1> operator*(const FT2 &, const Matrix<FT1> &);    // scalar-matrix multiplication
    template<typename FT1, typename FT2>
    Matrix<FT1> operator/(const Matrix<FT1> &, const FT2 &);    // matrix-scalar division
    template<typename FT1, typename FT2>
    Matrix<FT1> operator/(const FT2 &, const Matrix<FT1> &);    // scalar-matrix division


// multiplication
    template<typename FT>
    void mult(const Matrix<FT> &, const Matrix<FT> &,
              Matrix<FT> &);            // matrix-matrix multiplication (result has already been allocated)
    template<typename FT>
    void mult(const Matrix<FT> &, const std::vector<FT> &,
              std::vector<FT> &);    // matrix-vector multiplication (result has already been allocated)
    template<typename FT>
    Matrix<FT> mult(const Matrix<FT> &, const Matrix<FT> &);        // matrix-matrix multiplication
    template<typename FT>
    std::vector<FT> mult(const Matrix<FT> &, const std::vector<FT> &);    // matrix-vector multiplication

// element-wise multiplication / division
    template<typename FT>
    Matrix<FT> elem_mult(const Matrix<FT> &, const Matrix<FT> &);    // "*"
    template<typename FT>
    Matrix<FT> elem_divd(const Matrix<FT> &, const Matrix<FT> &);    // "/"
    template<typename FT>
    Matrix<FT> &elem_mult_eq(Matrix<FT> &, const Matrix<FT> &);        // "*="
    template<typename FT>
    Matrix<FT> &elem_divd_eq(Matrix<FT> &, const Matrix<FT> &);        // "/="

// transpose and conjugate transpose
    template<typename FT>
    Matrix<FT> transpose(const Matrix<FT> &);            // matrix transpose
    template<typename FT>
    Matrix<FT> conjugate_transpose(const Matrix<FT> &);    // matrix conjugate transpose

// transpose multiplication
    template<typename FT>
    Matrix<FT> transpose_mult(const Matrix<FT> &, const Matrix<FT> &);        // A^T * B
    template<typename FT>
    std::vector<FT> transpose_mult(const Matrix<FT> &, const std::vector<FT> &);    // A^T * b
    template<typename FT>
    Matrix<FT> mult_transpose(const Matrix<FT> &, const Matrix<FT> &);        // A * B^T
    template<typename FT>
    Matrix<FT> mult_transpose(const std::vector<FT> &, const std::vector<FT> &);    // a * b^T

// unit and diagonal matrix
    template<typename FT>
    Matrix<FT> identity(int, const FT &);            // Generate an identity matrix.
    template<typename FT>
    Matrix<FT> diagonal(const std::vector<FT> &);    // Generate a diagonal matrix by given its diagonal elements.
    template<typename FT>
    std::vector<FT> diagonal(const Matrix<FT> &);    // Get the diagonal entries of matrix.

// the trace of this matrix, i.e. the sum of the coefficients on the main diagonal.
// NOTE: the matrix can be any matrix, not necessarily square.
    template<typename FT>
    FT trace(const Matrix<FT> &);

// utilities
    template<typename FT>
    FT norm(const Matrix<FT> &);                // Compute Frobenius norm of matrix.
    template<typename FT>
    void swap(Matrix<FT> &, Matrix<FT> &);    // Swap two matrices.
    template<typename FT>
    std::vector<FT> sum(const Matrix<FT> &);    // Matrix's column vectors sum.
    template<typename FT>
    std::vector<FT> min(const Matrix<FT> &);    // Minimum of matrix's column vectors.
    template<typename FT>
    std::vector<FT> max(const Matrix<FT> &);    // Maximum of matrix's column vectors.
    template<typename FT>
    std::vector<FT> mean(const Matrix<FT> &);    // Matrix's column vectors mean.
    template<typename FT>
    Matrix<FT> abs(const Matrix<std::complex<FT> > &A);// Get magnitude of a std::complex matrix.
    template<typename FT>
    Matrix<FT> arg(const Matrix<std::complex<FT> > &A);// Get angle of a std::complex matrix.
    template<typename FT>
    Matrix<FT> real(const Matrix<std::complex<FT> > &A);// Get real part of a std::complex matrix.
    template<typename FT>
    Matrix<FT> imag(const Matrix<std::complex<FT> > &A);// Get imaginary part of a std::complex matrix.

// Convert real matrix to complex matrix.
    template<typename FT>
    Matrix<std::complex<FT> > complex_matrix(const Matrix<FT> &);

    template<typename FT>
    Matrix<std::complex<FT> > complex_matrix(const Matrix<FT> &, const Matrix<FT> &); // A for real, B for imaginary


//----------------------------- Utilities for std::vector ------------------------------//

// input and output
    template<typename FT>
    std::ostream &operator<<(std::ostream &, const std::vector<FT> &);// Overload of the output stream.
    template<typename FT>
    std::istream &operator>>(std::istream &, std::vector<FT> &);        // Overload of the input stream.

// arithmetic operators
    template<typename FT>
    std::vector<FT> operator-(const std::vector<FT> &);    // get negative vector
    template<typename FT>
    std::vector<FT> operator+(const std::vector<FT> &, const std::vector<FT> &);    // vector-vector addition
    template<typename FT>
    std::vector<FT> operator+(const std::vector<FT> &, const FT &);    // vector-scalar addition
    template<typename FT>
    std::vector<FT> operator+(const FT &, const std::vector<FT> &);    // scalar-vector addition
    template<typename FT>
    std::vector<FT> operator-(const std::vector<FT> &, const std::vector<FT> &);    // vector-vector subtraction
    template<typename FT>
    std::vector<FT> operator-(const std::vector<FT> &, const FT &);    // vector-scalar subtraction
    template<typename FT>
    std::vector<FT> operator-(const FT &, const std::vector<FT> &);    // scalar-vector subtraction
    template<typename FT1, typename FT2>
    std::vector<FT1> operator*(const std::vector<FT1> &, const FT2 &);    // complex vector-scalar multiplication
    template<typename FT1, typename FT2>
    std::vector<FT1> operator*(const FT2 &, const std::vector<FT1> &);    // scalar-complex vector multiplication
    template<typename FT1, typename FT2>
    std::vector<FT1> operator/(const std::vector<FT1> &, const FT2 &);    // complex vector-scalar division
    template<typename FT1, typename FT2>
    std::vector<FT1> operator/(const FT2 &, const std::vector<FT1> &);    // scalar-complex vector division


// element-wise multiplication / division
    template<typename FT>
    std::vector<FT> elem_mult(const std::vector<FT> &, const std::vector<FT> &);    // "*"
    template<typename FT>
    std::vector<FT> elem_divd(const std::vector<FT> &, const std::vector<FT> &);    // "/"
    template<typename FT>
    std::vector<FT> &elem_mult_eq(std::vector<FT> &, const std::vector<FT> &);    // "*="
    template<typename FT>
    std::vector<FT> &elem_divd_eq(std::vector<FT> &, const std::vector<FT> &);    // "/="

// dot product
    template<typename FT>
    FT operator*(const std::vector<FT> &, const std::vector<FT> &);    // Inner product for vectors.
    template<typename FT>
    FT dot(const std::vector<FT> &, const std::vector<FT> &);            // Inner product for vectors.
    /// custom
    template<typename FT>
    std::vector<FT> cross(const std::vector<FT> &, const std::vector<FT> &);       // Outer product (cross) for vectors.

// utilities
    template<typename FT>
    FT norm(const std::vector<FT> &);                    // Euclidean norm.
    template<typename FT>
    FT norm(const std::vector<std::complex<FT> > &);    // Euclidean norm.
    template<typename FT>
    void swap(std::vector<FT> &, std::vector<FT> &);// Swap two vectors.
    template<typename FT>
    std::vector<FT>
    linspace(FT, FT, int);// Generates a vector of n points linearly spaced between and including a and b.
    template<typename FT>
    FT sum(const std::vector<FT> &);    // vector sum.
    template<typename FT>
    FT min(const std::vector<FT> &);    // Minimum value of vector.
    template<typename FT>
    FT max(const std::vector<FT> &);    // Maximum value of vector.
    template<typename FT>
    FT mean(const std::vector<FT> &);

    template<typename FT>
    std::vector<FT> abs(const std::vector<std::complex<FT> > &);    // Get magnitude of a complex vector.
    template<typename FT>
    std::vector<FT> arg(const std::vector<std::complex<FT> > &);    // Get angle of a complex vector.
    template<typename FT>
    std::vector<FT> real(const std::vector<std::complex<FT> > &);    // Get real part of a complex vector.
    template<typename FT>
    std::vector<FT> imag(const std::vector<std::complex<FT> > &);    // Get imaginary part of a complex vector.

// Convert real vector to complex vector.
    template<typename FT>
    std::vector<std::complex<FT> > complex_vector(const std::vector<FT> &);

    template<typename FT>
    std::vector<std::complex<FT> >
    complex_vector(const std::vector<FT> &, const std::vector<FT> &); // A for real, B for imaginary


//-----------------------------  Matrix  Implementation ------------------------------------//

#include <cassert>

/**
* initialize
*/
    template<typename FT>
    void Matrix<FT>::init(int rows, int cols) {
        nRow = rows;
        nColumn = cols;
        nTotal = nRow * nColumn;

        data_ = new FT[nTotal];
        prow_ = new FT *[nRow];

        assert(data_ != NULL);
        assert(prow_ != NULL);

        FT *p = data_;
        for (int i = 0; i < nRow; ++i) {
            prow_[i] = p;
            p += nColumn;
        }
    }


/**
* copy matrix from normal array
*/
    template<typename FT>
    inline void Matrix<FT>::copy_from_array(const FT *v) {
        for (long i = 0; i < nTotal; ++i)
            data_[i] = v[i];
    }


/**
* set matrix by a scalar
*/
    template<typename FT>
    inline void Matrix<FT>::set_by_scalar(const FT &x) {
        for (long i = 0; i < nTotal; ++i)
            data_[i] = x;
    }


/**
* destroy the matrix
*/
    template<typename FT>
    void Matrix<FT>::destroy() {
        if (data_ == NULL)
            return;
        else
            delete[] data_;

        if (prow_ != NULL)
            delete[] prow_;
    }


/**
* constructors and destructor
*/
    template<typename FT>
    Matrix<FT>::Matrix()
            : data_(0), prow_(0), nRow(0), nColumn(0), nTotal(0) {
    }

    template<typename FT>
    Matrix<FT>::Matrix(const Matrix<FT> &A) {
        init(A.nRow, A.nColumn);
        copy_from_array(A.data_);
    }

    template<typename FT>
    Matrix<FT>::Matrix(int rows, int cols, const FT &x) {
        init(rows, cols);
        set_by_scalar(x);
    }

    template<typename FT>
    Matrix<FT>::Matrix(int rows, int cols, const FT *arrays) {
        init(rows, cols);
        copy_from_array(arrays);
    }

    template<typename FT>
    Matrix<FT>::~Matrix() {
        destroy();
    }


/**
* overload evaluate operator = from matrix to matrix
*/
    template<typename FT>
    Matrix<FT> &Matrix<FT>::operator=(const Matrix<FT> &A) {
        if (data_ == A.data_)
            return *this;

        if (nRow == A.nRow && nColumn == A.nColumn)
            copy_from_array(A.data_);
        else {
            destroy();
            init(A.nRow, A.nColumn);
            copy_from_array(A.data_);
        }

        return *this;
    }


/**
* overload evaluate operator = from scalar to matrix
*/
    template<typename FT>
    inline Matrix<FT> &Matrix<FT>::operator=(const FT &x) {
        set_by_scalar(x);
        return *this;
    }


/**
* overload operator [] for 0-offset access
*/
    template<typename FT>
    inline FT *Matrix<FT>::operator[](int row) {
        assert(0 <= row);
        assert(row < nRow);
        return prow_[row];
    }

    template<typename FT>
    inline const FT *Matrix<FT>::operator[](int row) const {
        assert(0 <= row);
        assert(row < nRow);
        return prow_[row];
    }

    template<typename FT>
    inline FT &Matrix<FT>::operator()(int row, int col) {
        assert(0 <= row);
        assert(row < nRow);
        assert(0 <= col);
        assert(col < nColumn);
        return prow_[row][col];
    }

    template<typename FT>
    inline const FT &Matrix<FT>::operator()(int row, int col) const {
        assert(0 <= row);
        assert(row < nRow);
        assert(0 <= col);
        assert(col < nColumn);
        return prow_[row][col];
    }


    template<typename FT>
    inline
    void Matrix<FT>::set(int row, int col, FT v) {
        assert(0 <= row);
        assert(row < nRow);
        assert(0 <= col);
        assert(col < nColumn);
        prow_[row][col] = v;
    }

    template<typename FT>
    inline
    const FT &Matrix<FT>::get(int row, int col) const {
        assert(0 <= row);
        assert(row < nRow);
        assert(0 <= col);
        assert(col < nColumn);
        return prow_[row][col];
    }


/**
* Clears the matrix
* This resets all values to 0 (zero)
*/
    template<typename FT>
    inline
    void Matrix<FT>::load_zero() {
        for (int i = 0; i < nRow; i++) {
            for (int j = 0; j < nColumn; j++) {
                prow_[i][j] = FT(0);
            }
        }
    }

/**
* Sets the matrix to identity
* This sets all coefficients of this matrix to be equal to
* nROW x nCOL identity matrix.
*/
    template<typename FT>
    inline
    void Matrix<FT>::load_identity(FT v /* = FT(1.0)*/) {
        for (int i = 0; i < nRow; i++) {
            for (int j = 0; j < nColumn; j++) {
                prow_[i][j] = (i == j) ? v : FT(0);
            }
        }
    }


/**
* type conversion functions
*/
    template<typename FT>
    inline Matrix<FT>::operator FT *() {
        return data_;
    }

    template<typename FT>
    inline Matrix<FT>::operator const FT *() const {
        return data_;
    }


    template<typename FT>
    inline int Matrix<FT>::rows() const {
        return nRow;
    }

    template<typename FT>
    inline int Matrix<FT>::cols() const {
        return nColumn;
    }


/**
* reallocate matrix's size
*/
    template<typename FT>
    Matrix<FT> &Matrix<FT>::resize(int rows, int cols) {
        if (rows == nRow && cols == nColumn)
            return *this;

        destroy();
        init(rows, cols);

        return *this;
    }


/**
* get the matrix's row vector
*/
    template<typename FT>
    std::vector<FT> Matrix<FT>::get_row(int row) const {
        assert(0 <= row);
        assert(row < nRow);
        std::vector<FT> tmp(nColumn);
        for (int j = 0; j < nColumn; ++j)
            tmp[j] = prow_[row][j];

        return tmp;
    }


/**
* get the matrix's column vector
*/
    template<typename FT>
    std::vector<FT> Matrix<FT>::get_column(int col) const {
        assert(0 <= col);
        assert(col < nColumn);
        std::vector<FT> tmp(nRow);
        for (int i = 0; i < nRow; ++i)
            tmp[i] = prow_[i][col];

        return tmp;
    }


/**
* set the matrix's row vector
*/
    template<typename FT>
    void Matrix<FT>::set_row(const std::vector<FT> &v, int row) {
        assert(0 <= row);
        assert(row < nRow);
        assert(v.size() == nColumn);
        for (int j = 0; j < nColumn; ++j)
            prow_[row][j] = v[j];
    }


/**
* set the matrix's column vector
*/
    template<typename FT>
    void Matrix<FT>::set_column(const std::vector<FT> &v, int col) {
        assert(0 <= col);
        assert(col < nColumn);
        assert(v.size() == nRow);
        for (int i = 0; i < nRow; ++i)
            prow_[i][col] = v[i];
    }


    template<typename FT>
    Matrix<FT> Matrix<FT>::transpose() const {
        Matrix<FT> t(nColumn, nRow);
        for (int r = 0; r < nRow; r++) {
            for (int c = 0; c < nColumn; c++) {
                t(c, r) = (*this)(r, c);
            }
        }
        return t;
    }


    template<typename FT>
    FT Matrix<FT>::trace() const {
        int range = std::min(nRow, nColumn);
        FT trac = 0.0;
        for (int i = 0; i < range; i++) {
            trac += (*this)[i][i];
        }
        return trac;
    }

/**
* compound assignment operators +=
*/
    template<typename FT>
    Matrix<FT> &Matrix<FT>::operator+=(const FT &x) {
        FT **rowPtr = prow_;
        FT *colPtr = 0;

        for (int i = 0; i < nRow; ++i) {
            colPtr = *rowPtr++;
            for (int j = 0; j < nColumn; ++j)
                *colPtr++ += x;
        }

        return *this;
    }

    template<typename FT>
    Matrix<FT> &Matrix<FT>::operator+=(const Matrix<FT> &rhs) {
        assert(nRow == rhs.rows());
        assert(nColumn == rhs.cols());

        FT **rowPtrL = prow_;
        FT *colPtrL = 0;
        FT **rowPtrR = rhs.prow_;
        const FT *colPtrR = 0;

        for (int i = 0; i < nRow; ++i) {
            colPtrL = *rowPtrL++;
            colPtrR = *rowPtrR++;
            for (int j = 0; j < nColumn; ++j)
                *colPtrL++ += *colPtrR++;
        }

        return *this;
    }


/**
* compound assignment operators -=
*/
    template<typename FT>
    Matrix<FT> &Matrix<FT>::operator-=(const FT &x) {
        FT **rowPtr = prow_;
        FT *colPtr = 0;

        for (int i = 0; i < nRow; ++i) {
            colPtr = *rowPtr++;
            for (int j = 0; j < nColumn; ++j)
                *colPtr++ -= x;
        }

        return *this;
    }

    template<typename FT>
    Matrix<FT> &Matrix<FT>::operator-=(const Matrix<FT> &rhs) {
        assert(nRow == rhs.rows());
        assert(nColumn == rhs.cols());

        FT **rowPtrL = prow_;
        FT *colPtrL = 0;
        FT **rowPtrR = rhs.prow_;
        const FT *colPtrR = 0;

        for (int i = 0; i < nRow; ++i) {
            colPtrL = *rowPtrL++;
            colPtrR = *rowPtrR++;
            for (int j = 0; j < nColumn; ++j)
                *colPtrL++ -= *colPtrR++;
        }

        return *this;
    }


/**
* compound assignment operators *=
*/
    template<typename FT>
    template<typename FT2>
    Matrix<FT> &Matrix<FT>::operator*=(const FT2 &x) {
        FT **rowPtr = prow_;
        FT *colPtr = 0;

        for (int i = 0; i < nRow; ++i) {
            colPtr = *rowPtr++;
            for (int j = 0; j < nColumn; ++j)
                *colPtr++ *= x;
        }

        return *this;
    }

// WARNING: this is element-by-element multiplication
    template<typename FT>
    Matrix<FT> &Matrix<FT>::operator*=(const Matrix<FT> &rhs) {
        assert(nRow == rhs.rows());
        assert(nColumn == rhs.cols());

        FT **rowPtrL = prow_;
        FT *colPtrL = 0;
        FT **rowPtrR = rhs.prow_;
        const FT *colPtrR = 0;

        for (int i = 0; i < nRow; ++i) {
            colPtrL = *rowPtrL++;
            colPtrR = *rowPtrR++;
            for (int j = 0; j < nColumn; ++j)
                *colPtrL++ *= *colPtrR++;
        }

        return *this;
    }


/**
* compound assignment operators /=
*/
    template<typename FT>
    template<typename FT2>
    Matrix<FT> &Matrix<FT>::operator/=(const FT2 &x) {
        FT **rowPtr = prow_;
        FT *colPtr = 0;

        for (int i = 0; i < nRow; ++i) {
            colPtr = *rowPtr++;
            for (int j = 0; j < nColumn; ++j)
                *colPtr++ /= x;
        }

        return *this;
    }

// WARNING: this is element-by-element division
    template<typename FT>
    Matrix<FT> &Matrix<FT>::operator/=(const Matrix<FT> &rhs) {
        assert(nRow == rhs.rows());
        assert(nColumn == rhs.cols());

        FT **rowPtrL = prow_;
        FT *colPtrL = 0;
        FT **rowPtrR = rhs.prow_;
        const FT *colPtrR = 0;

        for (int i = 0; i < nRow; ++i) {
            colPtrL = *rowPtrL++;
            colPtrR = *rowPtrR++;
            for (int j = 0; j < nColumn; ++j)
                *colPtrL++ /= *colPtrR++;
        }

        return *this;
    }


/**
* Overload the output stream function.
*/
    template<typename FT>
    std::ostream &operator<<(std::ostream &out, const Matrix<FT> &A) {
        int rows = A.rows();
        int cols = A.cols();

        out << "size: " << rows << " by " << cols << "\n";
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                if (std::abs(A[i][j]) < 1e-6)
                    out << 0 << "\t";
                else
                    out << A[i][j] << "\t";
            }
            out << "\n";
        }

        return out;
    }


/**
* Overload the input stream function.
*/
    template<typename FT>
    std::istream &operator>>(std::istream &in, Matrix<FT> &A) {
        int rows, cols;
        in >> rows >> cols;

        if (!(rows == A.rows() && cols == A.cols()))
            A.resize(rows, cols);

        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                in >> A[i][j];

        return in;
    }


/**
* get negative matrix
*/
    template<typename FT>
    Matrix<FT> operator-(const Matrix<FT> &A) {
        int rows = A.rows();
        int cols = A.cols();

        Matrix<FT> tmp(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                tmp[i][j] = -A[i][j];

        return tmp;
    }


/**
* matrix-scalar addition
*/
    template<typename FT>
    inline Matrix<FT> operator+(const Matrix<FT> &A, const FT &x) {
        Matrix<FT> tmp(A);
        return tmp += x;
    }

    template<typename FT>
    inline Matrix<FT> operator+(const FT &x, const Matrix<FT> &A) {
        return A + x;
    }


/**
* matrix-matrix addition
*/
    template<typename FT>
    inline Matrix<FT> operator+(const Matrix<FT> &A1, const Matrix<FT> &A2) {
        Matrix<FT> tmp(A1);
        return tmp += A2;
    }


/**
* matrix-scalar subtraction
*/
    template<typename FT>
    inline Matrix<FT> operator-(const Matrix<FT> &A, const FT &x) {
        Matrix<FT> tmp(A);
        return tmp -= x;
    }

    template<typename FT>
    inline Matrix<FT> operator-(const FT &x, const Matrix<FT> &A) {
        Matrix<FT> tmp(A);
        return -tmp += x;
    }


/**
* matrix-matrix subtraction
*/
    template<typename FT>
    inline Matrix<FT> operator-(const Matrix<FT> &A1, const Matrix<FT> &A2) {
        Matrix<FT> tmp(A1);
        return tmp -= A2;
    }

/**
* matrix-matrix multiplication
*/
    template<typename FT>
    Matrix<FT> operator*(const Matrix<FT> &A1, const Matrix<FT> &A2) {
        assert(A1.cols() == A2.rows());

        int rows = A1.rows();
        int cols = A2.cols();
        //	int K = A1.cols();

        Matrix<FT> tmp(rows, cols);
        //	for( int i=0; i<rows; ++i )
        //		for( int j=0; j<cols; ++j )
        //		{
        //            tmp[i][j] = 0;
        //			for( int k=0; k<K; ++k )
        //			    tmp[i][j] += A1[i][k] * A2[k][j];
        //		}

        mult(A1, A2, tmp);

        return tmp;
    }


/**
* matrix-vector multiplication
*/
    template<typename FT>
    std::vector<FT> operator*(const Matrix<FT> &A, const std::vector<FT> &b) {
        assert(A.cols() == b.size());

        int rows = A.rows();
        //	int cols = A.cols();

        std::vector<FT> tmp(rows);
        mult(A, b, tmp);

        return tmp;
    }

// matrix-scalar multiplication
    template<typename FT1, typename FT2>
    inline
    Matrix<FT1> operator*(const Matrix<FT1> &A, const FT2 &s) {
        Matrix<FT1> tmp(A);
        return tmp *= s;
    }

// scalar-matrix multiplication
    template<typename FT1, typename FT2>
    inline
    Matrix<FT1> operator*(const FT2 &s, const Matrix<FT1> &A) {
        return A * s;
    }

// matrix-scalar division
    template<typename FT1, typename FT2>
    inline
    Matrix<FT1> operator/(const Matrix<FT1> &A, const FT2 &s) {
        Matrix<FT1> tmp(A);
        return tmp /= s;
    }

// scalar-matrix division
    template<typename FT1, typename FT2>
    inline
    Matrix<FT1> operator/(const FT2 &s, const Matrix<FT1> &A) {
        int rows = A.rows();
        int clumns = A.cols();

        Matrix<FT1> tmp(rows, clumns);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < clumns; ++j)
                tmp[i][j] = s / A[i][j];

        return tmp;
    }

/**
* This is an optimized version of matrix multiplication,
* where the destination matrix has already been allocated.
*/
    template<typename FT>
    void mult(const Matrix<FT> &A, const Matrix<FT> &B, Matrix<FT> &C) {
        int M = A.rows();
        int N = B.cols();
        int K = A.cols();

        assert(B.rows() == K);

        C.resize(M, N);
        FT sum = 0;
        const FT *pRow, *pCol;

        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; ++j) {
                pRow = &A[i][0];
                pCol = &B[0][j];
                sum = 0;

                for (int k = 0; k < K; ++k) {
                    sum += (*pRow) * (*pCol);
                    pRow++;
                    pCol += N;
                }
                C[i][j] = sum;
            }
    }


/**
* This is an optimized version of matrix and vector multiplication,
* where the destination vector has already been allocated.
*/
    template<typename FT>
    void mult(const Matrix<FT> &A, const std::vector<FT> &b, std::vector<FT> &c) {
        int M = A.rows();
        int N = A.cols();

        assert(b.size() == N);

        c.resize(M);
        FT sum = 0;
        const FT *pRow, *pCol;

        for (int i = 0; i < M; i++) {
            pRow = &A[i][0];
            pCol = &b[0];
            sum = 0;

            for (int j = 0; j < N; ++j) {
                sum += (*pRow) * (*pCol);
                pRow++;
                pCol++;
            }
            c[i] = sum;
        }
    }


/**
* This is an optimized version of matrix multiplication,
* where the destination matrix has already been allocated.
*/
    template<typename FT>
    Matrix<FT> mult(const Matrix<FT> &A, const Matrix<FT> &B) {
        int M = A.rows();
        int N = B.cols();
        int K = A.cols();

        assert(B.rows() == K);

        Matrix<FT> C(M, N);
        FT sum = 0;
        const FT *pRow, *pCol;

        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; ++j) {
                pRow = &A[i][0];
                pCol = &B[0][j];
                sum = 0;

                for (int k = 0; k < K; ++k) {
                    sum += (*pRow) * (*pCol);
                    pRow++;
                    pCol += N;
                }
                C[i][j] = sum;
            }
        return C;
    }


/**
* This is an optimized version of matrix and vector multiplication,
* where the destination vector has already been allocated.
*/
    template<typename FT>
    std::vector<FT> mult(const Matrix<FT> &A, const std::vector<FT> &b) {
        int M = A.rows();
        int N = A.cols();

        assert(b.size() == N);

        std::vector<FT> c(M);
        FT sum = 0;
        const FT *pRow, *pCol;

        for (int i = 0; i < M; i++) {
            pRow = &A[i][0];
            pCol = &b[0];
            sum = 0;

            for (int j = 0; j < N; ++j) {
                sum += (*pRow) * (*pCol);
                pRow++;
                pCol++;
            }
            c[i] = sum;
        }
        return c;
    }


/**
* matrix-matrix element-wise multiplication
*/
    template<typename FT>
    inline Matrix<FT> elem_mult(const Matrix<FT> &A1, const Matrix<FT> &A2) {
        Matrix<FT> tmp(A1);
        return tmp *= A2;
    }

    template<typename FT>
    inline Matrix<FT> &elem_mult_eq(Matrix<FT> &A1, const Matrix<FT> &A2) {
        return A1 *= A2;
    }


/**
* matrix-matrix element-wise division
*/
    template<typename FT>
    inline Matrix<FT> elem_divd(const Matrix<FT> &A1, const Matrix<FT> &A2) {
        Matrix<FT> tmp(A1);
        return tmp /= A2;
    }

    template<typename FT>
    inline Matrix<FT> &elem_divd_eq(Matrix<FT> &A1, const Matrix<FT> &A2) {
        return A1 /= A2;
    }


/**
* matrix transpose
*/
    template<typename FT>
    Matrix<FT> transpose(const Matrix<FT> &A) {
        int rows = A.cols();
        int clumns = A.rows();

        Matrix<FT> tmp(rows, clumns);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < clumns; ++j)
                tmp[i][j] = A[j][i];

        return tmp;
    }


/**
* matrix conjugate transpose
*/
    template<typename FT>
    Matrix<FT> conjugate_transpose(const Matrix<FT> &A) {
        int rows = A.cols();
        int clumns = A.rows();

        Matrix<FT> tmp(rows, clumns);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < clumns; ++j)
                tmp[i][j] = conj(A[j][i]);

        return tmp;
    }


/**
* matrix-matrix tranpose multiplication: A^T * B.
*/
    template<typename FT>
    Matrix<FT> transpose_mult(const Matrix<FT> &A1, const Matrix<FT> &A2) {
        assert(A1.rows() == A2.rows());

        int rows = A1.cols();
        int cols = A2.cols();
        int K = A1.rows();

        Matrix<FT> tmp(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                for (int k = 0; k < K; ++k)
                    tmp[i][j] += A1[k][i] * A2[k][j];

        return tmp;
    }


/**
* matrix-vector tranpose multiplication: A^T * b.
*/
    template<typename FT>
    std::vector<FT> transpose_mult(const Matrix<FT> &A, const std::vector<FT> &v) {
        assert(A.rows() == v.size());

        int rows = A.rows();
        int cols = A.cols();

        std::vector<FT> tmp(cols);
        for (int i = 0; i < cols; ++i)
            for (int j = 0; j < rows; ++j)
                tmp[i] += A[j][i] * v[j];

        return tmp;
    }


/**
* matrix-matrix tranpose multiplication: A * B^T.
*/
    template<typename FT>
    Matrix<FT> mult_transpose(const Matrix<FT> &A1, const Matrix<FT> &A2) {
        assert(A1.cols() == A2.cols());

        int rows = A1.rows();
        int cols = A2.rows();
        int K = A1.cols();

        Matrix<FT> tmp(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                for (int k = 0; k < K; ++k)
                    tmp[i][j] += A1[i][k] * A2[j][k];

        return tmp;
    }


/**
* vector-vector tranpose multiplication: a * b^T.
*/
    template<typename FT>
    Matrix<FT> mult_transpose(const std::vector<FT> &a, const std::vector<FT> &b) {
        int rows = a.size();
        int cols = b.size();

        Matrix<FT> tmp(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                tmp[i][j] = a[i] * b[j];

        return tmp;
    }


/**
* Generate the identity matrix.
*/
    template<typename FT>
    Matrix<FT> identity(int N, const FT &x) {
        Matrix<FT> tmp(N, N);
        for (int i = 0; i < N; ++i)
            tmp[i][i] = x;

        return tmp;
    }


/**
* Get the diagonal entries of matrix.
*/
    template<typename FT>
    std::vector<FT> diagonal(const Matrix<FT> &A) {
        int nColumn = A.rows();
        if (nColumn > A.cols())
            nColumn = A.cols();

        std::vector<FT> tmp(nColumn);
        for (int i = 0; i < nColumn; ++i)
            tmp[i] = A[i][i];

        return tmp;
    }


/**
* Generate the diagonal of matrix by given its diagonal elements.
*/
    template<typename FT>
    Matrix<FT> diagonal(const std::vector<FT> &d) {
        int N = static_cast<int>(d.size());

        Matrix<FT> tmp(N, N);
        for (int i = 0; i < N; ++i)
            tmp[i][i] = d[i];

        return tmp;
    }


/**
* the trace of this matrix, i.e. the sum of the coefficients on the main diagonal.
* NOTE: the matrix can be any matrix, not necessarily square.
*/
    template<typename FT>
    FT trace(const Matrix<FT> &A) {
        int range = std::min(A.rows(), A.cols());
        FT trac = 0.0;
        for (int i = 0; i < range; i++) {
            trac += A[i][i];
        }
        return trac;
    }


/**
* Compute Frobenius norm of matrix.
*/
    template<typename FT>
    FT norm(const Matrix<FT> &A) {
        int m = A.rows();
        int n = A.cols();

        FT sum = 0;
        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                sum += A[i][j] * A[i][j];

        return std::sqrt(sum);
    }


/**
* Swap two matrices.
*/
    template<typename FT>
    void swap(Matrix<FT> &lhs, Matrix<FT> &rhs) {
        int m = lhs.rows();
        int n = lhs.cols();

        assert(m == rhs.rows());
        assert(n == rhs.cols());

        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                std::swap(lhs(i, j), rhs(i, j));
    }


/**
* Matrix's column vectors sum.
*/
    template<typename FT>
    std::vector<FT> sum(const Matrix<FT> &A) {
        int m = A.rows();
        int n = A.cols();
        std::vector<FT> s(n);

        for (int j = 0; j < n; ++j)
            for (int i = 0; i < m; ++i)
                s[j] += A[i][j];
        return s;
    }


/**
* Minimum of matrix's column vectors.
*/
    template<typename FT>
    std::vector<FT> min(const Matrix<FT> &A) {
        int m = A.rows();
        int n = A.cols();
        std::vector<FT> sum(n);

        for (int j = 0; j < n; ++j) {
            FT tmp = A[0][j];
            for (int i = 1; i < m; ++i)
                if (tmp > A[i][j])
                    tmp = A[i][j];
            sum[j] = tmp;
        }

        return sum;
    }


/**
* Maximum of matrix's column vectors.
*/
    template<typename FT>
    std::vector<FT> max(const Matrix<FT> &A) {
        int m = A.rows();
        int n = A.cols();
        std::vector<FT> sum(n);

        for (int j = 0; j < n; ++j) {
            FT tmp = A[0][j];
            for (int i = 1; i < m; ++i)
                if (tmp < A[i][j])
                    tmp = A[i][j];
            sum[j] = tmp;
        }

        return sum;
    }


/**
* Matrix's column vectors mean.
*/
    template<typename FT>
    inline std::vector<FT> mean(const Matrix<FT> &A) {
        return sum(A) / FT(A.rows());
    }


// Get magnitude of a std::complex matrix.
    template<typename FT>
    Matrix<FT> abs(const Matrix<std::complex<FT> > &A) {
        int m = A.rows(),
                n = A.cols();
        Matrix<FT> tmp(m, n);

        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                tmp[i][j] = std::abs(A[i][j]);

        return tmp;
    }

// Get angle of a std::complex matrix.
    template<typename FT>
    Matrix<FT> arg(const Matrix<std::complex<FT> > &A) {
        int m = A.rows(),
                n = A.cols();
        Matrix<FT> tmp(m, n);

        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                tmp[i][j] = arg(A[i][j]);

        return tmp;
    }

// Get real part of a std::complex matrix.
    template<typename FT>
    Matrix<FT> real(const Matrix<std::complex<FT> > &A) {
        int m = A.rows(),
                n = A.cols();
        Matrix<FT> tmp(m, n);

        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                tmp[i][j] = A[i][j].real();

        return tmp;
    }


// Get imaginary part of a std::complex matrix.
    template<typename FT>
    Matrix<FT> imag(const Matrix<std::complex<FT> > &A) {
        int m = A.rows(),
                n = A.cols();
        Matrix<FT> tmp(m, n);

        for (int i = 0; i < m; ++i)
            for (int j = 0; j < n; ++j)
                tmp[i][j] = A[i][j].imag();

        return tmp;
    }


// Convert real matrix to std::complex matrix.
    template<typename FT>
    Matrix<std::complex<FT> > complex_matrix(const Matrix<FT> &rA) {
        int rows = rA.rows();
        int cols = rA.cols();

        Matrix<std::complex<FT> > cA(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                cA[i][j] = rA[i][j];

        return cA;
    }

// Convert real matrix to std::complex matrix (A for real, B for imaginary)
    template<typename FT>
    Matrix<std::complex<FT> > complex_matrix(const Matrix<FT> &mR, const Matrix<FT> &mI) {
        int rows = mR.rows();
        int cols = mR.cols();

        assert(rows == mI.rows());
        assert(cols == mI.cols());

        Matrix<std::complex<FT> > cA(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                cA[i][j] = std::complex<FT>(mR[i][j], mI[i][j]);

        return cA;
    }


//-----------------------------  std::vector  Implementation ------------------------------------//


// Overload the output stream function.
    template<typename FT>
    std::ostream &operator<<(std::ostream &out, const std::vector<FT> &A) {
        out << "size: " << A.size() << "\n";
        for (std::size_t i = 0; i < A.size(); ++i)
            out << A[i] << "\t";
        out << "\n";
        return out;
    }


// Overload the input stream function.
    template<typename FT>
    std::istream &operator>>(std::istream &in, std::vector<FT> &A) {
        int size;
        in >> size;
        A.resize(size);
        for (int i = 0; i < size; ++i)
            in >> A[i];
        return in;
    }


// get negative vector
    template<typename FT>
    std::vector<FT> operator-(const std::vector<FT> &A) {
        std::vector<FT> tmp = A;
        for (std::size_t i = 0; i < tmp.size(); ++i)
            tmp[i] = -tmp[i];

        return tmp;
    }


// vector-vector addition
    template<typename FT>
    std::vector<FT> operator+(const std::vector<FT> &A, const std::vector<FT> &B) {
        assert(A.size() == B.size());

        std::vector<FT> tmp = A;
        for (std::size_t i = 0; i < B.size(); ++i)
            tmp[i] += B[i];

        return tmp;
    }

// vector-scalar addition.
    template<typename FT>
    std::vector<FT> operator+(const std::vector<FT> &A, const FT &s) {
        std::vector<FT> tmp = A;
        for (std::size_t i = 0; i < A.size(); ++i)
            tmp[i] += s;

        return tmp;
    }

// scalar-vector addition.
    template<typename FT>
    std::vector<FT> operator+(const FT &s, const std::vector<FT> &A) {
        std::vector<FT> tmp = A;
        for (std::size_t i = 0; i < A.size(); ++i)
            tmp[i] += s;

        return tmp;
    }

// vector-vector subtraction
    template<typename FT>
    std::vector<FT> operator-(const std::vector<FT> &A, const std::vector<FT> &B) {
        assert(A.size() == B.size());

        std::vector<FT> tmp = A;
        for (std::size_t i = 0; i < B.size(); ++i)
            tmp[i] -= B[i];

        return tmp;
    }

// vector-scalar subtraction.
    template<typename FT>
    std::vector<FT> operator-(const std::vector<FT> &A, const FT &s) {
        std::vector<FT> tmp = A;
        for (std::size_t i = 0; i < A.size(); ++i)
            tmp[i] -= s;

        return tmp;
    }

// scalar-vector subtraction.
    template<typename FT>
    std::vector<FT> operator-(const FT &s, const std::vector<FT> &A) {
        std::vector<FT> tmp = A;
        for (std::size_t i = 0; i < A.size(); ++i)
            tmp[i] = s - tmp[i];

        return tmp;
    }

// vector-scalar multiplication
    template<typename FT1, typename FT2>
    std::vector<FT1> operator*(const std::vector<FT1> &A, const FT2 &s) {
        std::vector<FT1> tmp = A;
        for (std::size_t i = 0; i < A.size(); ++i)
            tmp[i] *= s;

        return tmp;
    }

// vector-scalar multiplication
    template<typename FT1, typename FT2>
    std::vector<FT1> operator*(const FT2 &s, const std::vector<FT1> &A) {
        return A * s;
    }

// vector-scalar division
    template<typename FT1, typename FT2>
    std::vector<FT1> operator/(const std::vector<FT1> &A, const FT2 &s) {
        std::vector<FT1> tmp = A;
        for (std::size_t i = 0; i < A.size(); ++i)
            tmp[i] /= s;

        return tmp;
    }

// vector-scalar division
    template<typename FT1, typename FT2>
    std::vector<FT1> operator/(const FT2 &s, const std::vector<FT1> &A) {
        std::vector<FT1> tmp = A;
        for (std::size_t i = 0; i < A.size(); ++i)
            tmp[i] = s / tmp[i];

        return tmp;
    }

// vector-vector element-wise multiplication
    template<typename FT>
    inline std::vector<FT> elem_mult(const std::vector<FT> &v1, const std::vector<FT> &v2) {
        assert(v1.size() == v2.size());
        std::vector<FT> result = v1;
        for (std::size_t i = 0; i < v2.size(); ++i) {
            result[i] *= v2[i];
        }

        return result;
    }

    template<typename FT>
    inline std::vector<FT> &elem_mult_eq(std::vector<FT> &v1, const std::vector<FT> &v2) {
        assert(v1.size() == v2.size());
        for (std::size_t i = 0; i < v2.size(); ++i) {
            v1[i] *= v2[i];
        }

        return v1;
    }


// vector-vector element-wise division
    template<typename FT>
    inline std::vector<FT> elem_divd(const std::vector<FT> &v1, const std::vector<FT> &v2) {
        assert(v1.size() == v2.size());
        std::vector<FT> result = v1;
        for (std::size_t i = 0; i < v2.size(); ++i) {
            result[i] /= v2[i];
        }

        return result;
    }


    template<typename FT>
    inline std::vector<FT> &elem_divd_eq(std::vector<FT> &v1, const std::vector<FT> &v2) {
        assert(v1.size() == v2.size());
        for (std::size_t i = 0; i < v2.size(); ++i) {
            v1[i] /= v2[i];
        }

        return v1;
    }


// Inner product for vectors.
    template<typename FT>
    FT operator*(const std::vector<FT> &v1, const std::vector<FT> &v2) {
        assert(v1.size() == v2.size());

        FT sum = 0;
        typename std::vector<FT>::const_iterator itr1 = v1.begin();
        typename std::vector<FT>::const_iterator itr2 = v2.begin();

        while (itr1 != v1.end())
            sum += (*itr1++) * (*itr2++);

        return sum;
    }

// Inner product for vectors.
    template<typename FT>
    FT dot(const std::vector<FT> &v1, const std::vector<FT> &v2) {
        assert(v1.size() == v2.size());

        FT sum = 0;
        typename std::vector<FT>::const_iterator itr1 = v1.begin();
        typename std::vector<FT>::const_iterator itr2 = v2.begin();

        while (itr1 != v1.end())
            sum += (*itr1++) * (*itr2++);

        return sum;
    }

/// custom
// Outer product for vectors.
// only works for vectors of length 3
    template<typename FT>
    std::vector<FT> cross(const std::vector<FT> &v1, const std::vector<FT> &v2) {
        assert(v1.size() == v2.size());
        assert(v1.size() == 3);

        std::vector<FT> res(v1.size());

        res[0] = (v1[1] * v2[2]) - (v1[2] * v2[1]);
        res[1] = (v1[2] * v2[0]) - (v1[0] * v2[2]);
        res[2] = (v1[0] * v2[1]) - (v1[1] * v2[0]);

        return res;
    }


// std::vector's sum.
    template<typename FT>
    FT sum(const std::vector<FT> &v) {
        FT sum = 0;
        typename std::vector<FT>::const_iterator itr = v.begin();

        while (itr != v.end())
            sum += *itr++;

        return sum;
    }

// Minimum value of vector.
    template<typename FT>
    FT min(const std::vector<FT> &v) {
        FT m = v[0];
        for (std::size_t i = 1; i < v.size(); ++i)
            if (m > v[i])
                m = v[i];

        return m;
    }

// Maximum value of vector.
    template<typename FT>
    FT max(const std::vector<FT> &v) {
        FT M = v[0];
        for (std::size_t i = 1; i < v.size(); ++i)
            if (M < v[i])
                M = v[i];

        return M;
    }

// std::vector's mean.
    template<typename FT>
    FT mean(const std::vector<FT> &v) {
        return sum(v) / FT(v.size());
    }

// std::vector's norm in Euclidean space.
    template<typename FT>
    FT norm(const std::vector<FT> &v) {
        FT sum = 0;
        typename std::vector<FT>::const_iterator itr = v.begin();

        while (itr != v.end()) {
            sum += (*itr) * (*itr);
            itr++;
        }

        return FT(std::sqrt(sum));
    }

// std::vector's norm in Euclidean space.
    template<typename FT>
    FT norm(const std::vector<std::complex<FT> > &v) {
        FT sum = 0;
        typename std::vector<std::complex<FT> >::const_iterator itr = v.begin();

        while (itr != v.end()) {
            sum += std::norm(*itr++);
        }

        return FT(std::sqrt(sum));
    }

// Swap two vectors.
    template<typename FT>
    void swap(std::vector<FT> &lhs, std::vector<FT> &rhs) {
        typename std::vector<FT>::iterator itrL = lhs.begin(),
                itrR = rhs.begin();

        while (itrL != lhs.end())
            std::swap(*itrL++, *itrR++);
    }

// Generates a vector of n points linearly spaced between and including a and b.
    template<typename FT>
    std::vector<FT> linspace(FT a, FT b, int n) {
        if (n < 1)
            return std::vector<FT>();
        else if (n == 1)
            return std::vector<FT>(1, a);
        else {
            FT dx = (b - a) / (n - 1);

            std::vector<FT> tmp(n);
            for (int i = 0; i < n; ++i)
                tmp[i] = a + i * dx;

            return tmp;
        }
    }

// Get magnitude of a complex vector.
    template<typename FT>
    std::vector<FT> abs(const std::vector<std::complex<FT> > &v) {
        std::vector<FT> tmp(v.size());
        typename std::vector<FT>::iterator itrL = tmp.begin();
        typename std::vector<std::complex<FT> >::const_iterator itrR = v.begin();

        while (itrL != tmp.end())
            *itrL++ = std::abs(*itrR++);

        return tmp;
    }

// Get angle of a complex vector.
    template<typename FT>
    std::vector<FT> arg(const std::vector<std::complex<FT> > &v) {
        std::vector<FT> tmp(v.size());
        typename std::vector<FT>::iterator itrL = tmp.begin();
        typename std::vector<std::complex<FT> >::const_iterator itrR = v.begin();

        while (itrL != tmp.end())
            *itrL++ = std::arg(*itrR++);

        return tmp;
    }

// Get real part of a complex vector.
    template<typename FT>
    std::vector<FT> real(const std::vector<std::complex<FT> > &v) {
        std::vector<FT> tmp(v.size());
        typename std::vector<FT>::iterator itrL = tmp.begin();
        typename std::vector<std::complex<FT> >::const_iterator itrR = v.begin();

        while (itrL != tmp.end())
            *itrL++ = (*itrR++).real();

        return tmp;
    }

// Get imaginary part of a complex vector.
    template<typename FT>
    std::vector<FT> imag(const std::vector<std::complex<FT> > &v) {
        std::vector<FT> tmp(v.size());
        typename std::vector<FT>::iterator itrL = tmp.begin();
        typename std::vector<std::complex<FT> >::const_iterator itrR = v.begin();

        while (itrL != tmp.end())
            *itrL++ = (*itrR++).imag();

        return tmp;
    }

// Convert real vector to complex vector.
    template<typename FT>
    std::vector<std::complex<FT> > complex_vector(const std::vector<FT> &rv) {
        int N = rv.size();

        std::vector<std::complex<FT> > cv(N);
        typename std::vector<std::complex<FT> >::iterator itrL = cv.begin();
        typename std::vector<FT>::const_iterator itrR = rv.begin();

        while (itrR != rv.end())
            *itrL++ = *itrR++;

        return cv;
    }

// Convert real vector to complex vector (A for real, B for imaginary)
    template<typename FT>
    std::vector<std::complex<FT> > complex_vector(const std::vector<FT> &vR, const std::vector<FT> &vI) {
        int N = static_cast<int>(vR.size());

        assert(N == vI.size());

        std::vector<std::complex<FT> > cv(N);
        typename std::vector<std::complex<FT> >::iterator itrC = cv.begin();
        typename std::vector<FT>::const_iterator itrR = vR.begin(),
                itrI = vI.begin();

        while (itrC != cv.end())
            *itrC++ = std::complex<FT>(*itrR++, *itrI++);

        return cv;
    }
}

#endif // _MATH_MATRIX_H_
