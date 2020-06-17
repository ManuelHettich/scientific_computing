//=============================================================================
//                                                                            
//   Exercise code for the lecture
//   "Scientific Computing"
//   by Prof. Dr. Mario Botsch, Bielefeld University
//
//   Copyright (C) 2018  Computer Graphics Group, Bielefeld University.
//                                                                            
//=============================================================================
#pragma once
//=============================================================================

#include <math.h>
#include <iostream>
#include <assert.h>

#if defined(__SSE__)
#  include <xmmintrin.h>
#endif

#include <cstring>
#include "types.h"

//=============================================================================

namespace utils {

//== CLASS DEFINITION =========================================================

/// vector class
class Vector
{
public:

    /// default constructor
    Vector(unsigned int _n=0, const Scalar _s=0.0) : rows_(0), data_(0)
    {
        if (_n)
        {
            resize(_n);
            fill(_s);
        }
    }


    /// copy constructor
    Vector(const Vector& _v) : rows_(0), data_(0)
    {
        *this = _v;
    }


    /// destructor
    ~Vector()
    {
        deallocate();
    }


    /// cast to Scalar array
    operator const Scalar*() const { return data_; }
    operator       Scalar*()       { return data_; }


    /// get Scalar array
    const Scalar* data() const { return data_; }
    Scalar*       data()       { return data_; }


    /// resize vector
    void resize(unsigned int _n) 
    { 
        assert(_n);
        if (rows_ != _n)
        {
            deallocate();
            allocate(_n);
            rows_ = _n;
        }
    }    


    /// fill matrix with scalar _s
    void fill(const Scalar _s)
    {
        for (unsigned int i=0; i<rows_; ++i)
            (*this)(i) = _s;
    }


    /// assignment operator
    Vector& operator=(const Vector& _v)
    {
        resize(_v.size());
        memcpy(data_, _v.data_, rows_*sizeof(Scalar));
        return *this;
    }


    /// read & write element access
    Scalar& operator()(const unsigned int _i)
    {
        assert(_i < rows_);
        return data_[_i];
    }  


    /// read only element access
    const Scalar& operator()(const unsigned int _i) const
    {
        assert(_i < rows_);
        return data_[_i];
    }


    /// return vector's dimension
    unsigned int rows() const { return rows_; }
    unsigned int size() const { return rows_; }


    /// negate: -vector
    const Vector operator-()
    {
        Vector v(rows_);
        for (unsigned int i=0; i<rows_; ++i)
            v(i) = -(*this)(i);
        return v;
    }


    /// vector *= scalar
    Vector& operator*=(Scalar _s)
    {
        for (unsigned int i=0; i<rows_; ++i)
            (*this)(i) *= _s;
        return *this;
    }

    /// vector * scalar
    const Vector operator*(const Scalar _s) const
    {
        return Vector(*this) *= _s;
    }

    /// scalar * vector
    friend const Vector operator*(const Scalar _s, const Vector& _v)
    {
        return Vector(_v) *= _s;
    }


    /// vector += vector
    Vector& operator+=(const Vector& _v)
    {
        assert(size() == _v.size());
        for (unsigned int i=0; i<rows_; ++i)
            (*this)(i) += _v(i);
        return *this;
    }

    /// vector + vector
    const Vector operator+(const Vector& _v) const
    {
        assert(size() == _v.size());
        return Vector(*this) += _v;
    }


    /// vector -= vector
    Vector& operator-=(const Vector& _v)
    {
        assert(size() == _v.size());
        for (unsigned int i=0; i<rows_; ++i)
            (*this)(i) -= _v(i);
        return *this;
    }

    /// vector - vector
    const Vector operator-(const Vector& _v) const
    {
        assert(size() == _v.size());
        return Vector(*this) -= _v;
    }


    /// scalar product: vector * vector
    const Scalar operator*(const Vector& _v) const
    {
        assert(size() == _v.size());
        Scalar s(0.0);
        for (unsigned int i=0; i<rows_; ++i)
            s += (*this)(i) * _v(i);
        return s;
    }

    /// euclidean 2-norm
    const Scalar norm() const
    {
        return sqrt( (*this) * (*this) );
    }


private:

    // memory management (use SSE if it is enabled by -msse)
#if defined(__SSE__)
    void deallocate() {    if (data_) _mm_free(data_);    }
    void allocate(unsigned int _n) { data_ = (Scalar*)_mm_malloc(_n*sizeof(Scalar), 16); }
#else
    void deallocate() {    delete[] data_;    }
    void allocate(unsigned int _n) {    data_ = new Scalar[_n];    }
#endif


private:

    unsigned int  rows_;
    Scalar*       data_;
};


//== CLASS DEFINITION =========================================================


/// matrix class
class Matrix
{
public:

    /// construct with number of rows and columns.
    Matrix(unsigned int _rows=0, unsigned int _cols=0, Scalar _s=0.0)
        : rows_(0), cols_(0), data_(0)
    {
        if (_rows && _cols)
        {
            resize(_rows, _cols);
            fill(_s);
        }
    }


    /// copy constructor
    Matrix(const Matrix& _M)
    : rows_(0), cols_(0), data_(0)
    {
        *this = _M;
    }


    /// destructor
    ~Matrix()
    {
        deallocate();
    }


    /// cast to Scalar array
    operator const Scalar*() const { return data_; }
    operator       Scalar*()       { return data_; }


    /// get Scalar array
    const Scalar* data() const { return data_; }
    Scalar*       data()       { return data_; }


    /// assignment operator
    Matrix& operator=(const Matrix& _M)
    {
        resize(_M.rows(), _M.cols());
        memcpy(data_, _M.data_, rows_*cols_*sizeof(Scalar));
        return *this;
    }


    /// resize matrix
    void resize(unsigned int _rows, unsigned int _cols) 
    { 
        assert(_rows && _cols);
        if (rows_!=_rows || cols_!=_cols)
        {
            deallocate();
            allocate(_rows*_cols);
            rows_ = _rows; 
            cols_ = _cols; 
        }
    }


    /// fill matrix with scalar _s
    void fill(Scalar _s)
    { 
        for (unsigned int i=0; i<rows_; ++i)
            for (unsigned int j=0; j<cols_; ++j)
                (*this)(i,j) = _s;
    }


    /// set matrix to be the identity matrix
    void identity()
    { 
        fill(0.0);
        for (unsigned int i=0; i<std::min(rows_,cols_); ++i)
            (*this)(i,i) = 1.0;
    }    

    /// read & write access to matrix elements
    Scalar& operator()(const unsigned int _i, const unsigned int _j)
    {
        assert (_i < rows_ && _j < cols_);
        return data_[_i + _j * rows_];
    }

    /// read only access to matrix elements
    const Scalar& operator()(const unsigned int _i, const unsigned int _j) const
    {
        assert (_i < rows_ && _j < cols_);
        return data_[_i + _j * rows_];
    }


    /// number of rows
    unsigned int rows() const { return rows_; }

    /// number of columns
    unsigned int cols() const { return cols_; }


    /// negate: -matrix
    const Matrix operator-()
    {
        Matrix M(rows_, cols_);
        for (unsigned int i=0; i<rows_; ++i)
            for (unsigned int j=0; j<cols_; ++j)
                M(i,j) = -(*this)(i,j);
        return M;
    }


    /// matrix *= scalar
    Matrix& operator*=(const Scalar _s)
    {
        for (unsigned int i=0; i<rows_; ++i)
            for (unsigned int j=0; j<cols_; ++j)
                (*this)(i,j) *= _s;
        return *this;
    }

    /// matrix * scalar
    const Matrix operator*(const Scalar _s) const
    {
        return Matrix(*this) *= _s;
    }

    /// scalar * matrix
    friend const Matrix operator*(const Scalar _s, const Matrix& _M)
    {
        return Matrix(_M) *= _s;
    }


    /// matrix += matrix
    Matrix& operator+=(const Matrix& _M)
    {
        assert(rows()==_M.rows() && cols()==_M.cols());
        for (unsigned int i=0; i<rows_; ++i)
            for (unsigned int j=0; j<cols_; ++j)
                (*this)(i,j) += _M(i,j);
        return *this;
    }    

    /// matrix + matrix
    const Matrix operator+(const Matrix& _M) const
    {
        assert(rows()==_M.rows() && cols()==_M.cols());
        return Matrix(*this) += _M;
    }


    /// matrix -= matrix
    Matrix& operator-=(const Matrix& _M)
    {
        assert(rows()==_M.rows() && cols()==_M.cols());
        for (unsigned int i=0; i<rows_; ++i)
            for (unsigned int j=0; j<cols_; ++j)
                (*this)(i,j) -= _M(i,j);
        return *this;
    }    

    /// matrix - matrix
    const Matrix operator-(const Matrix& _M) const
    {
        assert(rows()==_M.rows() && cols()==_M.cols());
        return Matrix(*this) -= _M;
    }


    /// matrix * matrix
    const Matrix operator*(const Matrix& _M) const
    {
        assert(cols()==_M.rows());
        const unsigned int l=rows(), m=cols(), n=_M.cols();
        Matrix M(l,n,0.0);
        for (unsigned int i=0; i<l; ++i)
            for (unsigned int j=0; j<n; ++j)
                for (unsigned int k=0; k<m; ++k)
                    M(i,j) += (*this)(i,k) * _M(k,j);
        return M;
    }

    /// matrix *= matrix
    Matrix& operator*=(const Matrix& _M)
    {
        assert(cols()==_M.rows());
        return (*this = *this * _M);
    }


    /// matrix * vector
    const Vector operator*(const Vector& _x) const
    {
        Vector v(rows_);
        for (unsigned int i=0; i<rows_; ++i)
        {
            v(i) = 0.0;
            for (unsigned int j=0; j<cols_; ++j)
                v(i) += (*this)(i,j) * _x(j);
        }
        return v;
    }


    /// frobenius norm
    const Scalar norm() const
    {
        Scalar s(0.0);
        for (unsigned int i=0; i<rows_; ++i)
            for (unsigned int j=0; j<cols_; ++j)
                s += (*this)(i,j) * (*this)(i,j);
        return sqrt(s);
    }


    /// matrix transpose
    const Matrix transpose() const
    {
        Matrix T(cols_, rows_);
        for (unsigned int i=0; i<rows_; ++i)
            for (unsigned int j=0; j<cols_; ++j)
                T(j,i) = (*this)(i,j);
        return T;
    }


private:

    // memory management (use SSE if it is enabled by -msse)
#if defined(__SSE__)
    void deallocate() {    if (data_) _mm_free(data_);    }
    void allocate(unsigned int _n) { data_ = (Scalar*)_mm_malloc(_n*sizeof(Scalar), 16); }
#else
    void deallocate() {    delete[] data_;    }
    void allocate(unsigned int _n) {    data_ = new Scalar[_n];    }
#endif


private:

    unsigned int  rows_, cols_;
    Scalar*       data_;
};


//=============================================================================


/// output a Vector to a stream
inline std::ostream& 
operator<<(std::ostream& _os, const Vector& _v)
{
    for (unsigned int i=0; i<_v.size(); ++i)
        _os << _v(i) << ' ';
    _os << std::endl;
    return _os;
}


/// output a Matrix to a stream
inline std::ostream& 
operator<<(std::ostream& _os, const Matrix& _M)
{
    for (unsigned int i=0; i<_M.rows(); ++i)
    {
        for(unsigned int j=0; j<_M.cols(); ++j)
            _os << _M(i,j) << ' ';
            _os << std::endl;
    }
    return _os;
}


//=============================================================================
} // namespace
//=============================================================================
