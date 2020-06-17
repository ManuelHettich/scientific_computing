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

#include "types.h"
#include "sc_math.h"
#include <vector>

//=============================================================================

namespace utils {

//== CLASS DEFINITION =========================================================

/// Simple SparseMatrix class
class SparseMatrix
{
public:

    /// constructor
    SparseMatrix(int _nrows=0, int _ncols=0) : n_rows_(_nrows), n_cols_(_ncols) {}


    /// clear the matrix (if you have to fill it with new values)
    void clear()
    {
        values_.clear();
        rowptr_.clear();
        colind_.clear();
    }


    /// in case you have used the default constructor and want to resize the matrix
    void resize(int _nrows, int _ncols)
    {
        n_cols_ = _ncols;
        n_rows_ = _nrows;
        clear();
    }


    /// begin the next row
    void begin_row() 
    {
        if (rowptr_.empty() || rowptr_.back() != (int)values_.size())
            rowptr_.push_back(values_.size());
        
        assert(rowptr_.size() <= n_rows_ && "too many rows");
    }


    /// set value of column _j in current row to _value
    void add_value(unsigned int _j, Scalar _value)    
    {
        assert(_j >= 0 && _j < n_cols_);
        
        values_.push_back(_value);
        colind_.push_back(_j);
    }


    /// finish current row
    void end_row()
    {
        if (rowptr_.empty() || rowptr_.back() != (int)values_.size())
            rowptr_.push_back(values_.size());
        
        assert(rowptr_.size()-1 <= n_rows_ && "too many rows");
    }


    /// multiply by vector: y=A*_x
    const utils::Vector operator*(const utils::Vector& _x) const
    {
        assert(_x.size() == n_cols_);
        
        utils::Vector y(n_rows_, 0.0);
        
        for (size_t i=0; i<n_rows_; ++i)
            for (int p=rowptr_[i]; p<rowptr_[i+1]; ++p)
                y(i) += values_[p] * _x(colind_[p]);
        
        return y;
    }


    /// number of rows, columns, and non-zero elements
    const int rows()     const { return n_rows_; }
    const int cols()     const { return n_cols_; }
    const int nonzeros() const { return values_.size(); }

    void print() const
    {
        std::cout << "values: ";
        for (size_t i=0; i<values_.size(); ++i)
            std::cout << values_[i] << ", ";
        std::cout << std::endl;
        
        std::cout << "rowptr: ";
        for (size_t i=0; i<rowptr_.size(); ++i)
            std::cout << rowptr_[i] << ", ";
        std::cout << std::endl;

        std::cout << "colind: ";
        for (size_t i=0; i<colind_.size(); ++i)
            std::cout << colind_[i] << ", ";
        std::cout << std::endl;
    }


private:

    // SparseCholeskySolver has to get low-level access
    friend class SparseCholeskySolver;

    // ParallelCGSovler has to get low-level access
    friend class ParallelCGSolver;


    std::vector<Scalar>  values_;
    std::vector<int>     rowptr_;
    std::vector<int>     colind_;    
    unsigned int         n_rows_;
    unsigned int         n_cols_;
};


//=============================================================================
} // namespace
//=============================================================================

