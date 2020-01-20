#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cstddef>
#include <vector>
#include <ostream>
#include <iostream>

namespace dauphine
{
    class matrix
    {
    public:

        matrix(std::size_t nb_rows, std::size_t nb_cols);

    public:

        std::size_t m_nb_rows;
        std::size_t m_nb_cols;
        std::vector<double> m_data;
        std::size_t nb_rows() const;
        std::size_t nb_cols() const;
        std::vector<double> produit_mat_vect(std::vector<double>& v);
        
        
        void resize(std::size_t nb_rows,std::size_t nb_cols);
        
        double& operator()(std::size_t i, std::size_t j);
        matrix& operator+=(const matrix& rhs);
    };
    std::ostream& operator<<(std::ostream& out, matrix& m);
    matrix operator+(const matrix& lhs, const matrix& rhs);
    matrix operator+(double lhs, matrix& m);
    

}


#endif
