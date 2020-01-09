//
//  pde_solver.hpp
//  PDE_solver
//
//  Created by Florian on 08/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#ifndef pde_solver_hpp
#define pde_solver_hpp

#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>


namespace Solve
{

    class matrix_pde_case1
    {
    public:
        matrix_pde_case1(std::size_t maturity, double time_path, double sigma, double rate, double theta,  std::size_t number_path_mash, double S_max, double S_min);
    
    public: //private at the end of the project
        std::size_t m_underlying_nb_dimension;
        std::size_t m_maturity_nb_dimension;
        std::size_t number_path_underlying();
        std::size_t m_maturity();
        
        double underlying_max;
        double underling_min;
        double path_time;
        double time_product;
        
        std::vector<double> forward_coefficient(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min, double temp);
        std::vector<double> present_coefficient(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min, double temp);
        std::vector<double> backward_coefficient(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,  double temp);
        std::vector<double> matrix_lhs(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,  double theta);
        std::vector<double> matrix_rhs(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,  double theta);
        std::vector<double> LU_decomposition(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,  double theta);
        std::vector<double> Final_solve(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,  double theta);
        
        
    };
}
#endif /* pde_solver_hpp */
