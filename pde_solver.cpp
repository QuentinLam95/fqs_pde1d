//
//  pde_solver.cpp
//  PDE_solver
//
//  Created by Florian on 08/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#include "pde_solver.hpp"

namespace Solve
{
    matrix_pde_case1::matrix_pde_case1(BS_PDE* _pde, double _theta, std::size_t _space_dim, std::size_t _time_dim, double _S0, double _maturity)
    :pde(_pde), theta(_theta), space_dim(_space_dim), time_dim(_time_dim), S0(_S0), maturity(_maturity) {
		calculate_parameters();
	}
 
	void matrix_pde_case1::calculate_parameters() {
		double stdv = pde->standard_dev();

		x_max = log(S0) + 5 * stdv;
		x_min = log(S0) - 5 * stdv;
		dx = 10 * stdv / static_cast<double>(space_dim);
		dt =  maturity / static_cast<double>(time_dim);
		r = pde->get_right_boundary_type();
		l = pde->get_left_boundary_type();
	}
    
	/*
	std::size_t  matrix_pde_case1::number_path_underlying()

	{
		return m_underlying_nb_dimension;
	}

	std::size_t  matrix_pde_case1::m_maturity()

	   {
		   return m_maturity_nb_dimension;
	   }
	*/

        
    //if (side==true) A mettre dans un autre code
    //{
       // double temp=theta*maturity/time_path;
    //}
    //else
    //{
      //  double temp=-(1-theta)*maturity/time_path;
    //}
    
    std::vector<double> matrix_pde_case1::forward_coefficient(const double& temp)
    {
        double dx_2=pow(dx,2.0);
        std::vector<double>a_coef(space_dim-2);
        
        for(auto it = a_coef.begin(); it != a_coef.end(); ++it)
        {
			*it = temp *(-pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);             
        }
		return a_coef;
	}
    
    std::vector<double> matrix_pde_case1::present_coefficient(const double& temp)
    {
        double dx_2=pow(dx,2.0);
        std::vector<double>b_coef(space_dim-1);
		auto ptr = b_coef.begin();
		auto ftr = b_coef.end();

		if (l.compare("N") == 0) 
		{
			//present coeff + backward coeff because of neumann condition
			*ptr = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2) + temp * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);
		}
		else
		{
			*ptr = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2);
		}
		++ptr;

		for (auto iter = next(ptr,1); iter != b_coef.back(); ++iter)
		{
			*iter = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2);
			++iter;
		}

		if (r.compare("N") == 0)
		{
			//present coeff + forward coeff because of neumann condition
			*ftr = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2) + temp * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);
		}
		else
		{
			*ftr = 1 + temp * (pde->zero_coeff() + 2 * pde->diff_coeff() / dx_2);
		}

		return b_coef;
	}
    
    std::vector<double> matrix_pde_case1::backward_coefficient(const double& temp)
    {
		double dx_2 = pow(dx, 2.0);
		std::vector<double> c_coef(space_dim - 2);

		for (auto it = c_coef.begin(); it != c_coef.end(); ++it)
		{
			*it = temp * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2);
		}
		return c_coef;
    }
    
	void matrix_pde_case1::set_initial_conditions()
	{
		x_values.resize(space_dim - 1, 0.0);
		new_result.resize(space_dim - 1, 0.0);
		//filling x_values ; and initial condition in result vector
		for (auto x = x_values.begin(), auto r = new_result.begin(), double val = x_min + dx; x != x_values.end(); ++x, ++r, val += dx)
		{
			*x = val;
			*r = pde->init_cond(exp(val));
		}
	}

	std::vector<double> matrix_pde_case1::boundary_increment(const double& t)
	{
		std::vector<double> boundary_effect;
		boundary_effect.resize(space_dim -1, 0.0);
		auto it1 = boundary_effect.begin();
		auto it2 = boundary_effect.back();
		double left_b = pde->boundary_left(t);
		double right_b = pde->boundary_right(t, exp(x_max));

		if (l.compare("D") == 0)
		{
			*it1 = - left_b * dt * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2); //Boudary times the backward coefficient.
		}
		else
		{
			*it1 = left_b * dx * dt * (pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2); //Boudary times the backward coefficient.
		}

		if (r.compare("D") == 0)
		{
			*it2 = -right_b * dt * (-pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2); //Boudary times the forward coefficient.
		}
		else
		{
			*it2 = right_b * dx * dt * (-pde->conv_coeff() / (2 * dx) - pde->diff_coeff() / dx_2); //Boudary times the forward coefficient.
		}
	}
    
	dauphine::matrix matrix_pde_case1::transition_matrix(const double& temp)
	{
		std::vector<double> a = forward_coefficient(temp);
		std::vector<double> b = present_coefficient(temp);
		std::vector<double> c = backward_coefficient(temp);
		dauphine::matrix m(space_dim - 1, space_dim - 1);

		m(0, 0) = b[0];
		m(0, 1) = a[0];

		for (std::size_t i = 1; i != space_dim - 2 ; ++i)
		{
			m(i, i - 1) = c[i-1];
			m(i, i) = b[i];
			m(i, i + 1) = a[i];
		}

		m(space_dim - 2, space_dim - 3) = c[space_dim - 3];
		m(space_dim - 2, space_dim - 2) = b[space_dim - 2];

		return m ;
	}

	void matrix_pde_case1::Crout_Algo_Resolution()
	{
		double temp_lhs = theta * dt;
		double temp_rhs = - (1- theta) * dt;

		//creation of the left and right transition matrices
		dauphine::matrix M_lhs = transition_matrix(temp_lhs);
		dauphine::matrix M_rhs = transition_matrix(temp_rhs);

		// L - U decomposition of the left matrix for the Crout Algorithm 
		dauphine::matrix L(space_dim -1, space_dim - 1);
		dauphine::matrix U(space_dim - 1, space_dim - 1);
		L(0, 0) = M_lhs(0, 0);
		U(0, 0) = 1.0;
		U(0, 1) = M_lhs(0, 1) / L(0, 0);

		for (std::size_t i = 1; i != space_dim - 2; ++i)
		{
			L(i, i - 1) = M_lhs(i, i - 1);
			L(i, i) = M_lhs(i, i) - L(i, i - 1) * U(i - 1, i);
			U(i, i + 1) = M_lhs(i, i + 1) / L(i, i);
			U(i, i) = 1.0;
		}

		U(space_dim - 2, space_dim - 2) = 1.0;
		L(space_dim - 2, space_dim - 3) = M_lhs(space_dim - 2, space_dim - 3);
		L(space_dim - 2, space_dim - 2) = M_lhs(space_dim - 2, space_dim - 2) - L(space_dim - 2, space_dim - 3) * U(space_dim - 3, space_dim - 2);
		
		std::vector<double> v(space_dim - 1);
		std::vector<double> tmp(space_dim - 1);
		double t;

		for (std::size_t i = 1; i != time_dim; ++i)
		{
			old_result = new_result;
			t = maturity - i * dt; // think about && rvalue lvalue;
			tmp = boundary_increment(t);
			v = M_rhs.produit_mat_vect(M_rhs, old_result);
			std::transform(tmp.begin(), tmp.end(), v.begin(), std::plus<double>());

			new_result = LU_compute(L, U, v);
		}

	}

	std::vector<double> matrix_pde_case1::LU_compute(const dauphine::matrix& L, const dauphine::matrix& U, const std::vector<double>& b)
	{
		std::vector<double> x(b.size(),0.0);
		std::vector<double> z(b.size(),0.0);
		
		z[0] = b[0] / L(0, 0);

		for (std::size_t i = 1; i != b.size(); ++i)
		{
			z[i] = (b[i] - L(i, i - 1) * z[i - 1]) / L(i, i);
		}

		x[b.size() - 1] = z[b.size() - 1];
		for (auto i = b.size() - 1; i != 0; --i)
		{
			x[i] = z[i] - U(i, i + 1) * x[i + 1];
		}

		return x;
	}

	std::vector<double> matrix_pde_case1::get_price_curve()
	{
		//add test if resolution of the PDE has been done
		return new_result;
	}

	double matrix_pde_case1::get_price(const double& S)
	{
		//add test if resolution of the PDE has been done
		double x = log(S);
		int i = 0;
		while (x_values[i] < x)
		{
			++i
		}

		return new_result[i];
	}


        
    

}



