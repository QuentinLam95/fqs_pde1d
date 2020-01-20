//
//  main.cpp
//  Project C++
//
//  Created by Florian on 02/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#include <iostream>
#include <string>
#include <vector>
#include "Rates.hpp"
#include "Payoff.h"
#include "pde.hpp"
#include "pde_solver.hpp"
#include "matrix.hpp"

/*

namespace dauphine
{

    void Dirichlet_boundaries_conditions()
    {
        Boundaries_Dirichlet BD(0,200);

        double Payoff_0=BD.get_initial_cond();
        double Payoff_N=BD.get_final_cond();
        
        std::cout << "The lower Dirichlet condition is " <<Payoff_0 <<std::endl;
        std::cout << "The upper Dirichlet condition is " <<Payoff_N <<std::endl;
        
    };

    void Neumann_boundaries_conditions()
    {
        Boundaries_Neumann BN(0,100);

        double Payoff_0_df=BN.get_initial_cond();
        double Payoff_N_df=BN.get_final_cond();
    
        std::cout << "The lower Neumann condition is " << Payoff_0_df <<std::endl;
        std::cout << "The upper Neumann condition is " <<Payoff_N_df <<std::endl;
        
}

    }

void Rates_diffusion(std::string Model)
{
    if (Model=="Vasicek model"){
        std::vector<double> V=Vasicek_diffusion(1,365,1000,0.2,0.1,0.012,0.01);
        print(V);
    }
    
    else if(Model=="Cox Ingersoll Ross model"){
        std::vector<double> C=Cox_Ingersoll_Ross_diffusion(1,365,1000,0.2,0.1,0.012,0.01);
        print(C);
    }
    
    else if(Model=="Rendleman Batter model"){
        std::vector<double> RB=Rendleman_Batter_diffusion(1,365,1000,0.1,0.012,0.01);
        print(RB);
    }
    else{
        std::cout<<"Please select a rate model ";
    }
    
};

void test(int &i)   // i est une référence du paramètre constant.
{
    i = 2;    // Modifie le paramètre passé en référence.
    return;
}
*/

int main(int argc, char* argv[])
{
    // insert code here...
    //dauphine::Dirichlet_boundaries_conditions();
    //dauphine::Neumann_boundaries_conditions();
    //Rates_diffusion("Vasicek model");
	
	double S0 = 100.0;
	double K = 100.0;
	double sigma = 0.15;
	double maturity = 1;
	double r = 0.01;
	double theta = 0.5;
	std::size_t space_dim = 100;
	std::size_t time_dim = 50;
	std::string l_boundary_type = "D";
	std::string r_boundary_type = "D";
	PayOff* pay_off_call = new PayOffCall(K);
	VanillaOption* call_option = new VanillaOption(K, r, maturity, sigma, pay_off_call);
	BS_PDE* bs_pde = new BS_PDE(call_option, l_boundary_type, r_boundary_type);
	Solve::matrix_pde_case1* PDE_solve = new Solve::matrix_pde_case1(bs_pde, theta, space_dim, time_dim, S0, maturity);
	PDE_solve->Crout_Algo_Resolution();
	double price = PDE_solve->get_price(S0);


	std::cout << price << std::endl;


    //Discretization Dis(S_max, number_path_mash,S_min,maturity,time_path);
    //std::vector<double> V=Dis.get_diffusion_time();
    //Solve::matrix_pde_case1 M(10, 10, 0.1,0.05, 0.5, 5, 200, 0);
    //std::vector<double> V=M.Final_solve(10, 10, 0.1, 0.5, 5, 200, 0,0.5);
    //dauphine::matrix m(3,3);
    //std::vector<double> v(3);
    //std::vector<double> result=m.produit_mat_vect(m,v);
    //print(V);
    
    return 0;
    
}
