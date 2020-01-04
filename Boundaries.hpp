//
//  Boundaries.hpp
//  Project C++
//
//  Created by Florian on 02/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#ifndef Boundaries_hpp
#define Boundaries_hpp
#include <stdio.h>
#include <string>
#include <cmath>
#include <algorithm>


double Payoff (double S);
double Payoff_df (double S);
namespace dauphine
{
    class Boundaries_Dirichlet
    {
    public:
        
        explicit Boundaries_Dirichlet(double initial_cond,double final_cond);
        double get_initial_cond();
        double get_final_cond();
        ~Boundaries_Dirichlet();

        
    private:
        double BD_initial_cond;
        double BD_final_cond;

    };

    class Boundaries_Neumann
    {
    public:
        explicit Boundaries_Neumann(double initial_cond,double final_cond);
        double get_initial_cond();
        double get_final_cond();
        ~Boundaries_Neumann();
        
    private:
        double BN_initial_cond;
        double BN_final_cond;
    };
}

#endif /* Boundaries_hpp */
