//
//  Boundaries.cpp
//  Project C++
//
//  Created by Florian on 02/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#include "Boundaries.hpp"
#include <iostream>

double Payoff (double S)//Il faudra implémenter la fonction payoff
{
    return S;
};

double Payoff_df (double S)//Il faudra implémenter la fonction payoff
{
    return S;
};


namespace dauphine
{
    Boundaries_Dirichlet::Boundaries_Dirichlet(double initial_cond,double final_cond)
    :BD_initial_cond(initial_cond),
    BD_final_cond(final_cond)
    {
    };

    double Boundaries_Dirichlet::get_initial_cond()
    {
        return Payoff(BD_initial_cond);
    };

    double Boundaries_Dirichlet::get_final_cond()
    {
        return Payoff(BD_final_cond);
    }
    
    Boundaries_Dirichlet::~Boundaries_Dirichlet()
    {
        BD_initial_cond=0;
        BD_final_cond=0;
    }

    Boundaries_Neumann::Boundaries_Neumann(double initial_cond,double final_cond)
    :BN_initial_cond(initial_cond),
    BN_final_cond(final_cond)
    {
    };

    double Boundaries_Neumann::get_initial_cond()
    {
        return Payoff_df(BN_initial_cond);
    };

    double Boundaries_Neumann::get_final_cond()
    {
        return Payoff_df(BN_final_cond);
    }

    Boundaries_Neumann:: ~Boundaries_Neumann()
    {
        BN_initial_cond=0;
        BN_final_cond=0;
    }


}
