#include <iostream>
#include "closed_form.hpp"
#include <string>
#include "Boundaries.hpp"
#include "Rates.hpp"


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



int main(int argc, char* argv[])
{
    // insert code here...
    dauphine::Dirichlet_boundaries_conditions();
    dauphine::Neumann_boundaries_conditions();
    Rates_diffusion("Vasicek model");
    return 0;
    
}
