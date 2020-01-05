//
//  Rates.cpp
//  Rates
//
//  Created by Florian on 04/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#include "Rates.hpp"

std::default_random_engine generator;
std::normal_distribution<double> distribution(0.,1.);

std::vector<double> Brownian_Motion(const int nb_sim)
{
    std::vector<double> dBM(nb_sim);
    for(int i=0;i<nb_sim;i++){

        dBM[i]=distribution(generator);
        distribution.reset();
        
    }
    return dBM;
}

void print(std::vector<double> v)
{
    for(auto iter = v.begin(); iter != v.end(); ++iter)
    {
        std::cout<< *iter <<std::endl;
    }
    std::cout<< std::endl;
    
}

std::vector<double> Vasicek_diffusion(const int maturity,double number_paths,double nb_sim,double khappa_r,double theta_r,double sigma_r,double r0)
{
    std::vector<double> Vasicek(maturity*number_paths);
    double delta_t=1/number_paths;
    Vasicek[0]=r0;
    for(int i=1;i<maturity*number_paths;i++){
        
        std::vector<double> dBM=Brownian_Motion(nb_sim);
        Vasicek[i]=Vasicek[i-1]+khappa_r*(theta_r-Vasicek[i-1])*delta_t+sigma_r*
        std::accumulate(dBM.begin(), dBM.end(), 0.0)/dBM.size()*sqrt(delta_t);
    }
    return Vasicek;
}

std::vector<double> Rendleman_Batter_diffusion(const int maturity,double number_paths, double nb_sim, double theta_r,double sigma_r,double r0)
{
    std::vector<double> Rendleman_Batter(maturity*number_paths);
    double delta_t=1/number_paths;
    Rendleman_Batter[0]=r0;
    for(int i=1;i<maturity*number_paths;i++){
        
        std::vector<double> dBM=Brownian_Motion(nb_sim);
        Rendleman_Batter[i]=Rendleman_Batter[i-1]+theta_r*Rendleman_Batter[i-1]*delta_t+sigma_r*sqrt(delta_t)
        *std::accumulate(dBM.begin(), dBM.end(), 0.0)/dBM.size();
    }
    return Rendleman_Batter;
}

std::vector<double> Cox_Ingersoll_Ross_diffusion(const int maturity,double number_paths,double nb_sim,double khappa_r,double theta_r,double sigma_r,double r0)
{
    std::vector<double> Cox_Ingersoll_Ross(maturity*number_paths);
    double delta_t=1/number_paths;
    Cox_Ingersoll_Ross[0]=r0;
    for(int i=1;i<maturity*number_paths;i++){
        
        std::vector<double> dBM=Brownian_Motion(nb_sim);
        Cox_Ingersoll_Ross[i]=Cox_Ingersoll_Ross[i-1]+khappa_r*(theta_r-Cox_Ingersoll_Ross[i-1])*delta_t+
        sigma_r*sqrt(Cox_Ingersoll_Ross[i-1])*std::accumulate(dBM.begin(), dBM.end(), 0.0)/dBM.size()*sqrt(delta_t);
    }
    return Cox_Ingersoll_Ross;
}
