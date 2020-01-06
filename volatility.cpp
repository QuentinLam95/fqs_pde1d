//
//  volatility.cpp
//  Volatility
//
//  Created by Quentin Lam on 02/01/2020.
//  Copyright © 2020 Quentin Lam. All rights reserved.
//

#include "volatility.hpp"

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
    
std::vector<double> Heston_diffusion(const int maturity,double number_paths,double LT_vol,double reversion_rate,double vol_vol,double v0)

{
    std::vector<double> volatility_heston(maturity*number_paths);
    double delta_t=1/number_paths;
    volatility_heston[0]=v0;
    
    std::vector<double> dBM=Brownian_Motion(maturity*number_paths);
    
    for(int i=1;i<maturity*number_paths;i++)
    {
        double vol_max = std::max(volatility_heston[i-1], 0.0); //to avoid negative volatility
        volatility_heston[i]=volatility_heston[i-1]+reversion_rate*(LT_vol-vol_max)*delta_t+vol_vol*sqrt(vol_max* delta_t)*dBM[i];
        
    }
    return volatility_heston;
    
}

