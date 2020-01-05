//
//  Rates.hpp
//  Rates
//
//  Created by Florian on 04/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#ifndef Rates_hpp
#define Rates_hpp

#include <iostream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <string>
#include <random>
#include <vector>
#include <map>
#include <stdlib.h>
#include <math.h>

std::vector<double> Brownian_Motion(const int nb_sim);

std::vector<double> Vasicek_diffusion(const int maturity,double number_paths, double nb_sim, double khappa_r,double theta_r,
                                      double sigma_r,double r0);

std::vector<double> Rendleman_Batter_diffusion(const int maturity,double number_paths, double nb_sim, double theta_r,double sigma_r,double r0);

std::vector<double> Cox_Ingersoll_Ross_diffusion(const int maturity,double number_paths, double nb_sim, double khappa_r,double theta_r,
double sigma_r,double r0);


void print(std::vector<double> v); //Function to print values of a vector

#endif /* Rates_hpp */
