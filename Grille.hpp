//
//  Grille.hpp
//  Grille
//
//  Created by Florian on 04/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#ifndef Grille_hpp
#define Grille_hpp

#include <iostream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h>
#include <math.h>


    class Discretization
    {
        public:
        
         Discretization(const double &S_max, const double &number_path_mash,const double &S_min,const double &maturity, const double &time_path);
         std::vector<double> get_diffusion_underlying();
         std::vector<double> get_diffusion_time();
         ~Discretization();
        
        private:
        std::vector<double> log_ST;
        std::vector<double> vector_time;
        double number_path_underlying;
        double underlying_max;
        double underling_min;
        double path_time;
        double time_product;
    
        
    };
#endif /* Grille_hpp */
