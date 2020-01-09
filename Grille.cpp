//
//  Rates.cpp
//  Rates
//
//  Created by Florian on 04/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#include "Grille.hpp"



    Discretization::Discretization(const double &S_max, const double &number_path_mash,const double &S_min,const double &maturity, const double &time_path)
      {
        std::vector<double> log_St(number_path_mash+1);
        double d_mash=(S_max-S_min)/number_path_mash;
        log_St[0]=log(S_min);
        log_St[number_path_mash]=log(S_max);
          
          
          for(int i=1;i<number_path_mash;i++)
          {
              log_St[i]=exp(log_St[i-1])+d_mash;
              log_St[i]=log(log_St[i]);
          }
            log_ST=log_St;
          
          std::vector<double> time(time_path+1);
          double d_time=maturity/time_path;
          time[0]=0;
          time[time_path]=maturity;
          
          
          for(int t=1;t<time_path;t++)
          {
              time[t]=t*d_time;
          }
          vector_time=time;
          
    }

    std::vector<double> Discretization::get_diffusion_underlying()
    {
        return log_ST;
    }

    std::vector<double> Discretization::get_diffusion_time()
    {
        return vector_time;
    }
    Discretization::~Discretization()
    {
    }




