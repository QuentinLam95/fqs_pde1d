//
//  pde_solver.cpp
//  PDE_solver
//
//  Created by Florian on 08/01/2020.
//  Copyright © 2020 Florian. All rights reserved.
//

#include "pde_solver.hpp"
#include "Grille.hpp"
#include "closed_form.hpp"
#include "matrix.hpp"

namespace Solve
    {
    matrix_pde_case1::matrix_pde_case1(std::size_t maturity, double time_path, double sigma, double rate, double theta,  std::size_t number_path_mash, double S_max, double S_min)
    :m_underlying_nb_dimension(number_path_mash),m_maturity_nb_dimension(maturity)
    {
    }
    
    std::size_t  matrix_pde_case1::number_path_underlying()
    
    {
        return m_underlying_nb_dimension;
    }
    
    std::size_t  matrix_pde_case1::m_maturity()
       
       {
           return m_maturity_nb_dimension;
       }
    
    
    
    //if (side==true) A mettre dans un autre code
    //{
       // double temp=theta*maturity/time_path;
    //}
    //else
    //{
      //  double temp=-(1-theta)*maturity/time_path;
    //}
    
    std::vector<double> matrix_pde_case1::forward_coefficient(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min, double temp)
    
    {
        
        double dx=(S_max-S_min)/number_path_mash;
        double dx_2=pow(dx,2.0);
        double sigma_2=pow(sigma,2.0);
        std::vector<double>a_coefficient(number_path_mash);
        
        
         for(auto path_mash =0; path_mash!= number_path_mash; ++path_mash)
         {
             if (path_mash==0){
                 a_coefficient[path_mash]=-sigma_2/(2.0 * dx_2);
                 a_coefficient[path_mash]=temp*a_coefficient[path_mash];
             }
             else if(path_mash==number_path_mash-1)
             {
                 a_coefficient[path_mash]=rate-sigma_2/(2*dx_2)+
                 (0.5*sigma_2-rate)/dx;
                 a_coefficient[path_mash]=temp*a_coefficient[path_mash];
             }
             else
            {
                a_coefficient[path_mash]=(0.5*sigma_2-rate)/(2*dx)-0.5*sigma_2/dx_2;
                a_coefficient[path_mash]=temp*a_coefficient[path_mash];
            }
             
         }
        return a_coefficient;
    }
    
    std::vector<double> matrix_pde_case1::present_coefficient(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,double temp)
       {
           double dx=(S_max-S_min)/number_path_mash;
           double dx_2=pow(dx,2.0);
           double sigma_2=pow(sigma,2.0);
           std::vector<double>b_coefficient(number_path_mash);
            for(auto path_mash =0; path_mash!= number_path_mash; ++path_mash)
            {
                if (path_mash==0){
                    b_coefficient[path_mash]=-sigma_2/dx_2+(0.5*sigma_2-rate)/dx;
                    b_coefficient[path_mash]=1+temp*b_coefficient[path_mash];
                }
                else if(path_mash==number_path_mash-1)
                {
                    b_coefficient[path_mash]=sigma_2/dx_2-
                    (0.5*sigma_2-rate)/dx;
                    b_coefficient[path_mash]=1+temp*b_coefficient[path_mash];
                }
                else
               {
                   b_coefficient[path_mash]=rate+sigma_2/dx_2;
                   b_coefficient[path_mash]=1+temp*b_coefficient[path_mash];
               }

                
            }
           return b_coefficient;
       }
    
    std::vector<double> matrix_pde_case1::backward_coefficient(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,double temp)
       {
           double dx=(S_max-S_min)/number_path_mash;
           double dx_2=pow(dx,2.0);
           double sigma_2=pow(sigma,2.0);
           std::vector<double>c_coefficient(number_path_mash);
            for(auto path_mash =0; path_mash!= number_path_mash; ++path_mash)
            {
                if (path_mash==0){
                    c_coefficient[path_mash]=rate-(0.5*sigma_2-rate)/dx-sigma_2/(2*dx_2);
                    c_coefficient[path_mash]=temp*c_coefficient[path_mash];
                }
                else if(path_mash==number_path_mash-1)
                {
                    c_coefficient[path_mash]=-sigma_2/(2*dx_2);
                    c_coefficient[path_mash]=temp*c_coefficient[path_mash];
                }
                else
               {
                   c_coefficient[path_mash]=-0.5*sigma_2/dx_2-(0.5*sigma_2-rate)/(2*dx);
                   c_coefficient[path_mash]=temp*c_coefficient[path_mash];
               }
                
            }
           return c_coefficient;
       }
    
    
    
    
    std::vector<double> matrix_pde_case1::matrix_rhs(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,  double theta)
        {
            // Matrix of the right side, so with fn
                       double temp_rhs=-(1-theta)*maturity/time_path;
                       std::vector<double> a_rhs=matrix_pde_case1::forward_coefficient(maturity, time_path,sigma, rate, number_path_mash,S_max, S_min, temp_rhs);
                       std::vector<double> b_rhs=matrix_pde_case1::present_coefficient(maturity, time_path,sigma, rate, number_path_mash,S_max, S_min, temp_rhs);
                       std::vector<double> c_rhs=matrix_pde_case1::backward_coefficient(maturity, time_path,sigma, rate, number_path_mash,S_max, S_min, temp_rhs);
                       
                       dauphine::matrix m_rhs(number_path_mash,number_path_mash);
            
            
            
            Discretization vector_time_underlying(S_max, number_path_mash,S_min,maturity, time_path);
                       
                       
                   std::vector<double> St=vector_time_underlying.get_diffusion_underlying();
                   std::vector<double> Result(number_path_mash);
               
                       std::vector<double> Payoff_N(number_path_mash);
                   for(auto path_mash =0; path_mash!= number_path_mash; ++path_mash)
                       {
                           Payoff_N[path_mash]=dauphine::vanilla_payoff(exp(St[path_mash]),100., true);
                       }
                       
                   for(auto i =0; i!= number_path_mash; ++i)
                       {
                           if (i==0){
                               m_rhs(i,i)=b_rhs[i];
                               m_rhs(i,i+1)=a_rhs[i];
                           }
                           else
                           {
                               m_rhs(i,i-1)=c_rhs[i];
                               m_rhs(i,i)=b_rhs[i];
                               m_rhs(i,i+1)=a_rhs[i];
                           }
                           
                       }
                       Result=m_rhs.produit_mat_vect(m_rhs,Payoff_N);
                       
                       
                       std::vector<double> Payoff_Boundaries(number_path_mash);
                       for(auto i =0; i!= number_path_mash; ++i)
                       {
                           if (i==0 or i== number_path_mash){
                               Payoff_Boundaries[i]=c_rhs[i]*1+Result[i];
                           }
                           else{
                               Payoff_Boundaries[i]=a_rhs[i]*1+Result[i];
                           }
                       }
                       //std::cout<<m_rhs<<std::endl;
                       //std::cout<<U<<std::endl;
                       return Payoff_Boundaries;
                    }
    
        std::vector<double> matrix_pde_case1::LU_decomposition(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,  double theta)
        {
        
        
             // Matrix of the left side, so with fn
             double temp_lhs=theta*maturity/time_path;
             std::vector<double> a_lhs=matrix_pde_case1::forward_coefficient(maturity, time_path,sigma, rate, number_path_mash,S_max, S_min, temp_lhs);
             std::vector<double> b_lhs=matrix_pde_case1::present_coefficient(maturity, time_path,sigma, rate, number_path_mash,S_max, S_min, temp_lhs);
             std::vector<double> c_lhs=matrix_pde_case1::backward_coefficient(maturity, time_path,sigma, rate, number_path_mash,S_max, S_min, temp_lhs);
        
             
             dauphine::matrix m_lhs(number_path_mash,number_path_mash);
         for(auto i =0; i!= number_path_mash; ++i)
             {
                 if (i==0){
                     m_lhs(i,i)=b_lhs[i];
                     m_lhs(i,i+1)=a_lhs[i];
                 }
                 else
                 {
                     m_lhs(i,i-1)=c_lhs[i];
                     m_lhs(i,i)=b_lhs[i];
                     m_lhs(i,i+1)=a_lhs[i];
                 }
             }
                 
                 
        dauphine::matrix L(number_path_mash,number_path_mash);
        dauphine::matrix U(number_path_mash,number_path_mash);
        std::vector<double> z(number_path_mash);
        L(0,0)=m_lhs(0,0);
        U(0,1)=m_lhs(0,1)/L(0,0);
                 
        for(auto i =1; i!= number_path_mash-1; ++i)
        {
            L(i,i-1)=m_lhs(i,i-1);
            L(i,i)=m_lhs(i,i)-L(i,i-1)*U(i-1,i);
            U(i,i+1)=m_lhs(i,i+1)/L(i,i);
        }
                 
    L(number_path_mash-1,number_path_mash-2)=m_lhs(number_path_mash-1,number_path_mash-2);
    L(number_path_mash-1,number_path_mash-1)=m_lhs(number_path_mash-1,number_path_mash-1)-L(number_path_mash-1,number_path_mash-2)*U(number_path_mash-1,number_path_mash-1);
        
         std::vector<double> b=matrix_pde_case1::matrix_rhs(maturity, time_path, sigma, rate, number_path_mash, S_max, S_min,  theta);
        z[0]=b[0]/L(0,0);

        for(auto i =1; i!= number_path_mash; ++i)
        {
            z[i]=(b[i]-L(i,i-1)*z[i-1])/L(i,i);
        }
        
            //std::cout<<L<<std::endl;
            //std::cout<<U<<std::endl;
            std::vector<double> x(number_path_mash);
        
            x[number_path_mash]=z[number_path_mash];
        for(auto i =number_path_mash; i!= 0; --i)
        {
            x[i]=z[i]-U(i,i+1)*x[i+1];
        }
        return x;
    }
    
    std::vector<double> matrix_pde_case1::Final_solve(std::size_t maturity, double time_path,double sigma, double rate, double number_path_mash,double S_max, double S_min,  double theta)
    {
        dauphine::matrix Price(number_path_mash,maturity);
        std::vector<double> temp(number_path_mash);
          for(auto time =0; time!= maturity; ++time)
              {
                  temp=matrix_pde_case1::LU_decomposition(time, time_path, sigma, rate, number_path_mash,S_max, S_min, theta);
            for(auto mash =0; mash!= number_path_mash; ++mash)
                        {
                            Price(mash,time)=temp[mash];
              }
              }
            std::cout<<Price;
            return temp;
    }
}



