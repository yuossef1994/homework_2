//
//  nozzle_exact.cpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#include "nozzle_exact.hpp"
#include <iostream>
#include <math.h>
#include <fstream>

void nozzle_exact::initialization (double A_x, double P_0, double T_0, double Th_area)
{
    _A_x= A_x;
    _P_0= P_0;
    _T_0= T_0;
    _Th_area= Th_area;  
    A_par = _A_x/_Th_area;
    iter = 0;
    
    
}
void nozzle_exact::run_simulation(double M_initial)
{
    M = M_initial;
    //iteration loop
    do
    {
        iter = iter + 1;
        phi= (2/(gamma+1)) *(1+ M*M*((gamma-1)/2));
        F_M= pow( phi,(gamma+1)/(gamma-1)) - A_par*A_par*M*M;
        dF_M = 2*M*( pow(phi,2/(gamma-1)) - A_par*A_par);
        if (dF_M ==0)break;
        M = M - (F_M/dF_M);
     //displaying iteration on the output screen
        std::cout<< " iter number  "<< iter<< "  residual  "<<abs(F_M)<<"  mach number value  "<< M<<std::endl;
        
    }
    //defining a tolerence
    while(abs(F_M)> 1e-08);
}
void nozzle_exact::print_(double x)
{
    // calculating results
    epsi= 1 + ((gamma-1)*M*M)/2;
    T = _T_0/epsi;
    P= _P_0/pow(epsi,gamma/(gamma-1));
    rho = P/(R_air*T);
    u = M * sqrt(gamma*R_air*T);
  
   //write data into a file
  
   
   myfile<<std::setprecision(5)<<x<<","<<std::setprecision(5)<<u<<std::endl;
    
};
