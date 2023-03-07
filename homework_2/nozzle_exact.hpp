//
//  nozzle_exact.hpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#ifndef nozzle_exact_hpp
#define nozzle_exact_hpp

#include <stdio.h>
#include <fstream>

class nozzle_exact {
    
    
private:
   //private variables that has been used in the main fn has _Var to distinquish it from the public variable in the main
    double _A_x , _P_0, _T_0, _Th_area, T, P, rho, u; //initialize area, reference (pressure, temp), throttle area, temperature(x), pressure(x), rho(x), u(x)
    double M; //mach number
    double epsi; //to make calculations easier
    double gamma = 1.4, R_u=  8314, M_air= 28.96; //gamma air, universal gas const, molecular weight of air
    double R_air = R_u/M_air;
    double phi; //to make calculations easier
    double A_par; //A/A^*
    double F_M; //F(M) M is the mach number
    double dF_M; //d(F(M)/dM
    int iter = 0; //variable to track num of iterations
    std::fstream myfile; //creation of file to write the results
   
public:
    //constructor to open the file for write
    nozzle_exact(){
        
        myfile.open("/Users/cringedaddy/CFD class/hw_2/homework_2/homework_2/uexact.txt", std::ios::trunc | std::ios::out);
        myfile<<"#########x########   "<<"#########u_exact########"<<std::endl;
    };
    // destructor to close the file
    ~nozzle_exact()
    {
        myfile.close();
        
    };
    //initialization functions that takes Area(x),pressure, temperature, throttle area
    void initialization(double, double, double, double);
    //to run iterations
    void run_simulation(double);
    //to print the results in file
    void print_(double);
   
    
    
};




#endif /* nozzle_exact_hpp */
