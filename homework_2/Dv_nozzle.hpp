//
//  Dv_nozzle.hpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#ifndef Dv_nozzle_hpp
#define Dv_nozzle_hpp

#include <stdio.h>
#include <fstream>


class Dv_nozzle {
    
    
private:
    
    static int const _imax_f=51, _imax_c=50, _imin=0;
    int im_f, im_c;
    double A, B, C, D, E, F;  // just variables for calculations
    double k_2=0.3, k_4=0.1;
    double x_f[_imax_f];
    double avg_area;
    double avg_d_area;
    double x_c[_imax_c];
    double delta_t[_imax_c];
    double CFL=0.01;
    double area_f[_imax_f];
    double _area, _d_area;
    double mu[_imax_c];
    double epsilon_2;
    double _R[3]={0,0,0};
    double _R_iter[3]={0,0,0};
    double eigen_v_avg;
    double _delta_x;
    double _P_0,  _T_0, mach_in, mach_out, T_in, u_in, P_in, rho_in, et_in, T_out, u_out, P_out, rho_out, et_out;
    double _M_0;
    double _mach[_imax_c];
    double T[_imax_c];
    double et[_imax_c];
    double _V[3][_imax_c];
    double _U[3][_imax_c];
    double _F[3][_imax_f];
    double _D[3][_imax_f]; //JST dissipation scheme
    double gamma = 1.4, R_u=  8314, M_air= 28.96; //gamma air, universal gas const, molecular weight of air
    double R_air = R_u/M_air;
    double epsi;
    double rL2norm[3];
    double _rL2initial[3]={0,0,0};
    std::fstream L_vs_Iter; //creation of file to write the results
    
    std::fstream rho_vs_x; //creation of file to write the results
    std::fstream u_vs_x; //creation of file to write the results
    std::fstream p_vs_x; //creation of file to write the results
    
public:
    //constructor to open the file for write
    Dv_nozzle(){
        
     
        
        rho_vs_x.open("/Users/cringedaddy/CFD class/hw_2/homework_2/homework_2/rho.txt", std::ios::trunc | std::ios::out);
        
        u_vs_x.open("/Users/cringedaddy/CFD class/hw_2/homework_2/homework_2/u.txt", std::ios::trunc | std::ios::out);
        p_vs_x.open("/Users/cringedaddy/CFD class/hw_2/homework_2/homework_2/p.txt", std::ios::trunc | std::ios::out);
        L_vs_Iter.open("/Users/cringedaddy/CFD class/hw_2/homework_2/homework_2/L2.txt", std::ios::trunc | std::ios::out);
        L_vs_Iter<<"#iter#"<<"##L2_1#"<<"##L2_2#"<<"##L2_3#"<<std::endl;
        
    };
    // destructor to close the file
    ~Dv_nozzle()
    {
        
        rho_vs_x.close();
        u_vs_x.close();
        p_vs_x.close();
        L_vs_Iter.close();
        
    };

   
    void set_geometry(int,double,double);
    void initialization(double,double);
    void set_boundary_cond();
    void euler_explicit();
    double area(double);
    double d_area(double);
    void time_step();
    void print_res();
    void rL2initial();
    double L2norm(int);
};




#endif /* Dv_nozzle_hpp */
