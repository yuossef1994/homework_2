//
//  main.cpp
//  homework_2
//
//  Created by Youssef Z on 2/26/23.
//

#include <iostream>
#include "Dv_nozzle.hpp"
#include "nozzle_exact.hpp"
#include <math.h>

int main(int argc, const char * argv[]) {
    
    /*
    double A_x, T_0, P_0, Th_area, M_initial;
    T_0=600;
    P_0=300000;
    int max_faces=161;
    
    double x_min =-1;
    double x_max =1;
    int im_f_exact = max_faces-1;
   
    double *x= new double();
    double *x_f= new double();
   
   
   
   
   for( int i =0; i <= im_f_exact; i++)
   {
     x_f[i]=   x_min + (float(i)/float(im_f_exact))*(x_max- x_min) ;
       
   }
   
 
   
    
    nozzle_exact exact;
    
    for (int i=0;i<=im_f_exact;i++)
    {
        A_x = 0.2 + 0.4*(1 + sin(3.14*(x_f[i]-0.5)));
        Th_area = 0.2 + 0.4*(1 + sin(3.14*(x_f[(im_f_exact/2)]-0.5)));
        
        if ( i < 2)M_initial=0.5;
          
        else if (i == 2) M_initial=1;
        else M_initial = 5;
        
            
        exact.initialization(A_x, P_0, T_0, Th_area);
        exact.run_simulation(M_initial);
        exact.print_(x[i]);
        
    }
    
    
    */
    
    
    
    
    
    int max_iter=10000000;
    double L2norm;
    Dv_nozzle nozzle;
    nozzle.set_geometry(200, 1, -1);
    nozzle.initialization(300000, 600);
    nozzle.set_boundary_cond();
    nozzle.time_step();
    nozzle.euler_explicit();
    nozzle.rL2initial();
    
    for(int i=1;i<=max_iter;i++)
    {
        
        nozzle.set_boundary_cond();
        nozzle.time_step();
        nozzle.euler_explicit();
       
        L2norm= nozzle.L2norm(i);
        
        
        std::cout<< " iter number  "<<i <<"  L2 norm is  " <<L2norm<<std::endl;
        
        
        if(L2norm < 1e-8){
            
            std::cout<< "solution converged"<<std::endl;
            break;
        }
       
    }
    
    
    nozzle.print_res();
    
    
    return 0;
}
