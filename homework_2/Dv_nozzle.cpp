//
//  Dv_nozzle.cpp
//  hw_1
//
//  Created by Youssef Z on 2/11/23.
//

#include "Dv_nozzle.hpp"
#include <iostream>
#include <math.h>
#include <fstream>
using namespace std;



void Dv_nozzle::set_geometry (int imax, double x_max, double x_min)
{
    
    //imax is number of cells
    //_imax is number of faces
    
     im_f = _imax_f-1;
     im_c = _imax_c-1;
    
    
    
    
    for( int i =_imin; i <= im_f; i++)
    {
      x_f[i]=   x_min + (float(i)/float(im_f))*(x_max- x_min) ;
        
    }
    
  
    for( int j =_imin; j <= im_c; j++)
    {
        x_c[j] =   (x_f[j] + x_f[j+1])/2;
        
    }
   
  
}
    void Dv_nozzle::initialization (double P_0, double T_0)
{
        
        
        
        _P_0= P_0;
        _T_0= T_0;
        
        
        
        
        
        
        
        
        
        // 0 is rho, 1 is velocity,  2 is pressure
        //values calculated at cell as cell average
        for (int i= _imin; i<= im_c; i++)
        {
            
            _M_0 = 0.85*x_c[i] + 1;
            epsi= 1 + ((gamma-1)*_M_0*_M_0)/2;
            T[i] = _T_0/epsi;
            
            _V[1][i] = _M_0 * sqrt(gamma*R_air*abs(T[i]));
            _V[2][i]= _P_0/pow(epsi,gamma/(gamma-1));
            _V[0][i]= _V[2][i]/(R_air*T[i]);
            et[i] = (R_air/(gamma-1))*T[i] + 0.5* _V[1][i]*_V[1][i];
            _U[0][i]= _V[2][i]/(R_air*T[i]);
            _U[1][i]= _V[0][i]*_V[1][i];
            _U[2][i]=_V[0][i]*et[i];
            
            
            
            _mach[i]=_M_0;
            
        }
        for (int i= _imin+1; i<= im_c-1; i++)
        {
            
                mu[i]= abs(    ( _V[2][i+1] -2*_V[2][i]  + _V[2][i-1]  )  / ( _V[2][i+1] +2*_V[2][i]  + _V[2][i-1]  )   );
                
           
            
        }
        
       
        
        
        // values calculated at the faces
        for (int i= 1; i<= im_f-1; i++)
        {
            
            
            A= ( _V[0][i]*_V[1][i] + _V[0][i-1]*_V[1][i-1] )/2;
            B= ( (_V[0][i]*_V[1][i]*_V[1][i] + _V[2][i]) + (_V[0][i-1]*_V[1][i-1]*_V[1][i-1] + _V[2][i-1]) )/2;
            
            C= (gamma/(gamma-1))*_V[1][i]*_V[2][i] + (_V[0][i]*_V[1][i]*_V[1][i]*_V[1][i])/2;
            D=(gamma/(gamma-1))*_V[1][i-1]*_V[2][i-1] + (_V[0][i-1]*_V[1][i-1]*_V[1][i-1]*_V[1][i-1])/2;
            
            _F[0][i] = A;
            _F[1][i]= B;
            _F[2][i]= (C + D)/2;
           
        }
        
        
        
        
        
        
        
        
        
        //damping
        for (int i= _imin+2; i<= im_f-2; i++)
        {
            // damping make a different loop here cuz two values are extrapolated
            eigen_v_avg = ( abs(_V[1][i+1]) + abs(_V[1][i]) + sqrt(gamma*R_air*abs(T[i+1])) + sqrt(gamma*R_air*abs(T[i])) )/(2);
            
            if(i== _imin+2)
            {
                epsilon_2 = k_2*std::max(std::max(mu[i-1],mu[i]) ,mu[i+1]);
            }
            else if (i==im_f-2)
            {
                //epsilon_2 = k_2*(std::max(std::max(mu[i-2],mu[i-1]),mu[i]));
                epsilon_2 = k_2*(std::max(mu[i-1]),mu[i]);
            }
             else
             {
                 //epsilon_2 = k_2*std::max(std::max(std::max(mu[i-2],mu[i-1]),mu[i]) ,mu[i+1]);
                 epsilon_2 = k_2*std::max(mu[i],mu[i+1]);
              };
            
            
            
            
            
            
            _D[0][i]= -eigen_v_avg*epsilon_2*( _U[0][i] - _U[0][i-1]) +  eigen_v_avg*std::max(0.0,(k_4-epsilon_2))*(_U[0][i+1] -3*_U[0][i]+3*_U[0][i-1]-_U[0][i-2]);
            _D[1][i]= -eigen_v_avg*epsilon_2*( _U[1][i] - _U[1][i-1]) +  eigen_v_avg*std::max(0.0,(k_4-epsilon_2))*(_U[1][i+1] -3*_U[1][i]+3*_U[1][i-1]-_U[1][i-2]);
            _D[2][i]= -eigen_v_avg*epsilon_2*( _U[2][i] - _U[2][i-1]) +  eigen_v_avg*std::max(0.0,(k_4-epsilon_2))*(_U[2][i+1] -3*_U[2][i]+3*_U[2][i-1]-_U[2][i-2]);
            
        }
            
       
        
        
        
        
        
    }

void Dv_nozzle::set_boundary_cond()
{
    //inflow boundary conditions
    mach_in = 0.5*(3.0*_mach[0] - _mach[1] );
    epsi= 1 + ((gamma-1)*mach_in*mach_in)/2;
    T_in = _T_0/epsi;
    
    
    u_in =  mach_in * sqrt(gamma*R_air*abs(T_in));
    
    P_in= _P_0/pow(epsi,gamma/(gamma-1));
    rho_in= P_in/(R_air*T_in);
    
    et_in = (R_air/(gamma-1))*T_in + 0.5*u_in*u_in;
    
    _F[0][0] = rho_in*u_in ;
    _F[1][0]= rho_in*u_in*u_in + P_in;
    _F[2][0]= (gamma/(gamma-1))*u_in*P_in + (rho_in*u_in*u_in*u_in)/2;
  
   
    
    
    
    
    // damping
    //extrapolate face at 1
    _D[0][1]=2*_D[0][2] - _D[0][3];
    _D[1][1]=2*_D[1][2] - _D[1][3];
    _D[2][1]=2*_D[2][2] - _D[2][3];
   //extrapolate face at 0
    _D[0][0]=2*_D[0][1] - _D[0][2];
    _D[1][0]=2*_D[1][1] - _D[1][2];
    _D[2][0]=2*_D[2][1] - _D[2][2];
   
//******************************************
  // outflow boundary conditions
    P_out = 120000;
    /*
    rho_out = 0.5*(3*_V[0][im_c]-_V[0][im_c-1]);
    u_out = 0.5*(3*_V[1][im_c]-_V[1][im_c-1]);
    T_out = 0.5*(3*T[im_c]-T[im_c-1]);
                               
    et_out = 0.5*(3*et[im_c]-et[im_c-1]);
   */
   // P_out = 0.5*(3*_V[2][im_c]-_V[2][im_c-1]);
   // rho_out = 0.5*(3*_U[0][im_c]-_U[0][im_c-1]);
    
    
   // u_out = (0.5*(3*_U[1][im_c]-_U[1][im_c-1]))/rho_out;
    
    /*
    T_out = 0.5*(3*T[im_c]-T[im_c-1]);
                               
    et_out = 0.5*(3*et[im_c]-et[im_c-1]);
    
    rho_out = _V[0][im_c];
    u_out = _V[1][im_c];
    T_out = T[im_c];
                               
    et_out = et[im_c];
    */
    //rho_out = _V[0][im_c];
   // u_out = _V[1][im_c];
    rho_out = 0.5*(3*_V[0][im_c]-_V[0][im_c-1]);
    u_out = 0.5*(3*_V[1][im_c]-_V[1][im_c-1]);
    _F[0][im_f] = rho_out*u_out ;
    _F[1][im_f]= rho_out*u_out*u_out + P_out;
    _F[2][im_f]= (gamma/(gamma-1))*u_out*P_out + (rho_out*u_out*u_out*u_out)/2;
    
    // damping
    //extrapolate face at im_f-1
    _D[0][im_f-1]=2*_D[0][im_f-2] - _D[0][im_f-3];
    _D[1][im_f-1]=2*_D[1][im_f-2] - _D[1][im_f-3];
    _D[2][im_f-1]=2*_D[2][im_f-2] - _D[2][im_f-3];
   //extrapolate face at imf
    _D[0][im_f]=2*_D[0][im_f-1] - _D[0][im_f-2];
    _D[1][im_f]=2*_D[1][im_f-1] - _D[1][im_f-2];
    _D[2][im_f]=2*_D[2][im_f-1] - _D[2][im_f-2];
    
    
}

void Dv_nozzle::euler_explicit()
{
    
    for (int i=0;i<=im_c;i++)
    {
        _delta_x=abs(x_f[i+1]-x_f[i]);
        avg_area = ((area(x_f[i+1])+area(x_f[i]))/2);
        avg_d_area = ((d_area(x_f[i+1])+d_area(x_f[i]))/2);
        
    _U[0][i]= ( (-_F[0][i+1]- _D[0][i+1]   ) *area(x_f[i+1]) + (   _F[0][i] + _D[0][i]   )*area(x_f[i])    )*
               
               
              ( delta_t[i]/(avg_area *_delta_x)) + _U[0][i];
           
        
    _U[1][i]= (  (-_F[1][i+1]-_D[1][i+1]    )*area(x_f[i+1]) + (   _F[1][i] + _D[1][i]   )*area(x_f[i]) +_V[2][i]*avg_d_area*_delta_x    )*(delta_t[i]/(avg_area*_delta_x)) + _U[1][i];
        
        
    _U[2][i]= (( -_F[2][i+1] -_D[2][i+1]) *area(x_f[i+1]) + (_F[2][i] + _D[2][i]  )*area(x_f[i]))*(delta_t[i]/(avg_area *_delta_x)) + _U[2][i];
           
           
    }
    
    for (int i=0;i<=im_c;i++)
    {
       
        _V[0][i]=_U[0][i];
        _V[1][i]= _U[1][i]/_V[0][i];
        
        et[i] = _U[2][i]/ _V[0][i];
        T[i] = (et[i] - 0.5* _V[1][i]*_V[1][i])*((gamma-1)/R_air);
        _V[2][i]=_V[0][i]*R_air*T[i];
        _mach[i]=_V[1][i]/sqrt(gamma*R_air*abs(T[i]));
       
    }
    
    for (int i= _imin+1; i<= im_c-1; i++)
    {
        
            mu[i]= abs(    ( _V[2][i+1] -2*_V[2][i]  + _V[2][i-1]  )  / ( _V[2][i+1] +2*_V[2][i]  + _V[2][i-1]  )   );
            
       
        
    }
    
    
    
    
    
    
    
    
    //calculate the flux
    for (int i=1;i<=im_c;i++)
    {
       
        A= ( _V[0][i]*_V[1][i] + _V[0][i-1]*_V[1][i-1] )/2;
        B= ( (_V[0][i]*_V[1][i]*_V[1][i] + _V[2][i]) + (_V[0][i-1]*_V[1][i-1]*_V[1][i-1] + _V[2][i-1]) )/2;
        
        C= (gamma/(gamma-1))*_V[1][i]*_V[2][i] + (_V[0][i]*_V[1][i]*_V[1][i]*_V[1][i])/2;
        D=(gamma/(gamma-1))*_V[1][i-1]*_V[2][i-1] + (_V[0][i-1]*_V[1][i-1]*_V[1][i-1]*_V[1][i-1])/2;
        
        _F[0][i] = A;
        _F[1][i]= B;
        _F[2][i]= (C + D)/2;
    }
    //calculate the damping term
    for (int i= _imin+2; i<= im_f-2; i++)
    {
        // damping make a different loop here cuz two values are extrapolated

        eigen_v_avg = ( abs(_V[1][i+1]) + abs(_V[1][i]) + sqrt(gamma*R_air*abs(T[i+1])) + sqrt(gamma*R_air*abs(T[i])) )/(2*CFL);
        if(i== _imin+2)
        {
            epsilon_2 = k_2*std::max(std::max(mu[i-1],mu[i]) ,mu[i+1]);
        }
        else if (i==im_f-2)
        {
            epsilon_2 = k_2*(std::max(std::max(mu[i-2],mu[i-1]),mu[i]));
        }
        
         else
         {
             epsilon_2 = k_2*std::max(std::max(std::max(mu[i-2],mu[i-1]),mu[i]) ,mu[i+1]);
         };
        
        
        _D[0][i]= -eigen_v_avg*epsilon_2*( _U[0][i] - _U[0][i-1]) +  eigen_v_avg*std::max(0.0,(k_4-epsilon_2))*(_U[0][i+1] -3*_U[0][i]+3*_U[0][i-1]-_U[0][i-2]);
        _D[1][i]= -eigen_v_avg*epsilon_2*( _U[1][i] - _U[1][i-1]) +  eigen_v_avg*std::max(0.0,(k_4-epsilon_2))*(_U[1][i+1] -3*_U[1][i]+3*_U[1][i-1]-_U[1][i-2]);
        _D[2][i]= -eigen_v_avg*epsilon_2*( _U[2][i] - _U[2][i-1]) +  eigen_v_avg*std::max(0.0,(k_4-epsilon_2))*(_U[2][i+1] -3*_U[2][i]+3*_U[2][i-1]-_U[2][i-2]);
    
        
    }
    
}
double Dv_nozzle::area(double x)
{
    _area = 0.2 + 0.4*(1 + sin(3.14159265359*(x-0.5)));
    
    
    return _area;
    
}
double Dv_nozzle::d_area(double x)
{
    _d_area = 0.4*3.14159265359*cos(3.14159265359*(x-0.5));
    
    
    return _d_area;
    
}
void Dv_nozzle::time_step()
{
   
    for (int i =0;i<=im_c;i++)
    {
        _delta_x=abs(x_f[i+1]-x_f[i]);
        delta_t[i] = (CFL *_delta_x)/(abs(_V[1][i])  + sqrt(gamma*R_air*abs(T[i]))) ;
    }
    
  
    
}

void Dv_nozzle::print_res(){
    rho_vs_x<<"######### x ######   "<<"#########rho######   "<<std::endl;
    u_vs_x<<"######### x ######   "<<"#########ux######   "<<std::endl;
    p_vs_x<<"######### x ######   "<<"#########P######   "<<std::endl;
    for (int i= 0;i<=im_c;i++)
    {
        
        rho_vs_x<<std::setprecision(5)<<x_c[i]<<","<<std::setprecision(5)<< _V[0][i]<<std::endl;
        
        
        
        u_vs_x<<std::setprecision(5)<<x_c[i]<<","<<std::setprecision(5)<< _V[1][i]<<std::endl;
        
        
        
        p_vs_x<<std::setprecision(5)<<x_c[i]<<","<<std::setprecision(5)<< _V[2][i]<<std::endl;
        
        
        
    }
    
    
    
}
double Dv_nozzle::L2norm(int iter){
    
    
    rL2norm[0]=0;
    rL2norm[1]=0;
    rL2norm[2]=0;
    
    for (int i=0;i<=im_c;i++)
    {
        _delta_x=abs(x_f[i+1]-x_f[i]);
        _R_iter[0]= (_F[0][i+1]+_D[0][i+1]) *area(x_f[i+1]) -(_F[0][i] + _D[0][i]   )*area(x_f[i])  ;
        
        rL2norm[0] = rL2norm[0] +_R_iter[0]*_R_iter[0];
        
        
        _R_iter[1]= ( _F[1][i+1]+_D[1][i+1]    )*area(x_f[i+1]) -(_F[1][i] + _D[1][i]   )*area(x_f[i]) -_V[2][i]*( ((d_area(x_f[i+1])+d_area(x_f[i]))/2)  )*_delta_x  ;
        
        rL2norm[1] = rL2norm[1] +  _R_iter[1]* _R_iter[1];
        
        _R_iter[2]= ( _F[2][i+1] +_D[2][i+1]    )*area(x_f[i+1])- (_F[2][i] + _D[2][i]  )*area(x_f[i]);
        
        rL2norm[2] = rL2norm[2] +  _R_iter[2]* _R_iter[2];
    }
    double imax_c = _imax_c*1.0;
    
    rL2norm[0]=sqrt(rL2norm[0]/imax_c)/_rL2initial[0];
    rL2norm[1]=sqrt(rL2norm[1]/imax_c)/_rL2initial[1];
    rL2norm[2]=sqrt(rL2norm[2]/imax_c)/_rL2initial[2];
    
    
    std::cout<< " R eqn1 "<< _R_iter[0]<<std::endl;
    std::cout<< " R eqn2 "<< _R_iter[1]<<std::endl;
    std::cout<< " R eqn3 "<< _R_iter[2]<<std::endl;
    std::cout<< " L eqn1 "<< rL2norm[0]<<std::endl;
    std::cout<< " L eqn2 "<< rL2norm[1]<<std::endl;
    std::cout<< " L eqn3 "<< rL2norm[2]<<std::endl;
    std::cout<< " L inti eqn1 "<< _rL2initial[0]<<std::endl;
    std::cout<< " L inti eqn2 "<< _rL2initial[1]<<std::endl;
    std::cout<< " L inti eqn3 "<< _rL2initial[2]<<std::endl;
   // rL2norm[0] = sqrt(rL2norm[0]/imax_c)/_rL2initial[0];
  //  rL2norm[1] = sqrt(rL2norm[1]/imax_c)/_rL2initial[1];
  //  rL2norm[2] = sqrt(rL2norm[2]/imax_c)/_rL2initial[2];
    L_vs_Iter<<std::setprecision(5)<<iter<<","<<std::setprecision(5)<<rL2norm[0]<<","<<std::setprecision(5)<<rL2norm[1]<<","<<std::setprecision(5)<<rL2norm[2]<<std::endl;
    
    return std::max(rL2norm[0],std::max(rL2norm[1],rL2norm[2]));
    
    
}
void Dv_nozzle::rL2initial()


{
    _rL2initial[0]=0;
    _rL2initial[1]=0;
    _rL2initial[2]=0;
    for (int i=0;i<=im_c;i++)
    {
        _R[0]= (_F[0][i+1]+ _D[0][i+1]) *area(x_f[i+1]) - (_F[0][i] + _D[0][i])*area(x_f[i])  ;
        
        
        _rL2initial[0] = _rL2initial[0] + _R[0]*_R[0];
        
        
        _R[1]= ( _F[1][i+1]+_D[1][i+1]    )*area(x_f[i+1]) - (_F[1][i] + _D[1][i] )*area(x_f[i]) -_V[2][i]*(((d_area(x_f[i+1])+d_area(x_f[i]))/2)  )*_delta_x  ;
        
        _rL2initial[1] = _rL2initial[1] + _R[1]*_R[1];
        
        _R[2]= ( _F[2][i+1] +_D[2][i+1]    )*area(x_f[i+1]) -(_F[2][i] + _D[2][i] )*area(x_f[i]);
        
        _rL2initial[2] = _rL2initial[2] + _R[2]*_R[2];
    }
    double imax_c = _imax_c*1.0;
    _rL2initial[0] = sqrt(_rL2initial[0]/imax_c);
    _rL2initial[1] = sqrt(_rL2initial[1]/imax_c);
    _rL2initial[2] = sqrt(_rL2initial[2]/imax_c);
    
  
    
    
    
}
