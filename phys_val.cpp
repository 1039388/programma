#include <iostream>
#include <cmath>
#include <math.h>
#include <chrono>
#include <fstream>
#include <filesystem> 
#include <string>
#include "vcpkg\installed\x64-windows\include\boost\math\quadrature\trapezoidal.hpp"
#include "vcpkg\installed\x64-windows\include\boost\math\quadrature\tanh_sinh.hpp"
#include "vcpkg\installed\x64-windows\include\boost\math\quadrature\exp_sinh.hpp"
#include "vcpkg\installed\x64-windows\include\boost\math\quadrature\gauss.hpp"
#include "vcpkg\installed\x64-windows\include\boost\math\differentiation\finite_difference.hpp"
#include "vcpkg\installed\x64-windows\include\boost\math\quadrature\gauss_kronrod.hpp"
#include "gnuplot-iostream.h"
#include "vcpkg\installed\x64-windows\include\boost\numeric\odeint.hpp"
#include "vcpkg\installed\x64-windows\include\boost\math\complex.hpp"
#include <complex>
#define M_PI 3.14159265358979323846

using namespace std::complex_literals;



double pl=1;   //плотность
double mu=1;


std::complex<double> sigma=1.38*10000*pow(10,12)/1.113; // 3-4 множители перевод в СГС 1.38*10000(Ом*см)^-1
double tau=2.5*pow(10,-14);
double c_0=3*pow(10,10);



double R=1.5;
double B_s=1;
double M=4;
double q=4;
double k=4;
double bbbeta=0.96;
double p_0=bbbeta*2*R*R/(2*(bbbeta+R*R-1));
// double p_0=0.6501;   


double L=350;
double N_c=pow(10,13);
double B_v_=3000;
double B_s_=B_v_*R;
double m_i=940*pow(10,6)*1.6*pow(0.1,12); // масса покоя в эрг

double w_0=B_s_/(L*sqrt(N_c*m_i/c_0/c_0));
double c=c_0/w_0;
//608053020
double tol = 1e-8;
int max_refine = 15;
std::complex<double> dzeta(std::complex<double> omega) {
    using namespace std::complex_literals;
std::complex<double> epsilon=1.0+1.0i*6.55*pow(10,18)/w_0/omega;
//  return ((1.0 - 1.0i)* std::sqrt(omega * mu /(8 * M_PI * sigma)));
//  return std::sqrt(mu /(1.0+1.0i*4.0 * M_PI * sigma/(1.0-1.0i*omega*tau)));
return sqrt(mu/epsilon);
//    return 0.0;
}

double f(double psi){
    if(psi >= 0 && psi <= 1){
        return 1-pow(psi,k);
    }
    else{return 0;
    }
}
double p(double psi){
    return p_0*f(psi)/(2*R*R);
}
double B_v(double z){
    return (1+(M-1)*pow(sin(M_PI*z/2),q))/R;
}
double B_v2(double z){
    return 1.0/(B_v(z)*B_v(z));
}
double d_B_v(double z){
    return (M-1)*q*M_PI/2*cos(M_PI*z/2)*pow(sin(M_PI*z/2),q-1)/R;
}
double b_v(double z){
    return (1+(M-1)*pow(sin(M_PI*z/2),q))/R;
}
double B(double z,double psi ){
    double bb_v=b_v(z);
    if(bb_v>1){return bb_v;}
    else{return sqrt((bb_v*bb_v*R*R-p(psi)*2*(R*R))/(R*R-p(psi)*2*(R*R)));}
}
double a_v(double z){
    return sqrt(2/B_v(z));
} 
double d_a_v(double z){
    return -sqrt(2/B_v(z))/(2*B_v(z))*d_B_v(z);
} 
double dd_a_v(double z){
    return boost::math::differentiation::finite_difference_derivative(d_a_v,z);
} 
double a(double z){
     boost::math::quadrature::gauss_kronrod<double, 31> integrator;
        auto f_fixed_z= [&](double psi) { return 1/B(z, psi); };
        return sqrt(2.0*integrator.integrate(f_fixed_z, 0,1));
}
double a2(double z){
     boost::math::quadrature::gauss_kronrod<double, 31> integrator;
        auto f_fixed_z= [&](double psi) { return 1/B(z, psi); };
        return 2.0*integrator.integrate(f_fixed_z, 0,1);
}
double d_a(double z){
    return boost::math::differentiation::finite_difference_derivative(a,z);
}
double d_a2(double z){
    return boost::math::differentiation::finite_difference_derivative(a2,z);
}
double p_sred(double z,double psi){
    return p(psi)*(1-B(z,psi)/B_s);
}
double p_avg(double z){
     boost::math::quadrature::gauss_kronrod<double, 31> integrator;
            auto f_fixed_z= [&](double psi) { return p_sred(z,psi)/B(z, psi); };
    double res= 2.0/(a2(z))*integrator.integrate(f_fixed_z,0,1);
    if (res<0
    ){return 0;}else{return res;};
}double r_w(double z){
    return 1.01*a(z);
}
int main(){
    std::vector<std::pair<double,double>> B_VV,B_VV1,B_VV2, pp,pp1,pp2,aa,aa1,aa2,rw,rw1,rw2;
    for(double i=0;i<=1000;i++){
        k=2;
        B_VV.push_back(std::make_pair(i/1000,B_v(i/1000)));
        pp.push_back(std::make_pair(i/1000,p(i/1000)));
        aa.push_back(std::make_pair(i/1000,a(i/1000)));
        rw.push_back(std::make_pair(i/1000,r_w(i/1000)));
        k=4;
        B_VV1.push_back(std::make_pair(i/1000,B_v(i/1000)));
        pp1.push_back(std::make_pair(i/1000,p(i/1000)));
        aa1.push_back(std::make_pair(i/1000,a(i/1000)));
        rw1.push_back(std::make_pair(i/1000,r_w(i/1000)));
        k=608053020;
        B_VV2.push_back(std::make_pair(i/1000,B_v(i/1000)));
        pp2.push_back(std::make_pair(i/1000,p(i/1000)));
        aa2.push_back(std::make_pair(i/1000,a(i/1000)));
        rw2.push_back(std::make_pair(i/1000,r_w(i/1000)));
    }
    Gnuplot gp2;
gp2<<"set terminal wxt 3 title 'a, r_w'\n";
gp2<<"plot '-' with lines title 'k=2', '-' with lines title 'k=4', '-' with lines title 'k=∞','-' with lines title 'k=2', '-' with lines title 'k=4', '-' with lines title 'k=∞'\n";gp2.send1d(pp);gp2.send1d(pp1);gp2.send1d(pp2);//gp2.send1d(rw);gp2.send1d(rw1);gp2.send1d(rw2);

}