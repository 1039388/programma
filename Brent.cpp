#include <iostream>
#include <cmath>
#include <math.h>
#include <chrono>
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

double pl=1;//плотность
double mu=1;
double sigma=1.38*10000;    
double c=3*pow(10,10);

double R=1.5;
double B_s=1;
double M=4;
double q=8;
double bbbeta=0.64;
double p_0=bbbeta*2*R*R/(2*(bbbeta+R*R-1));
// double p_0=0.6501;   
double k=2;
// double dzeta=;
//608053020
double tol = 1e-8;
int max_refine = 15;
std::complex<double> dzeta(std::complex<double> omega) {
    using namespace std::complex_literals;
  return (1.0* std::sqrt(omega * mu / (8 * M_PI * sigma)) - 1.0i* std::sqrt(omega * mu / (8 * M_PI * sigma)));
//   return 0.0;
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
}
double right_func1(double z){
    return 2.0*(p_avg(z)/(B_v(z)*B_v(z))*dd_a_v(z)/a_v(z));
}
double right_func2(double z){
    return (0.5*pow((d_B_v(z)/B_v(z)+d_a2(z)/a2(z)),2)*(1.0-p_avg(z)/(B_v(z)*B_v(z))));
}
double right_func(double z){
    return (2.0*(p_avg(z)/(B_v(z)*B_v(z))*dd_a_v(z)/a_v(z))+0.5*pow((d_B_v(z)/B_v(z)+d_a2(z)/a2(z)),2)*(1.0-p_avg(z)/(B_v(z)*B_v(z))));
}
double w2()
{
    boost::math::quadrature::gauss_kronrod<double, 81> integrator;
    // return (integrator.integrate(right_func1,-1.0,1.0)+integrator.integrate(right_func2,-1.0,1.0))/integrator.integrate(B_v2,-1.0,1.0)/pl;
    return (integrator.integrate(right_func,-1.0,1.0)/integrator.integrate(B_v2,-1.0,1.0)/pl);

}
double r_w(double z){
    return a(z)*1.2;
}
double lambda(double z){
   // return INFINITY;
    return (r_w(z)*r_w(z)+a(z)*a(z))/(r_w(z)*r_w(z)-a(z)*a(z));
}
std::complex<double> lambda_1(double z, std::complex<double> omega){
    using namespace std::complex_literals;
    return ((omega*r_w(z)/c+1.0i*dzeta(omega))*r_w(z)*r_w(z)+(omega*r_w(z)/c-1.0i*dzeta(omega))*a2(z))/((omega*r_w(z)/c+1.0i*dzeta(omega))*r_w(z)*r_w(z)-(omega*r_w(z)/c-1.0i*dzeta(omega))*a2(z));
}
double w_2() {

    Gnuplot gp; // Creates a Gnuplot instance
   std::vector<std::pair<double, double>> w_b,w_b1,w_b2,w_b3,nol_b;
   std::vector<std::pair<double, double>> plot_data1,plot_data2,plot_data3,plot_data4,plot_data5,plot_data6;
   int nomer=0;
   while(p_0<0.9999){

    //  gp << "set terminal qt 1 title 'B(z,1)'\n";    gp << "plot '-' with lines\n";    gp.send1d(plot_data2);
    // gp << "set terminal qt 2 title 'a_v(z)'\n";    gp << "plot '-' with lines\n";    gp.send1d(plot_data3);
    // gp << "set terminal qt 0 title 'a(z)'\n"; gp << "plot '-' with lines\n";  gp.send1d(plot_data4);
    // gp << "set terminal qt 4 title 'B_v'\n";    gp << "plot '-' with lines\n";    gp.send1d(plot_data5);
  //  gp << "set terminal qt 5 title 'OTLADKA'\n";    gp << "plot '-' with lines\n";    gp.send1d(plot_data6);
   //gp << "pause -1\n"; 
    // Keep windows open
    k=1;
    double w=w2();
    double beta=p_0*(R*R-1)/(R*R-p_0) ;
//   std::cout << std::endl<<"w^2="<< w << std::endl;
//   std::cout << std::endl<< "beta=" << beta << std::endl;
    w_b.push_back({beta,w});

    k=2;
     w=w2();
    w_b1.push_back({beta,w});
    k=4;
     w=w2();
    w_b2.push_back({beta,w});
    k=608053020;
     w=w2();
    w_b3.push_back({beta,w});

    nol_b.push_back({beta,0});

    nomer++;
    p_0=p_0+0.01;
    std::cout<<beta<<std::endl;
}
 for (double z = -1.1; z <= 1.11; z += 0.01) {
     //   plot_data1.push_back({z,right_func(x)});
        //  plot_data2.push_back({z,B(z,0)});
        //  plot_data3.push_back({z,a_v(z)});
        // plot_data4.push_back({z,a(z)});
        plot_data5.push_back({z,B_v(z)});
        plot_data4.push_back({z,B(z,0)});
      //   plot_data6.push_back({x,1/2*pow((d_B_v(x)/B_v(x)+d_a2(x)/a2(x)),2)*(1-p_avg(x)/(B_v(x)*B_v(x)))});
        // plot_data6.push_back({z,right_func2(z)+right_func1(z)});
 }
//gp << "set terminal qt\n";

//   gp << "set terminal wxt 3 title 'B_v(z)'\n"; gp << "plot '-' with lines\n"; gp.send1d(plot_data5);
 //gp << "set terminal wxt 1 title 'B'\n"; gp << "plot '-' with lines\n";  gp.send1d(plot_data4);

   gp << "set terminal wxt 0 title 'w^2(β)'\n";
    gp << "set xlabel 'beta'\n";
    gp << "set ylabel 'w^2'\n";
  //gp << "set yrange [-0.1:0.1]\n"; 
    //gp << "set xrange [0:1]\n";         // Set the Y-axis limit
gp << "set xtics 0.1\n";           // Set major tic spacing on X-axis
gp << "set mxtics 5\n";            // Set number of minor tic intervals on X-axis
gp << "set grid xtics mxtics\n";   // Enable grid for major and minor X-tics
gp << "set key outside right top\n";  // Position legend outside the plot
gp << "set border linewidth 1.5\n";       // Жирная рамка
gp << "set style line 1 linewidth 2\n"; // Стиль для первой линии
gp << "set style line 2 linewidth 5\n"; // Стиль для второй линии
//gp<<"set border lw 2\n"; gp <<"set linestyle lw 2 "; 
gp << "plot '-' with lines linewidth 2 title 'k=1', '-' with lines linewidth 2 title 'k=2', '-' with lines linewidth 2 title 'k=4', '-' with lines linewidth 2 title 'k=inf', 0 with lines linewidth 2 title '0'\n";
gp.send1d(w_b);
gp.send1d(w_b1);
gp.send1d(w_b2);
gp.send1d(w_b3);
return 0;
}
// Define the state type: a vector of two elements (y0, y1)

typedef std::array<double, 2> state_type;
double w=1;
std::complex<double> w1=-1.0i*sqrt(0.0137);
// Parameters of your ODE: q(x)y'' + p(x)y' + d(x)y + a(x) = 0
// Define your coefficient functions here.
double q1(double z) {
    return lambda(z)+1.0-2.0*p_avg(z)/(B_v(z)*B_v(z));
}
std::complex<double> q2(double z) {
    return lambda_1(z,w1)+1.0-2.0*p_avg(z)/(B_v(z)*B_v(z));
}
double p1(double z) {
    auto f=[](double z1){return lambda(z1)-2.0*p_avg(z1)/(B_v(z1)*B_v(z1));};
    return boost::math::differentiation::finite_difference_derivative(f,z);
}
double d1(double z) {
    auto f=[](double z1){return -((d_B_v(z1)/B_v(z1)+d_a2(z1)/a2(z1))*(1.0-p_avg(z1)/(B_v(z1)*B_v(z1))));};
    return boost::math::differentiation::finite_difference_derivative(f,z);
}
double a1(double z) {
    return w*pl/(B_v(z)*B_v(z)  )-right_func(z);
}
std::complex<double> a_2(double z){
        return w1*w1*pl/(B_v(z)*B_v(z)  )-right_func(z);
}
// Define the ODE system (the RHS of our system of equations)
void my_ode_system( const state_type &y, state_type &dydx, double x ) {
    dydx[0] = y[1]; // y0' = y1
    dydx[1] = -( p1(x) * y[1] + d1(x) * y[0] + a1(x) ) / q1(x); // y1' = -(...)/q(x)
}
double Qfunc(double z){
    return d_B_v(z)/B_v(z)+d_a2(z)/a2(z);
}
double z_s(){
    return 2/M_PI*std::asin(pow((R-1)/(M-1),1/q));
}
double d_phi_z_s(double phi_z_s,double delta){
return (Qfunc(z_s()+delta)-Qfunc(z_s()-delta))/(lambda(z_s())+1)*phi_z_s;
}

// An observer to print the state at each step
void ideal(){
    double Z_s=z_s();
    double phi_z_s;
    std::vector<std::pair<double,double>> res,res1,dres,dres1,magf,magf_v;
        Gnuplot gp;

    double delta=0.001/100;
    double beta=p_0*(R*R-1)/(R*R-p_0);
    std::cout<< beta<<std::endl;
    auto my_observer=[&]( const state_type &y, const double x ) {
    res.push_back(std::make_pair(x,y[0]));
    dres.push_back(std::make_pair(x,y[1]));
    magf.push_back(std::make_pair(x,B(x,0.2)));
    magf_v.push_back(std::make_pair(x,B_v(x)));
    if(Z_s<res.back().first&&Z_s<res[res.size()-2].first){
        phi_z_s=y[0];
    }
    };
    auto my_observer1=[&]( const state_type &y, const double x ) {
    res1.push_back(std::make_pair(x,y[0]));
    dres1.push_back(std::make_pair(x,y[1]));
    };
    
    dres.push_back(std::make_pair(1.0,1.0));

    while(fabs(dres.back().second)>delta){
    w+=0.01/50;
   res.clear();
   res1.clear();
   dres.clear();
   dres1.clear();
    boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;
    // Initial conditions: {y(x0), y'(x0)}
    state_type y0 = { 1, 0};
    // Integration range and initial step size
    double x0 = 0;   // Start of the interval
    double x1 = 1;     // End of the interval
    double dx = 0.001; // Initial step size
    // Perform the integration
    integrate_const( stepper, my_ode_system, y0, x0, x1, dx, my_observer );
       y0 = { 1, 0};
      dx=-dx;     // End of the interval
    integrate_const( stepper, my_ode_system, y0, x0, -1.0, -dx, my_observer1 );
    std::cout<<"w^2: " <<w<< " phi'="<<dres.back().second <<std::endl;
}
std::cout<<std::endl<<"beta: " <<beta<<std::endl;
std::cout<<"w^2: " <<w<< std::endl;
std::cout<<"z_s: " <<Z_s<< std::endl;
std::cout<<"phi(z_s): " <<phi_z_s<< std::endl;
std::cout<<"[d_phi]|z_s: " <<d_phi_z_s(phi_z_s,0.05)<< std::endl;

    
gp << "set terminal wxt 1  title 'y'\n";  ; 
 gp << "plot '-' with lines, '-' with lines\n";  gp.send1d(res); gp.send1d(res1);

   gp << "set terminal wxt 2  title 'dy'\n";  ; 
gp << "plot '-' with lines\n";  gp.send1d(dres);
//    gp << "set terminal wxt 3  title 'B'\n";  ; 
// gp << "plot '-' with lines\n";  gp.send1d(magf);
//    gp << "set terminal wxt 4  title 'B'\n";  ; 
// gp << "plot '-' with lines\n";  gp.send1d(magf_v);
}

std::complex<double> q11(double z) {
    return 1;
}
double p11(double z) {
    return 0;
}
double d11(double z) {
    return 1;
}
double a11(double z) {
    return 0;
}
typedef std::vector<std::complex<double>> state_type1;
// Define the ODE system (the RHS of our system of equations)
void my_ode_system1( const state_type1 &y, state_type1 &dydx, double x ) {
    dydx[0] = y[1]; // y0' = y1
    dydx[1] = -( p1(x) * y[1] + d1(x) * y[0] + a_2(x) ) / q2(x); 
}
// std::pair<double,double> in_app(){
// return 0;
// }
bool compareSecond( std::pair<double, double> a,  std::pair<double, double> b) {
    return abs(a.second) < abs(b.second);
}
int main(){

    double Z_s=z_s();
        double phi_z_s;
double delta=0.00001;
double r=1.00000000000;
double phase=M_PI/4;
double phase_step=M_PI/5;
double r_step=0.1;
double previous_max=1;
int counter=0;
int swich=0;
int swich1=0;

    w1=r*exp(1.0i*phase);
std::vector<std::pair<double,std::pair<double,double>>> bc_arr;
    std::vector<std::pair<double,std::complex<double>>> res,res1,magf,magf_v;
    std::vector<std::pair<double,double>> real_res,im_res,dres,dres1;
        Gnuplot gp;
    double beta=p_0*(R*R-1)/(R*R-p_0);
    
    std::cout<< beta<<std::endl;
    std::complex<double> bc=1.0+1.0i;
    std::complex<double> previous_bc=bc*10.0;

    while((std::abs(bc))/abs(previous_max)>delta ){

std::cout<<std::endl<<"dzeta: " <<dzeta(w1);
        res.clear();dres.clear();dres1.clear();real_res.clear();im_res.clear();
    auto my_observer=[&]( const state_type1 &y, const double x ) {
    res.push_back(std::make_pair(x,y[0]));
    dres.push_back(std::make_pair(x,std::real(y[1])));
    dres1.push_back(std::make_pair(x,std::imag(y[1])));

    real_res.push_back(std::make_pair(x,std::real(y[0])));
    im_res.push_back(std::make_pair(x,std::imag(y[0])));
     if(Z_s<res.back().first&&Z_s<res[res.size()-2].first){
        phi_z_s=std::real(y[0]);
    }
    };
    //dres.push_back(std::make_pair(1.0,1.0));

    boost::numeric::odeint::runge_kutta_dopri5<state_type1> stepper;
    // Initial conditions: {y(x0), y'(x0)}
    state_type1 y0 = { 1, 0};
    // Integration range and initial step size
    double x0 = 0;   // Start of the interval
    double x1 = 1;     // End of the interval
    double dx = 0.01; // Initial step size
    // Perform the integration
    integrate_const( stepper, my_ode_system1, y0, x0, x1, dx, my_observer );
    bc=real_res.back().second+1.0i*im_res.back().second;
    std::cout << "|bc|:"<<std::abs(bc) << "real: " <<std::real(bc)<< " Imag: " <<std::imag(bc)<<" w="<< w1<< ",phase=Pi*"<<phase/M_PI<<"r=" <<std::setprecision(10)<<r<<std::endl;
    std::cout<< "delta="<< std::abs(bc)-std::abs(previous_bc)<<std::endl;
    

    if (swich==0){
        
        if(std::abs(previous_bc)<std::abs(bc)){
          //  std::cout << "prev:  "<<std::abs(previous_bc) << "bc:   "<<std::abs(bc); 
            if(abs(phase_step)<pow(0.1,7)){swich=1; counter=0;previous_bc=bc;}else{
            phase-=phase_step;
            phase_step/=3;
            phase+=phase_step;
            counter++;
        if(counter>=5&&(phase_step)<1e-4){
            phase-=phase_step;
                phase_step=-phase_step*pow(2,counter-1);
                phase+=phase_step;
                counter=0;
            }
        }
        }else{
        phase+=phase_step;
        previous_bc=bc;
       
        previous_max= sqrt(pow(std::max_element(real_res.begin(), real_res.end(),compareSecond)->second,2)+pow(std::max_element(im_res.begin(), im_res.end(),compareSecond)->second,2));

        counter=0;

    }
    }else{
        previous_max= sqrt(pow(std::max_element(real_res.begin(), real_res.end(),compareSecond)->second,2)+pow(std::max_element(im_res.begin(), im_res.end(),compareSecond)->second,2));

        if(std::abs(previous_bc)<std::abs(bc)){
            r-=r_step;
            r_step/=2;
            r+=r_step;
            counter++;
            if(counter>=5&&(abs(r_step)<1e-4)){
                r-=r_step;
                r_step=-r_step*pow(2,counter-1);
                r+=r_step;
                counter=0;
            }
            if(abs(r_step)<1e-8){  std::cout <<std::endl << "Dela"<<std::endl; break;}
        }else{
        r+=r_step;
        counter=0;
        previous_bc=bc;
        previous_max= sqrt(pow(std::max_element(real_res.begin(), real_res.end(),compareSecond)->second,2)+pow(std::max_element(im_res.begin(), im_res.end(),compareSecond)->second,2));
    }
    }

    w1=r*exp(1.0i*phase);
    
    std::cout << std::endl <<"R_step: "<< r_step <<" Phase_step: "<< phase_step << " Prev: "<< std::abs(previous_bc);
    }
    std::cout<<"w^2: " <<w1*w1<< " phi'="<<dres.back().second <<std::endl;

std::cout<<std::endl<<"beta: " <<beta<<std::endl;
std::cout<<std::endl<<"dzeta: " <<dzeta(w1)<<std::endl;
std::cout<<"z_s: " <<Z_s<< std::endl;
 std::cout<<"phi(z_s): " <<phi_z_s<< std::endl;
 std::cout<<"[d_phi]|z_s: " <<d_phi_z_s(phi_z_s,0.05)<< std::endl;
std::cout <<"|d_phi|1|:"<< sqrt(pow(dres.back().second,2)+pow(dres1.back().second,2));
    Gnuplot gp1;
 gp << "set terminal wxt 1  title 'dphi'\n";
gp << "plot '-' with lines title 'Re(dphi)', '-' with lines title 'Im(dphi)'\n";  gp.send1d(dres);gp.send1d(dres1); 
gp1 << "set terminal wxt 2  title 'phi'\n";  ;  //gp1 << "set ytics 1\n"; 
gp1 << "plot '-' with lines title 'Re(phi)', '-' with lines title 'Im(phi)'\n";gp1.send1d(real_res);gp1.send1d(im_res);  
std::cout << std::endl<<"real:"<< std::real( (1.0* std::sqrt(w1 * mu / (8 * M_PI * sigma)) - 1.0i* std::sqrt(w1 * mu / (8 * M_PI * sigma))))<<std::endl<<"Imag:"<<std::imag( (1.0* std::sqrt(w1 * mu / (8 * M_PI * sigma)) - 1.0i* std::sqrt(w1 * mu / (8 * M_PI * sigma)))) ;
std::cout << std::endl<< "w^2="<<w1*w1 << "  w:"<<w1 << "r="<<r<<",phase=Pi*"<<phase/M_PI;
std::vector<std::pair<double,double>> aaa,xi,rw;

double size=1/static_cast<double>(res.size());

double xz=-size;
for(int j=0;j<res.size();j++){
   //std::cout << "i:" << j <<"z="<<xz<<std::endl;
    xz+=size;
    double az=a(xz);
    double bv=B_v(xz);
    double rww=r_w(xz);
aaa.push_back(std::make_pair(xz,az));
rw.push_back(std::make_pair(xz,rww));
xi.push_back(std::make_pair(xz,(real_res[j].second/bv/az)));
}
Gnuplot gp2;
gp2<<"set terminal wxt 3 title 'a'\n";
gp2<<"plot '-' with lines title 'a', '-' with lines title 'xi', '-' with lines title 'r_w'\n";gp2.send1d(aaa);gp2.send1d(xi);gp2.send1d(rw);
// std::vector<std::pair<double,double>> lab1,lab2,lab21;
// for(double i=0;i<1000;i++){
//     lab1.push_back(std::make_pair(i/1000.0,lambda(i/1000.0)));
//     lab2.push_back(std::make_pair(i/1000.0,std::real(lambda_1(i/1000.0,w1))));
//     lab21.push_back(std::make_pair(i/1000.0,std::imag(lambda_1(i/1000.0,w1))));
// }
//  gp << "set terminal wxt 3  title 'lamba'\n";
// gp << "plot '-' with lines title 'lambda', '-' with lines title 'Re(lam)', '-' with lines title 'Im(lam)'\n";  gp.send1d(lab1);gp.send1d(lab2); gp.send1d(lab21); 
// std::cout <<std::endl <<lambda(0.5)<<std::endl <<lambda(0.7) << std::endl << std::real(lambda_1(0.5,0.0001)) << std::endl << std::imag(lambda_1(0.5,w1));
}
