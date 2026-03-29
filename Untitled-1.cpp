#include <iostream>
#include <cmath>
#include <math.h>
#include <chrono>
#include <fstream>
#include <filesystem> 
#include <string>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/complex.hpp>
#include "gnuplot-iostream.h"
#include <complex>
#define M_PI 3.14159265358979323846

using namespace std::complex_literals;





double pl=1;   //плотность
double mu=1;


std::complex<double> sigma=1.38*10000*pow(10,12)/1.113; // 3-4 множители перевод в СГС 1.38*10000(Ом*см)^-1
double tau=2.5*pow(10,-14);
double c_0=3*pow(10,10);



double R=1.1;
double B_s=1.0;
double M=8.0;
double q=4.0;
double k=4.0;
double bbbeta=0.9;
double p_0=bbbeta*2*R*R/(2*(bbbeta+R*R-1));
// double p_0=0.6501;   

double a_nenorm_0=1.0;
double a2_nenorm_0=1.0;



double L=350.0;

double RR_w=1.5;
double a_0=20.0;
// double r_0=45;
double r_0=RR_w*a_0;
// double RR_w=r_0/a_0;
double N_c=pow(10,13);
double B_v_=3000.0;
double B_s_=B_v_*R;
double m_i=940*pow(10,6)*1.6*pow(0.1,12); // масса покоя протона в эрг

// double w_0=B_s_/(L*sqrt(N_c*m_i/c_0/c_0));
double w_0=B_v_/(L*sqrt(N_c*1.6726*pow(0.1,24)));

double c=c_0;
double zeta_0=c/w_0/r_0;
//608053020
double tol = 1e-8;
int max_refine = 15;

std::complex<double> dzeta(std::complex<double> omega) {
    using namespace std::complex_literals;

std::complex<double> epsilon=-2.2*pow(10,5)+1.0i*6.55*pow(10,18)/w_0/omega;
//  return ((1.0 - 1.0i)* std::sqrt(omega * mu /(8 * M_PI * sigma)));
//  return std::sqrt(mu /(1.0+1.0i*4.0 * M_PI * sigma/(1.0-1.0i*omega*tau)));
// return sqrt(mu/epsilon);
   return 0.0;
}

double f(double psi){
    if(psi >= 0 && psi <= 1){
        return 1-pow(psi,k);
    }
    else{return 0;
    }
}
double p(double psi){
    return p_0*f(psi)/(2.0*R*R);
}
double B_v(double z){
    return (1.0+(M-1.0)*pow(sin(M_PI*z/2.0),q))/R;
}
double B_v2(double z){
    return 1.0/(B_v(z)*B_v(z));
}
double d_B_v(double z){
    return (M-1.0)*q*M_PI/2.0*cos(M_PI*z/2.0)*pow(sin(M_PI*z/2.0),q-1)/R;
}
double b_v(double z){
    return (1.0+(M-1.0)*pow(sin(M_PI*z/2.0),q))/R;
}
double B(double z,double psi ){
    double bb_v=b_v(z);
    if(bb_v>1.0){return bb_v;}
    else{return sqrt((bb_v*bb_v-p(psi)*2.0)/(1.0-p(psi)*2.0));}
}
double a_v(double z){
    return sqrt(2.0/B_v(z));
} 
double d_a_v(double z){
    return -sqrt(2.0/B_v(z))/(2.0*B_v(z))*d_B_v(z);
} 
double dd_a_v(double z){
    return boost::math::differentiation::finite_difference_derivative(d_a_v,z);
} 
double a_nenorm(double z){
     boost::math::quadrature::gauss_kronrod<double, 31> integrator;
        auto f_fixed_z= [&](double psi) { return 1.0/B(z, psi); };
        return sqrt(2.0*integrator.integrate(f_fixed_z, 0,1));
}
double a(double z){
    return a_nenorm(z)/a_nenorm_0;
}
double a2_nenorm(double z){
     boost::math::quadrature::gauss_kronrod<double, 31> integrator;
        auto f_fixed_z= [&](double psi) { return 1.0/B(z, psi); };
        return 2.0*integrator.integrate(f_fixed_z, 0,1);
}

double a2(double z){
    return a2_nenorm(z)/a2_nenorm_0;
}
double d_a(double z){
    return boost::math::differentiation::finite_difference_derivative(a,z);
}
double d_a2(double z){
    return boost::math::differentiation::finite_difference_derivative(a2,z);
}
double p_sred(double z,double psi){
    return p(psi)*(1.0-B(z,psi)/B_s);
}
double p_avg(double z){
     boost::math::quadrature::gauss_kronrod<double, 31> integrator;
            auto f_fixed_z= [&](double psi) { return p_sred(z,psi)/B(z, psi); };
    double res= 2.0/(a2(z))*integrator.integrate(f_fixed_z,0,1.0);
    if (res<0.0
    ){return 0;}else{return res;};
}
double right_func1(double z){//не используется
    return 2.0*(p_avg(z)/(B_v(z)*B_v(z))*dd_a_v(z)/a_v(z));
}
double right_func2(double z){//не используется
    return (0.5*pow((d_B_v(z)/B_v(z)+d_a2(z)/a2(z)),2)*(1.0-p_avg(z)/(B_v(z)*B_v(z))));
}
double right_func(double z){// в решении уравнения используется только эта функция
    return (2.0*(p_avg(z)/(B_v(z)*B_v(z))*dd_a_v(z)/a_v(z))+0.5*pow((d_B_v(z)/B_v(z)+d_a2(z)/a2(z)),2)*(1.0-p_avg(z)/(B_v(z)*B_v(z))));
}
double w2()//теперь не используетсяы
{
    boost::math::quadrature::gauss_kronrod<double, 81> integrator;
    // return (integrator.integrate(right_func1,-1.0,1.0)+integrator.integrate(right_func2,-1.0,1.0))/integrator.integrate(B_v2,-1.0,1.0)/pl;
    return (integrator.integrate(right_func,-1.0,1.0)/integrator.integrate(B_v2,-1.0,1.0)/pl);

}
double r_w(double z){
    return 1.0;
}
std::string r_w_str = "r_w_1.0000000001_a(z)";
double lambda(double z){
   // return INFINITY;
    return (r_w(z)*r_w(z)*RR_w*RR_w+a(z)*a(z))/(r_w(z)*r_w(z)*RR_w*RR_w-a(z)*a(z));
}
std::complex<double> lambda_0(double z, std::complex<double> omega)// не нормированная, уже не используется
{
    using namespace std::complex_literals;
    return ((omega*r_w(z)/c+1.0i*dzeta(omega))  *r_w(z)*r_w(z)+(omega*r_w(z)/c-1.0i*dzeta(omega))*a2(z))
    /((omega*r_w(z)/c+1.0i*dzeta(omega))*r_w(z)*r_w(z)-(omega*r_w(z)/c-1.0i*dzeta(omega))*a2(z));
}
std::complex<double> lambda_1(double z, std::complex<double> omega)// нормиованная
{
    using namespace std::complex_literals;
    return ((omega*r_w(z)+1.0i*zeta_0*dzeta(omega))*pow(RR_w*r_w(z),2)+(omega*r_w(z)-1.0i*zeta_0*dzeta(omega))*a2(z))
    /((omega*r_w(z)+1.0i*zeta_0*dzeta(omega))*pow(RR_w*r_w(z),2)-(omega*r_w(z)-1.0i*zeta_0*dzeta(omega))*a2(z));
}

double w_2() {// была нужна раньше 

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
std::complex<double> w1=-1.0i*sqrt(0.0137);// эти строки использовались раннее
// Parameters of your ODE: q(x)y'' + p(x)y' + d(x)y + a(x) = 0
// Define your coefficient functions here.
double q1(double z) {
    return lambda(z)+1.0-2.0*p_avg(z)/(B_v(z)*B_v(z));
}
std::complex<double> q2(double z) {
    return lambda_1(z,w1)+1.0-2.0*p_avg(z)/(B_v(z)*B_v(z));
}
double p1(double z){
    auto f = [](double z1){
         return lambda(z1)
               - 2.0*p_avg(z1)/(B_v(z1)*B_v(z1));
    };
    return boost::math::differentiation::finite_difference_derivative(f, z);
}
std::complex<double> p2(double z)
{
    auto f_re = [](double z1) {
        return std::real(lambda_1(z1, w1)
               - 2.0*p_avg(z1)/(B_v(z1)*B_v(z1)));
    };

    auto f_im = [](double z1) {
        return std::imag(lambda_1(z1, w1));
    };

    double dre = boost::math::differentiation::finite_difference_derivative(f_re, z);
    double dim = boost::math::differentiation::finite_difference_derivative(f_im, z);

        return std::complex<double>(dre, dim);
}
double d1(double z) {
    auto f=[](double z1){return -((d_B_v(z1)/B_v(z1)+d_a2(z1)/a2(z1))*(1.0-p_avg(z1)/(B_v(z1)*B_v(z1))));};
    return boost::math::differentiation::finite_difference_derivative(f,z);
}
double a1(double z) {
    return w*pl/(B_v(z)*B_v(z))-right_func(z);
}
std::complex<double> a_2(double z){
        return w1*w1*pl/(B_v(z)*B_v(z))-right_func(z);
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
    return 2.0/M_PI*std::asin(pow((R-1)/(M-1),1/q));
}
double d_phi_z_s(double phi_z_s,double delta){
return (Qfunc(z_s()+delta)-Qfunc(z_s()-delta))/(lambda(z_s())+1)*phi_z_s;
}

void ideal(){//использовалось раннее
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
    double x1 = 1.0;     // End of the interval
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
void my_ode_system1( const state_type1 &y, state_type1 &dydx, double x ) { //решаемое уравнение
    dydx[0] = y[1]; // y0' = y1
    dydx[1] = -( p2(x) * y[1] + d1(x) * y[0] + a_2(x) ) / q2(x); 
}
// std::pair<double,double> in_app(){
// return 0;
// }
bool compareSecond( std::pair<double, double> a,  std::pair<double, double> b) {
    return abs(a.second) < abs(b.second);
}

//начальные приближения
double r=0.11; 

double phase=M_PI/2;

void reshatel(int resuis,double dx) // первый аргумент задаёт название папки для сохранения данных, второй шаг для решения уравнения
{
r_0=a_0*RR_w;
// a_nenorm_0=a_nenorm(0);
// a2_nenorm_0=a2_nenorm(0);
RR_w=RR_w*a(0);
std::cout <<"RR_w: " <<RR_w*a(0)<< std::endl;
    std::cout <<"k: "<< k<<std::endl;
    std::cout <<"M: "<< M<<std::endl;
    std::cout <<"q: "<< q<<std::endl;
    std::cout <<"a(0): "<< a(0)<<std::endl;
    std::cout <<"r_w(0): "<<r_w(0)<<std::endl;
    std::cout <<"w_0: "<< w_0<<std::endl;
    std::cout <<"zeta_0: "<< zeta_0<<std::endl;
    std::cout <<"Lambda : "<< lambda(0)<<std::endl;
    double Z_s=z_s();
        double phi_z_s;
double delta=0.00001; // задаёт точность зануления на правой границе 
double phase_step_min=pow(0.1,4);
double r_step_min=pow(0.1,4);// сделаны чтобы обрывать бесконечные уменьшения шага
double phase_step=-M_PI/30;
double r_step=0.1;
double previous_max=1;
int counter=0;
int counter1=0;
int swich=0; // переключатель для пристрелки по фазе и радиусу
int swich1=0;  // переключатель для повторной пристрелки 

    w1=r*exp(1.0i*phase);
std::vector<std::pair<double,std::pair<double,double>>> bc_arr;
    std::vector<std::pair<double,std::complex<double>>> res,res1,magf,magf_v;
    std::vector<std::pair<double,double>> real_res,im_res,dres,dres1;
        Gnuplot gp;
    double beta=p_0*(R*R-1)/(R*R-p_0);
    
    std::cout<<"Beta: " <<beta<<std::endl;
    std::complex<double> bc=1.0+1.0i;
    std::complex<double> previous_bc=bc*10.0;

    while((std::abs(bc))/abs(previous_max)>delta ){//начало цикла для решения уравнения

std::cout<<std::endl<<"dzeta: " <<dzeta(w1);
        res.clear();dres.clear();dres1.clear();real_res.clear();im_res.clear();
  
        auto my_observer=[&]( const state_type1 &y, const double x ) {// запись решения в данной точке
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

    boost::numeric::odeint::runge_kutta_fehlberg78<state_type1> stepper;// метод решения
    // Initial conditions: {y(x0), y'(x0)}
    state_type1 y0 = { 1, 0};
    // Integration range and initial step size
    double x0 = 0;   // Start of the interval
    double x1 = 1;     // End of the interval
     // Initial step size
    // Perform the integration
    integrate_const( stepper, my_ode_system1, y0, x0, x1, dx, my_observer );// само решение
    bc=dres.back().second+1.0i*dres1.back().second;
    std::cout << "|bc|:"<<std::abs(bc) << "real: " <<std::real(bc)<< " Imag: " <<std::imag(bc)<<" w="<< w1<< ",phase=Pi*"<<phase/M_PI<<"r=" <<std::setprecision(10)<<r<<std::endl;
    std::cout<< "delta="<< std::abs(bc)-std::abs(previous_bc)<<std::endl<<swich1<<std::endl;
    
// Далее сам метод пристрелки(он описан в pdf файлах с результатами) и запись решения в файлы
    if (swich==0){
        
        if(std::abs(previous_bc)<std::abs(bc)){
          //  std::cout << "prev:  "<<std::abs(previous_bc) << "bc:   "<<std::abs(bc); 
            if(abs(phase_step)<phase_step_min){swich=1; counter=0;previous_bc=bc;}else{
            phase-=phase_step;
            phase_step/=2;
            phase+=phase_step;
            counter++;
        if(counter>=5&&(phase_step<phase_step_min*100)){
            phase-=phase_step;
                phase_step=-phase_step*pow(2,counter-3);
                phase+=phase_step;
                counter=0;
            }
        }
        }else {
        phase+=phase_step;
        previous_bc=bc;
       
        previous_max= sqrt(pow(std::max_element(dres.begin(), dres.end(),compareSecond)->second,2)+pow(std::max_element(dres1.begin(), dres1.end(),compareSecond)->second,2));

        counter=0;

    }
    }else{
                previous_max= sqrt(pow(std::max_element(dres.begin(), dres.end(),compareSecond)->second,2)+pow(std::max_element(dres1.begin(), dres1.end(),compareSecond)->second,2));

        if(std::abs(previous_bc)<std::abs(bc)){
            r-=r_step;
            r_step/=2;
            r+=r_step;
            counter++;
            if(counter>=5&&(abs(r_step)<r_step_min*100)){
                r-=r_step;
                r_step=-r_step*pow(2,counter-2);
                r+=r_step;
                counter=0;
            }
            if(abs(r_step)<r_step_min){  if(swich1==5){std::cout <<std::endl << "Dela"<<std::endl; break;}else{r_step_min/=10;phase_step_min/=10;swich1++;swich=0;phase_step*=100;r_step*=100;counter=0;previous_bc=bc;}}
        }else{
        r+=r_step;
        counter=0;
        previous_bc=bc;
        previous_max= sqrt(pow(std::max_element(dres.begin(), dres.end(),compareSecond)->second,2)+pow(std::max_element(dres1.begin(), dres1.end(),compareSecond)->second,2));
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
//     Gnuplot gp1;
//  gp << "set terminal wxt 1  title 'dphi'\n";
// gp << "plot '-' with lines title 'Re(dphi)', '-' with lines title 'Im(dphi)'\n";  gp.send1d(dres);gp.send1d(dres1); 
// gp1 << "set terminal wxt 2  title 'phi'\n";  ;  //gp1 << "set ytics 1\n"; 
// gp1 << "plot '-' with lines title 'Re(phi)', '-' with lines title 'Im(phi)'\n";gp1.send1d(real_res);gp1.send1d(im_res);  
std::cout << std::endl<<"real:"<< std::real( (1.0* std::sqrt(w1 * mu / (8 * M_PI * sigma)) - 1.0i* std::sqrt(w1 * mu / (8 * M_PI * sigma))))<<std::endl<<"Imag:"<<std::imag( (1.0* std::sqrt(w1 * mu / (8 * M_PI * sigma)) - 1.0i* std::sqrt(w1 * mu / (8 * M_PI * sigma)))) ;
std::cout << std::endl<< "w^2="<<w1*w1 << "  w:"<<w1 << "r="<<r<<",phase=Pi*"<<phase/M_PI;
std::cout << std::endl<<"sIgma"<< sigma;

std::vector<std::pair<double,double>> aaa,xi,rw,lambuda;
std::vector<std::pair<double,std::complex<double>>> lambuda1;

double size=1/static_cast<double>(res.size());

double xz=0;
for(int j=0;j<res.size();j++){
   //std::cout << "i:" << j <<"z="<<xz<<std::endl;
    double az=a(xz);
    double bv=B_v(xz);
    double rww=r_w(xz);
    double lambdaa=lambda(xz);
    std::complex<double> lambdaa1=lambda_1(xz,30);
aaa.push_back(std::make_pair(xz,az));
rw.push_back(std::make_pair(xz,rww));
xi.push_back(std::make_pair(xz,(real_res[j].second/bv/az)));
lambuda.push_back(std::make_pair(xz,lambdaa));
lambuda1.push_back(std::make_pair(xz,lambdaa1));
    // std::cout << "a:" << az<<"z:"<<xz<<std::endl;
    // std::cout << "a(0):" << a(0)<<std::endl;
    xz+=size;
}



//std::string name1="_beta="+std::to_string(beta)+"_rw="+r_w_str+"_M="+std::to_string(M)+"_q="+std::to_string(q)+"_k="+std::to_string(k)+"R="+std::to_string(R)+".txt";
std::string name1=".txt";
std::string name_real="Re_phi"+name1;
std::string name_imag="Im_phi"+name1;
std::string name_dreal="Re_dphi"+name1;
std::string name_dimag="Im_dphi"+name1;
std::string name_a="a"+name1;
std::string name_xi="xi"+name1;
std::string name_rw="rw"+name1;
std::string name_w="omega_dzeta"+name1;
std::string name_lambda="lambda"+name1;
std::string name_lambda1="lambda1"+name1;
std::string resis;
if(resuis==0){
     resis="_new_ideal";
}else if(resuis==1){
    
     resis="_resist_+_";
}else if(resuis==2){
     resis="_resist_-_";
}else if(resuis==3){
    
     resis="_resist_-+_";
}else if(resuis==4){
     resis="_resist_--_";

}else if(resuis==5){
     resis="_resist_--_fehlberg_-_fixed_d";

}else if(resuis==6){
     resis="_new_ideal_+";
}
std::string folder1 = "M="+std::to_string(M)+"_"+"q="+std::to_string(q)+"_"+"k="+std::to_string(k)+"_"+"R="+std::to_string(R);
std::string folder2 =r_w_str+resis;
std::string folder3 ="_beta="+std::to_string(beta);
std::filesystem::create_directories(folder1 + "/" + folder2 + "/" + folder3);


{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_real);
    for (const auto& [x, y] : real_res) {
        file << x << " " << y << '\n';
    }
    file.close();}
{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_imag);
    for (const auto& [x, y] : im_res) {
        file << x << " " << y << '\n';
    }
    file.close();}
{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_dreal);
    for (const auto& [x, y] : dres) {
        file << x << " " << y << '\n';
    }
    file.close();}
{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_dimag);
    for (const auto& [x, y] : dres1) {
        file << x << " " << y << '\n';
    }
    file.close();}
{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_a);
    for (const auto& [x, y] : aaa) {
        file << x << " " << y << '\n';
    }
    file.close();}
{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_xi);
    for (const auto& [x, y] : xi) {
        file << x << " " << y << '\n';
    }
    file.close();}
{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_rw);
    for (const auto& [x, y] : rw) {
        file << x << " " << y << '\n';
    }
    file.close();}
{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_w);
        file << "omega:" <<w1 << " dzeta:"<< dzeta(w1) << " bc:" << sqrt(pow(dres.back().second,2)+pow(dres1.back().second,2));
    file.close();}

{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_lambda);
    for (const auto& [x, y] : lambuda) {
        file << x << " " << y << '\n';
    }
    file.close();}
{std::ofstream file(folder1 + "/" + folder2 + "/" +folder3 + "/" + name_lambda1);
    for (const auto& [x, y] : lambuda1) {
        file << x << " " << y << '\n';
    }
    file.close();}



// Gnuplot gp2;
// gp2<<"set terminal wxt 3 title 'a'\n";
// gp2<<"plot '-' with lines title 'a', '-' with lines title 'xi', '-' with lines title 'r_w'\n";gp2.send1d(aaa);gp2.send1d(xi);gp2.send1d(rw);
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

// решение при заданной w
void reshatel_0(){

double Z_s=z_s();
double phi_z_s;
double delta=0.0001;
// double r=1.1400000000000;
// double phase=M_PI/2*0;
double phase_step=-M_PI/3;
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
    bc=dres.back().second+1.0i*dres1.back().second;
    std::cout << "|bc|:"<<std::abs(bc) << "real: " <<std::real(bc)<< " Imag: " <<std::imag(bc)<<" w="<< w1<< ",phase=Pi*"<<phase/M_PI<<"r=" <<std::setprecision(10)<<r<<std::endl;
    std::cout<< "delta="<< std::abs(bc)-std::abs(previous_bc)<<std::endl;
  
    
    std::cout << std::endl <<"R_step: "<< r_step <<" Phase_step: "<< phase_step << " Prev: "<< std::abs(previous_bc);
    
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
std::cout << std::endl<<"sIgma"<< sigma;



// std::vector<std::pair<double,double>> aaa,xi,rw;

// double size=1/static_cast<double>(res.size());

// double xz=-size;
// for(int j=0;j<res.size();j++){
//    //std::cout << "i:" << j <<"z="<<xz<<std::endl;
//     xz+=size;
//     double az=a(xz);
//     double bv=B_v(xz);
//     double rww=r_w(xz);
// aaa.push_back(std::make_pair(xz,az));
// rw.push_back(std::make_pair(xz,rww));
// xi.push_back(std::make_pair(xz,(real_res[j].second/bv/az)));
// }

// Gnuplot gp2;
// gp2<<"set terminal wxt 3 title 'a'\n";
// gp2<<"plot '-' with lines title 'a', '-' with lines title 'xi', '-' with lines title 'r_w'\n";gp2.send1d(aaa);gp2.send1d(xi);gp2.send1d(rw);
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


int main(){

r=0.95; 

 phase=M_PI/2;
 r_w_str="st_RR_w"+std::to_string(RR_w);
 
 for(double i=0.250;i<=0.99;i+=0.025){// в виде подобного цикла обычно считаются графики, которые я вам показываю
    RR_w=1.5;
    bbbeta=i;   
    p_0=bbbeta*2.0*R*R/(2.0*(bbbeta+R*R-1.0));

    reshatel(0, 0.003);

}
}