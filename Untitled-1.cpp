#include <iostream>
#include <cmath>
#include <math.h>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <filesystem> 
#include <string>
#include <stdexcept>
#include <boost/math/constants/constants.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/exp_sinh.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/differentiation/finite_difference.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/math/interpolators/barycentric_rational.hpp>
// #include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/interpolators/pchip.hpp>
#include <boost/math/interpolators/makima.hpp>
#include <boost/math/complex.hpp>
#include "gnuplot-iostream.h"
#include <complex>
#include <execution>
#include <memory> 
#include <algorithm> 
#include <vector>
#define M_PI 3.14159265358979323846
using namespace std::complex_literals;
typedef std::vector<std::complex<double>> state_type1;

class FileInterpolator {
private:
    std::unique_ptr<boost::math::interpolators::barycentric_rational<double>> interp;
    // using interpolator_type = boost::math::interpolators::makima<std::vector<double>>;
    // std::unique_ptr<interpolator_type> interp;
    bool ready = false;
    double z_min = 0.0;
    double z_max = 1.438794035387821;

public:
    FileInterpolator() = default;

    // Метод загрузки данных из файла
    void load(const std::string& filename) {
        std::vector<double> x, y;
        std::ifstream file(filename);
        double xi, yi;
        
        if (!file.is_open()) {
            std::cerr << "open error " << filename << std::endl;
            return;
        }
double z;
        while (file >> xi >> yi >> z) {
            x.push_back(xi);
            y.push_back(yi);
        }
        file.close();

        if (x.size() > 2) {
            // Создаем интерполятор (порядок аппроксимации 3)
            interp = std::make_unique<boost::math::interpolators::barycentric_rational<double>>(
                std::move(x), std::move(y), 3
            );
            ready = true;
        }
        // if (x.size() >= 5) {
        //     // В конструктор Makima передаются только векторы x и y.
        //     interp = std::make_unique<interpolator_type>(std::move(x), std::move(y));
        //     ready = true;
        // } else if (x.size() > 0) {
        //     std::cerr << "Error: Not enough points for Makima (need at least 5)" << std::endl;
        // }
    }
    

    // Получить значение функции в точке z
    double val(double z) const {
        if (z < z_min) z = z_min;
    if (z > z_max) z = z_max;

        return ready ? (*interp)(z) : 0.0;
    }

    // Получить точную аналитическую производную в точке z
    double prime(double z) const {
        return ready ? interp->prime(z) : 0.0;
    }

    // Проверка, загружен ли файл
    bool is_ready() const { return ready; }
};

struct SharedInterpolationData {
    FileInterpolator Bv_data;
    FileInterpolator dBv_data;
    FileInterpolator ddav_av_data;
    FileInterpolator dBv_Bv_data;

    // аргументы со значениями по умолчанию, можно менять указывая их в аргументах или не менять, вовсе ничего не указывая
    void load_all(const std::string& papka_int= R"(C:\Users\MJ\Desktop\ballon\pole\clean\)", 
                  const std::string& B_v = "tblbv_1_2.dat",
                  const std::string& dB_v = "tblbvD1_1_2.dat",
                  const std::string& dB_v_Bv = "tblAvD1Av_1_2.dat",
                  const std::string& ddav_av = "tblAvD2Av_1_2.dat") 
    {
        Bv_data.load(papka_int + B_v);
        dBv_data.load(papka_int + dB_v);
        ddav_av_data.load(papka_int + ddav_av);
        dBv_Bv_data.load(papka_int + dB_v_Bv);
        
        std::cout << "All interpolation data loaded from " << papka_int << " successfully." << std::endl;
    }
};


struct PlasmaModel_A1{
double B_s_, c, zeta_0, r_0;
double p_0, w_0;

double pl=1;   //плотность
double mu=1;
std::complex<double> sigma=1.38*10000*pow(10,12)/1.113; // 3-4 множители перевод в СГС 1.38*10000(Ом*см)^-1
double tau=2.5*pow(10,-14);
double c_0=3*pow(10,10);
double m_i=940*pow(10,6)*1.6*pow(0.1,12); // масса покоя протона в эрг

double L=350.0;
double N_c=pow(10,13);
double B_v_=3000.0;


double RR_w=2.25;
double a_0=20.0;

double R=1.1;
double B_s=1.0;
double M=8.0;
double q=4.0;
double k=4.0;
double bbbeta=0.9;



// double w_0=B_s_/(L*sqrt(N_c*m_i/c_0/c_0));

//608053020
PlasmaModel_A1() { update(); }

void update() {
    p_0 = bbbeta * 2.0 * R * R / (2.0 * (bbbeta + R * R - 1.0));
    B_s_=B_v_*R;
    w_0=B_v_/(L*sqrt(N_c*1.6726*pow(0.1,24)));
    r_0=RR_w*a_0;
    c=c_0;
    zeta_0=c/w_0/r_0;
    last_z_ddav = last_z_an = last_z_a2n = -1.0;
    last_z = last_z_da2 = last_z_pavg = -1.0;
        // Загрузка файлов и инициализация интерполяторов здесь (один раз!)
        // load_interpolation_data();
}


std::complex<double> dzeta(std::complex<double> omega) const{
    using namespace std::complex_literals;

std::complex<double> epsilon=-2.2*pow(10,5)+1.0i*6.55*pow(10,18)/w_0/omega;
//  return ((1.0 - 1.0i)* std::sqrt(omega * mu /(8 * M_PI * sigma)));
//  return std::sqrt(mu /(1.0+1.0i*4.0 * M_PI * sigma/(1.0-1.0i*omega*tau)));
// return sqrt(mu/epsilon);
   return 0.0;
}

double f (double psi)const{
    if(psi >= 0 && psi <= 1){
        return 1-pow(psi,k);
    }
    else{return 0;
    }
}
double p (double psi)const{
    return p_0*f(psi)/(2.0*R*R);
}

double B_v (double z)const{
    return (1.0+(M-1.0)*pow(sin(M_PI*z/2.0),q))/R;
    // return Bv_int(z);
}
double B_v2 (double z)const{
    return 1.0/(B_v(z)*B_v(z));
    // return 1.0/(Bv_int(z)*Bv_int(z));
}
double d_B_v (double z)const{
    return (M-1.0)*q*M_PI/2.0*cos(M_PI*z/2.0)*pow(sin(M_PI*z/2.0),q-1)/R;
    // return DBv_int(z);
}
double b_v (double z)const{
    return (1.0+(M-1.0)*pow(sin(M_PI*z/2.0),q))/R;
    // return Bv_int(z);

}
double B (double z,double psi )const{
    double bb_v=b_v(z);

    if(bb_v>1.0){return bb_v;}
    else{double pp=p(psi)*2.0;
        return sqrt((bb_v*bb_v-pp)/(1.0-pp));}
}
double a_v (double z)const{
    return sqrt(2.0/B_v(z));
} 
double d_a_v (double z)const{
    return -sqrt(2.0/B_v(z))/(2.0*B_v(z))*d_B_v(z);
} 
mutable double last_z_ddav = -1.0, last_ddav_val = 0.0;
double dd_a_v(double z) const {
        if (std::abs(z - last_z_ddav) < 1e-12) return last_ddav_val;
        auto func = [this](double x) { return this->d_a_v(x); };
        last_ddav_val = boost::math::differentiation::finite_difference_derivative(func, z);
        last_z_ddav = z;
        return last_ddav_val;
    }
mutable double last_z_an = -1.0;
mutable double last_an_val = 0.0;
double a_nenorm(double z) const {
    if (std::abs(z - last_z_an) < 1e-12) return last_an_val;

    boost::math::quadrature::gauss_kronrod<double, 31> integrator;
    auto f_fixed_z = [this, z](double psi) { return 1.0 / B(z, psi); };
    last_an_val = sqrt(2.0 * integrator.integrate(f_fixed_z, 0, 1));
    last_z_an = z;
    return last_an_val;
}
double a (double z)const{
    return a_nenorm(z);
}
mutable double last_z_a2n = -1.0, last_a2n_val = 0.0;
double a2_nenorm(double z) const {
        if (std::abs(z - last_z_a2n) < 1e-12) return last_a2n_val;
        boost::math::quadrature::gauss_kronrod<double, 31> integrator;
        auto f_fixed_z = [this, z](double psi) { return 1.0 / B(z, psi); };
        last_a2n_val = 2.0 * integrator.integrate(f_fixed_z, 0, 1);
        last_z_a2n = z;
        return last_a2n_val;
    }
double a2 (double z)const{
    return a2_nenorm(z);
}
mutable double last_z = -1.0;
mutable double last_da = 0.0;
double d_a (double z)const{
    if (std::abs(z - last_z) < 1e-12) return last_da; // Возвращаем кэш
    auto func = [this](double x) { return this->a(x); };
    last_z = z;
    last_da = boost::math::differentiation::finite_difference_derivative(func, z);
    return last_da;
}

mutable double last_z_da2 = -1.0;
mutable double last_da2_val = 0.0;
double d_a2(double z) const {
    if (std::abs(z - last_z_da2) < 1e-12) return last_da2_val; 
    
    auto func = [this](double x) { return this->a2(x); };
    last_da2_val = boost::math::differentiation::finite_difference_derivative(func, z);
    last_z_da2 = z;
    return last_da2_val;
}
double p_sred (double z,double psi)const{
    return p(psi)*(1.0-B(z,psi)/B_s);
}
mutable double last_z_pavg = -1.0, last_pavg_val = 0.0;
double p_avg(double z) const {
        if (std::abs(z - last_z_pavg) < 1e-12) return last_pavg_val;
        boost::math::quadrature::gauss_kronrod<double, 31> integrator;
        auto f_fixed_z = [this, z](double psi) { return p_sred(z, psi) / B(z, psi); };
        double res = 2.0 / (a2(z)) * integrator.integrate(f_fixed_z, 0, 1.0);
        last_pavg_val = (res < 0.0) ? 0.0 : res;
        last_z_pavg = z;
        return last_pavg_val;
    }
};
struct PlasmaModel_A1_int{
double B_s_, c, zeta_0, r_0;
double p_0, w_0;

double pl=1;   //плотность
double mu=1;
std::complex<double> sigma=1.38*10000*pow(10,12)/1.113; // 3-4 множители перевод в СГС 1.38*10000(Ом*см)^-1
double tau=2.5*pow(10,-14);
double c_0=3*pow(10,10);
double m_i=940*pow(10,6)*1.6*pow(0.1,12); // масса покоя протона в эрг

double L=350.0;
double N_c=pow(10,13);
double B_v_=3000.0;


double RR_w=2.25;
double a_0=20.0;

double R=1.1;
double B_s=1.0;
double M=8.0;
double q=4.0;
double k=4.0;
double bbbeta=0.9;



// double w_0=B_s_/(L*sqrt(N_c*m_i/c_0/c_0));

//608053020

// Ссылки на общие данные, чтобы не читать файлы на каждом потоке
    const FileInterpolator& Bv_data;
    const FileInterpolator& dBv_data;
    const FileInterpolator& ddav_av_data;
    const FileInterpolator& dBv_Bv_data;

    PlasmaModel_A1_int(const SharedInterpolationData& shared_data) 
        : Bv_data(shared_data.Bv_data), 
          dBv_data(shared_data.dBv_data), 
          ddav_av_data(shared_data.ddav_av_data), 
          dBv_Bv_data(shared_data.dBv_Bv_data) 
    { 
        update(); 
    }

void update() {
    p_0 = bbbeta * 2.0 * R * R / (2.0 * (bbbeta + R * R - 1.0));
    B_s_=B_v_*R;
    w_0=B_v_/(L*sqrt(N_c*1.6726*pow(0.1,24)));
    r_0=RR_w*a_0;
    c=c_0;
    zeta_0=c/w_0/r_0;
    last_z_ddav = last_z_an = last_z_a2n = -1.0;
    last_z = last_z_da2 = last_z_pavg = -1.0;
        // load_interpolation_data();
}


std::complex<double> dzeta(std::complex<double> omega) const{
    using namespace std::complex_literals;

std::complex<double> epsilon=-2.2*pow(10,5)+1.0i*6.55*pow(10,18)/w_0/omega;
//  return ((1.0 - 1.0i)* std::sqrt(omega * mu /(8 * M_PI * sigma)));
//  return std::sqrt(mu /(1.0+1.0i*4.0 * M_PI * sigma/(1.0-1.0i*omega*tau)));
// return sqrt(mu/epsilon);
   return 0.0;
}

double f (double psi)const{
    if(psi >= 0 && psi <= 1){
        return 1-pow(psi,k);
    }
    else{return 0;
    }
}
double p (double psi)const{
    return p_0*f(psi)/(2.0*R*R);
}

double B_v (double z)const{
    if (Bv_data.is_ready()) {
            return Bv_data.val(z) / R; 
        }else{throw std::runtime_error("Error: Bv_data is not ready. Check file: " );}
    // return (1.0+(M-1.0)*pow(sin(M_PI*z/2.0),q))/R;
    // return Bv_int(z);
}
double B_v2 (double z)const{
    double bv=B_v(z);
    return 1.0/(bv*bv);
}
double d_B_v (double z)const{
    if (dBv_data.is_ready()) {
            return dBv_data.val(z) / R; 
        }else{throw std::runtime_error("Error: dBv_data is not ready. Check file: " );}
}
double b_v (double z)const{
    
    return B_v(z);

}
double B (double z,double psi )const{
    double bb_v=b_v(z);

    if(bb_v>1.0){return bb_v;}
    else{double pp=p(psi)*2.0;
        return sqrt((bb_v*bb_v-pp)/(1.0-pp));}
}
double a_v (double z)const{
    return sqrt(2.0/B_v(z));
} 
double d_a_v (double z)const{
    return -sqrt(2.0/B_v(z))/(2.0*B_v(z))*d_B_v(z);
} 
mutable double last_z_ddav = -1.0, last_ddav_val = 0.0;
double dd_a_v(double z) const {
        if (std::abs(z - last_z_ddav) < 1e-12) return last_ddav_val;
        auto func = [this](double x) { return this->d_a_v(x); };
        last_ddav_val = boost::math::differentiation::finite_difference_derivative(func, z);
        last_z_ddav = z;
        return last_ddav_val;
    }
mutable double last_z_an = -1.0;
mutable double last_an_val = 0.0;
double a_nenorm(double z) const {
    if (std::abs(z - last_z_an) < 1e-12) return last_an_val;

    boost::math::quadrature::gauss_kronrod<double, 31> integrator;
    auto f_fixed_z = [this, z](double psi) { return 1.0 / B(z, psi); };
    last_an_val = sqrt(2.0 * integrator.integrate(f_fixed_z, 0, 1));
    last_z_an = z;
    return last_an_val;
}
double a (double z)const{
    return a_nenorm(z);
}
mutable double last_z_a2n = -1.0, last_a2n_val = 0.0;
double a2_nenorm(double z) const {
        if (std::abs(z - last_z_a2n) < 1e-12) return last_a2n_val;
        boost::math::quadrature::gauss_kronrod<double, 31> integrator;
        auto f_fixed_z = [this, z](double psi) { return 1.0 / B(z, psi); };
        last_a2n_val = 2.0 * integrator.integrate(f_fixed_z, 0, 1);
        last_z_a2n = z;
        return last_a2n_val;
    }
double a2 (double z)const{
    return a2_nenorm(z);
}
mutable double last_z = -1.0;
mutable double last_da = 0.0;
double d_a (double z)const{
    if (std::abs(z - last_z) < 1e-12) return last_da; // Возвращаем кэш
    auto func = [this](double x) { return this->a(x); };
    last_z = z;
    last_da = boost::math::differentiation::finite_difference_derivative(func, z);
    return last_da;
}

mutable double last_z_da2 = -1.0;
mutable double last_da2_val = 0.0;
double d_a2(double z) const {
    if (std::abs(z - last_z_da2) < 1e-12) return last_da2_val; 
    
    auto func = [this](double x) { return this->a2(x); };
    last_da2_val = boost::math::differentiation::finite_difference_derivative(func, z);
    last_z_da2 = z;
    return last_da2_val;
}
double p_sred (double z,double psi)const{
    return p(psi)*(1.0-B(z,psi)/B_s);
}
mutable double last_z_pavg = -1.0, last_pavg_val = 0.0;
double p_avg(double z) const {
        if (std::abs(z - last_z_pavg) < 1e-12) return last_pavg_val;
        boost::math::quadrature::gauss_kronrod<double, 31> integrator;
        auto f_fixed_z = [this, z](double psi) { return p_sred(z, psi) / B(z, psi); };
        double res = 2.0 / (a2(z)) * integrator.integrate(f_fixed_z, 0, 1.0);
        last_pavg_val = (res < 0.0) ? 0.0 : res;
        last_z_pavg = z;
        return last_pavg_val;
    }
    double ddav_av(double z) const{
 if (ddav_av_data.is_ready()) {
            return ddav_av_data.val(z); 
        }else{throw std::runtime_error("Error: ddav_av_data is not ready. Check file: " );}
    }
    double dBv_bv(double z) const{
         if (dBv_Bv_data.is_ready()) {
            return -2.0*dBv_Bv_data.val(z); 
        }else{throw std::runtime_error("Error: dBv_Bv_data is not ready. Check file: ");}
    }
};

struct Wall_St {
    // double RR_w = 2.25;

    double r_w(double z) const {
        return 1.0; 
    }

};
template <class Model>
struct Wall_Pr {
    const Model& m; 
    
    Wall_Pr(const Model& model) : m(model) {}

    double r_w(double z) const {
        return m.a(z); 
    }

};


template <class Model, class Wall>
struct LoDestroEquation {
    const Model& m;    
    const Wall& w;    
 std::complex<double> w1; 
 double a2_0_cached;
    mutable double last_z = -1e18;
    mutable std::complex<double> cached_q2, cached_p2, cached_a2;
    mutable double cached_d1;
    void set_omega(std::complex<double> new_w) {
        if (this->w1 == new_w) return; 
        this->w1 = new_w;
        this->last_z = -1e18; 
    }
    LoDestroEquation(const Model& model, const Wall& wall, std::complex<double> omega) 
        : m(model), w(wall), w1(omega), a2_0_cached(model.a2(0)) {}

    std::complex<double> lambda_1(double z) const {
        using namespace std::complex_literals;
        
        double rw = w.r_w(z);
        double a2_val = m.a2(z);
        
        double rw_eff = m.RR_w * rw; 
        double rw_eff2 = rw_eff * rw_eff;
        
        std::complex<double> term_plus = w1 * rw + 1.0i * m.zeta_0 * m.dzeta(w1);
        std::complex<double> term_minus = w1 * rw - 1.0i * m.zeta_0 * m.dzeta(w1);

        std::complex<double> numerator = term_plus * rw_eff2*a2_0_cached + term_minus * a2_val;
        std::complex<double> denominator = term_plus * rw_eff2*a2_0_cached - term_minus * a2_val;

        return numerator / denominator;
    }
double lambda(double z){
    return (w.r_w(z)*w.r_w(z)*m.RR_w*m.RR_w*a2_0_cached+m.a2(z))/(w.r_w(z)*w.r_w(z)*m.RR_w*m.RR_w*a2_0_cached-m.a2(z));
}
double right_func(double z) const {
        double bv = m.B_v(z);
        double av = m.a_v(z);
        double pavg = m.p_avg(z);
        double dbv = m.d_B_v(z);
        double da2 = m.d_a2(z);
        double a2 = m.a2(z);
        double ddav = m.dd_a_v(z);

        double term1 = 2.0 * (pavg / (bv * bv) * ddav / av);
        double term2 = 0.5 * std::pow((dbv / bv + da2 / a2), 2) * (1.0 - pavg / (bv * bv));

        return (term1 + term2);
    }
    // Обновление всех коэффициентов уравнения для конкретного z
    void update_cache(double z) const {
        if (std::abs(z - last_z) < 1e-13) return;

        double bv = m.B_v(z);
        double inv_bv2 = 1.0 / (bv * bv);
        double p_avg_val = m.p_avg(z);

        cached_q2 = lambda_1(z) + 1.0 - 2.0 * p_avg_val * inv_bv2;

        auto f_p2_re = [this](double z1) {
            double b = m.B_v(z1);
            return std::real(lambda_1(z1) - 2.0 * m.p_avg(z1) / (b * b));
        };
        auto f_p2_im = [this](double z1) {
            return std::imag(lambda_1(z1));
        };

        cached_p2 = std::complex<double>(boost::math::differentiation::finite_difference_derivative(f_p2_re, z), boost::math::differentiation::finite_difference_derivative(f_p2_im, z));

 
        auto f_d1 = [this](double z1) {
            double b = m.B_v(z1);
            double a2 = m.a2(z1);
            return -((m.d_B_v(z1) / b + m.d_a2(z1) / a2) * (1.0 - m.p_avg(z1) / (b * b)));
        };
        cached_d1 = boost::math::differentiation::finite_difference_derivative(f_d1, z);

        // 5. Свободный член a2 (правая часть)
        cached_a2 = w1 * w1 * m.pl * inv_bv2 - right_func(z);

        last_z = z;
    }

     
    void operator()(const state_type1 &y, state_type1 &dydx, double x) const {
        update_cache(x);
        
        dydx[0] = y[1]; // y' = v
        // y'' = -(p2*y' + d1*y + a2) / q2
        dydx[1] = -(cached_p2 * y[1] + cached_d1 * y[0] + cached_a2) / cached_q2;
    }
};

template <class Model, class Wall>
struct LoDestroEquation_int2 {
    const Model& m;    
    const Wall& w;    
 std::complex<double> w1; 
 double a2_0_cached;
    mutable double last_z = -1e18;
    mutable std::complex<double> cached_q2, cached_p2, cached_a2;
    mutable double cached_d1;
    void set_omega(std::complex<double> new_w) {
        if (this->w1 == new_w) return; 
        this->w1 = new_w;
        this->last_z = -1e18; 
    }
    LoDestroEquation_int2(const Model& model, const Wall& wall, std::complex<double> omega) 
        : m(model), w(wall), w1(omega), a2_0_cached(model.a2(0)) {}

    std::complex<double> lambda_1(double z) const {
        using namespace std::complex_literals;
        
        double rw = w.r_w(z);
        double a2_val = m.a2(z);
        
        double rw_eff = m.RR_w * rw; 
        double rw_eff2 = rw_eff * rw_eff;
        
        std::complex<double> term_plus = w1 * rw + 1.0i * m.zeta_0 * m.dzeta(w1);
        std::complex<double> term_minus = w1 * rw - 1.0i * m.zeta_0 * m.dzeta(w1);

        std::complex<double> numerator = term_plus * rw_eff2*a2_0_cached + term_minus * a2_val;
        std::complex<double> denominator = term_plus * rw_eff2*a2_0_cached - term_minus * a2_val;

        return numerator / denominator;
    }
double lambda(double z){
    return (w.r_w(z)*w.r_w(z)*m.RR_w*m.RR_w*a2_0_cached+m.a2(z))/(w.r_w(z)*w.r_w(z)*m.RR_w*m.RR_w*a2_0_cached-m.a2(z));
}
double right_func(double z) const {
        double bv2 = m.B_v2(z);
        double ddav_av = m.ddav_av(z);
        double pavg = m.p_avg(z);
        double da2 = m.d_a2(z);
        double a2 = m.a2(z);
        double dbv_bv=m.dBv_bv(z);

        double term1 = 2.0 * (pavg / (bv2) * ddav_av);
        double term2 = 0.5 * std::pow((dbv_bv+ da2 / a2), 2) * (1.0 - pavg / (bv2));

        return (term1 + term2);
    }
    // Обновление всех коэффициентов уравнения для конкретного z
    void update_cache(double z) const {
        if (std::abs(z - last_z) < 1e-13) return;

        double bv2 = m.B_v2(z);
        // double inv_bv2 = 1.0 / (bv * bv);
        double p_avg_val = m.p_avg(z);

        cached_q2 = lambda_1(z) + 1.0 - 2.0 * p_avg_val * bv2;

        auto f_p2_re = [this](double z1) {
            double bv2 = m.B_v2(z1);
            return std::real(lambda_1(z1) - 2.0 * m.p_avg(z1)*bv2);
        };
        auto f_p2_im = [this](double z1) {
            return std::imag(lambda_1(z1));
        };

        cached_p2 = std::complex<double>(boost::math::differentiation::finite_difference_derivative(f_p2_re, z), boost::math::differentiation::finite_difference_derivative(f_p2_im, z));

 
        auto f_d1 = [this](double z1) {
            double bv2 = m.B_v2(z1);
            double dbv_bv = m.dBv_bv(z1);
            double a2 = m.a2(z1);
            return -((dbv_bv + m.d_a2(z1) / a2) * (1.0 - m.p_avg(z1)*bv2));
        };
        cached_d1 = boost::math::differentiation::finite_difference_derivative(f_d1, z);

        // 5. Свободный член a2 (правая часть)
        cached_a2 = w1 * w1 * m.pl * bv2 - right_func(z);

        last_z = z;
    }

     
    void operator()(const state_type1 &y, state_type1 &dydx, double x) const {
        update_cache(x);
        
        dydx[0] = y[1]; // y' = v
        // y'' = -(p2*y' + d1*y + a2) / q2
        dydx[1] = -(cached_p2 * y[1] + cached_d1 * y[0] + cached_a2) / cached_q2;
    }
};


bool compareSecond(const std::pair<double, double>& a, const std::pair<double, double>& b) {
    return std::abs(a.second) < std::abs(b.second);
}
// double r=0.11;
// double phase=M_PI/2.1;
template <template<class, class> class EquationType, class Model, class Wall>
void reshatel(Model& model, Wall& wall, int resuis, double dx,double r_init, double phase_init) {
    using namespace std::complex_literals;
    double r=r_init;
    double phase=phase_init;
    std::complex<double> w1=r*exp(1.0i*phase);
    model.update();
    EquationType<Model, Wall> equation(model, wall, w1);
    std::cout <<"RR_w: " <<model.RR_w<< std::endl;
    std::cout <<"k: "<< model.k<<std::endl;
    std::cout <<"M: "<< model.M<<std::endl;
    std::cout <<"q: "<< model.q<<std::endl;
    std::cout <<"a(0): "<< model.a(0)<<std::endl;
    std::cout <<"r_w(0): "<<wall.r_w(0)<<std::endl;
    std::cout <<"w_0: "<< model.w_0<<std::endl;
    std::cout <<"zeta_0: "<< model.zeta_0<<std::endl;
    std::cout <<"Lambda : "<< equation.lambda(0)<<std::endl;
    std::cout <<"Lambda_1 : "<< equation.lambda_1(0)<<std::endl;
    // double Z_s=z_s();
    //     double phi_z_s;
double delta=0.00000000; // задаёт точность зануления на правой границе 
double phase_step_min=pow(0.1,4);
double r_step_min=pow(0.1,4);// сделаны чтобы обрывать бесконечные уменьшения шага
double phase_step=-M_PI/30;
double r_step=0.1;
double previous_max=1;
int counter=0;
int counter1=0;
int swich=0; // переключатель для пристрелки по фазе и радиусу
int swich1=0;  // переключатель для повторной пристрелки 

std::vector<std::pair<double,std::pair<double,double>>> bc_arr;
    std::vector<std::pair<double,std::complex<double>>> res,res1,magf,magf_v;
    std::vector<std::pair<double,double>> real_res,im_res,dres,dres1;
    double beta=model.p_0*(model.R*model.R-1)/(model.R*model.R-model.p_0);
    
    std::cout<<"Beta: " <<beta<<std::endl;
    std::complex<double> bc=1.0+1.0i;
    std::complex<double> previous_bc=bc*10.0;

    while((std::abs(bc))/abs(previous_max)>delta ){//начало цикла для решения уравнения

std::cout<<std::endl<<"dzeta: " <<model.dzeta(w1);
        res.clear();dres.clear();dres1.clear();real_res.clear();im_res.clear();
  
        auto my_observer=[&]( const state_type1 &y, const double x ) {// запись решения в данной точке
    res.push_back(std::make_pair(x,y[0]));
    dres.push_back(std::make_pair(x,std::real(y[1])));
    dres1.push_back(std::make_pair(x,std::imag(y[1])));

    real_res.push_back(std::make_pair(x,std::real(y[0])));
    im_res.push_back(std::make_pair(x,std::imag(y[0])));
    //  if(Z_s<res.back().first&&Z_s<res[res.size()-2].first){
    //     phi_z_s=std::real(y[0]);
    // }
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
    // integrate_const( stepper, equation, y0, x0, 0.5, 0.005, my_observer );// само решение
    integrate_const( stepper, equation, y0, x0, x1, dx, my_observer );// само решение
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
    equation.set_omega(w1);
    std::cout << std::endl <<"R_step: "<< r_step <<" Phase_step: "<< phase_step << " Prev: "<< std::abs(previous_bc);
    }
    std::cout<<"w^2: " <<w1*w1<< " phi'="<<dres.back().second <<std::endl;

std::cout<<std::endl<<"beta: " <<model.bbbeta<<std::endl;
std::cout<<std::endl<<"dzeta: " <<model.dzeta(w1)<<std::endl;
// std::cout<<"z_s: " <<Z_s<< std::endl;
//  std::cout<<"phi(z_s): " <<phi_z_s<< std::endl;
//  std::cout<<"[d_phi]|z_s: " <<d_phi_z_s(phi_z_s,0.05)<< std::endl;
std::cout <<"|d_phi|1|:"<< sqrt(pow(dres.back().second,2)+pow(dres1.back().second,2));
//     Gnuplot gp1;
//  gp << "set terminal wxt 1  title 'dphi'\n";
// gp << "plot '-' with lines title 'Re(dphi)', '-' with lines title 'Im(dphi)'\n";  gp.send1d(dres);gp.send1d(dres1); 
// gp1 << "set terminal wxt 2  title 'phi'\n";  ;  //gp1 << "set ytics 1\n"; 
// gp1 << "plot '-' with lines title 'Re(phi)', '-' with lines title 'Im(phi)'\n";gp1.send1d(real_res);gp1.send1d(im_res);  
// std::cout << std::endl<<"real:"<< std::real( (1.0* std::sqrt(w1 * mu / (8 * M_PI * sigma)) - 1.0i* std::sqrt(w1 * mu / (8 * M_PI * sigma))))<<std::endl<<"Imag:"<<std::imag( (1.0* std::sqrt(w1 * mu / (8 * M_PI * sigma)) - 1.0i* std::sqrt(w1 * mu / (8 * M_PI * sigma)))) ;
std::cout << std::endl<< "w^2="<<w1*w1 << "  w:"<<w1 << "r="<<r<<",phase=Pi*"<<phase/M_PI;
// std::cout << std::endl<<"sIgma"<< sigma;

std::vector<std::pair<double,double>> aaa,xi,rw,lambuda;
std::vector<std::pair<double,std::complex<double>>> lambuda1;

double size=1/static_cast<double>(res.size());

double xz=0;
for(int j=0;j<res.size();j++){
   //std::cout << "i:" << j <<"z="<<xz<<std::endl;
    double az=model.a(xz);
    double bv=model.B_v(xz);
    double rww=wall.r_w(xz);
    double lambdaa=equation.lambda(xz);
    std::complex<double> lambdaa1=equation.lambda_1(xz);
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
}else if(resuis==7){
     resis="_GDL_ideal";
}else if(resuis==8){
     resis="_GDL_res";
}
std::string folder1 = "M=" + std::to_string(model.M) + "_" +
                      "q=" + std::to_string(model.q) + "_" +
                      "k=" + std::to_string(model.k) + "_" +
                      "R=" + std::to_string(model.R);

// Исправленная склейка строк (не хватало плюсов и to_string):
std::string folder2 = "R_w=" + std::to_string(model.RR_w) + resis;
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
        file << "omega:" <<w1 << " dzeta:"<< model.dzeta(w1) << " bc:" << sqrt(pow(dres.back().second,2)+pow(dres1.back().second,2));
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

}


// #pragma omp threadprivate(r, phase)
int main() {
    SharedInterpolationData shared_data;
    std::string papka_int = R"(C:\Users\MJ\Desktop\ballon\pole\clean\)";
    shared_data.load_all(papka_int);
double start_beta = 0.2;
    double end_beta = 0.98;
    double step = 0.05;
    int steps = static_cast<int>((end_beta - start_beta) / step);
double r1 = 0.8;
        double phase1 = M_PI / 2.0;
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i <= steps; ++i) {
        double current_beta = start_beta + i * step;

        // Передаем ссылку на shared_data в модель
        PlasmaModel_A1_int model(shared_data);
        model.RR_w = 2.25;
        model.bbbeta = current_beta;
        model.R = 3.2;
        model.B_v_ = 3671.324;
        model.update();

        Wall_St wall;

        #pragma omp critical(print)
        std::cout << "\nStarting Thread " << omp_get_thread_num() << " for beta: " << current_beta << std::endl;

        // Передаем начальные r и phase как аргументы!
        reshatel<LoDestroEquation_int2>(model, wall, 7, 0.01, 0.8, M_PI / 2.0);
    }
    return 0;

}