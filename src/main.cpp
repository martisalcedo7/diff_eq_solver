#include <iostream>
#include <vector>
#include <matplotlib-cpp-master/matplotlibcpp.h>

#define STATES_N 4

namespace plt = matplotlibcpp;

void f_dot(double *dot, double t, ){

}

double x_dot(double t, double x, double y, double z, double w){
    // Equation
    // v = i*R + 1/C *1/s * i + L * s * i
    return y;
}

double y_dot(double t, double x, double y, double z, double w){
    // Equation
    // v = i*R + 1/C *1/s * i + L * s * i
    // static const double v = 15.0; //volts
    // static const double r = 10.0; //ohms
    // static const double l = 0.15; //henry
    // static const double c = 50.0e-6; //farads
    static const double g = 9.81;
    static const double m1 = 5;
    static const double m2 = 5;
    static const double k1 = 20000;
    static const double k2 = 10000;
    static const double d = 50;
    static const double f = 0;
    return -g+(m2/m1)*g-(k1/m1)*x+(k2/m1)*z-(d/m1)*y+f/m1;
}

double z_dot(double t, double x, double y, double z, double w){
    // Equation
    // v = i*R + 1/C *1/s * i + L * s * i
    return w;
}

double w_dot(double t, double ){
    // Equation
    // v = i*R + 1/C *1/s * i + L * s * i
    // static const double v = 15.0; //volts
    // static const double r = 10.0; //ohms
    // static const double l = 0.15; //henry
    // static const double c = 50.0e-6; //farads
    static const double g = 9.81;
    static const double m1 = 5;
    static const double m2 = 5;
    static const double k1 = 20000;
    static const double k2 = 10000;
    static const double d = 50;
    static const double f = 0;
    return -(m2/m1)*g+(k1/m1)*x-(k2/m1)*z+(d/m1)*y+f/m1-(k2/m2)*z;
}

void solver(double *states, const double T, const unsigned long int steps, const int states_number){
    
    double time = 0.0;
    // double x_k1, x_k2, x_k3, x_k4 = 0.0;
    // double y_k1, y_k2, y_k3, y_k4 = 0.0;
    // double z_k1, z_k2, z_k3, z_k4 = 0.0;
    // double w_k1, w_k2, w_k3, w_k4 = 0.0;

    double rk[4][states_number] = {};

    for(unsigned long int step=0; step < (steps - 1); ++step){
        // std::cout << step << std::endl;
        time = T * step;

        rk[0][state] = x_dot(time, *(states + step));
        rk[1][state] = x_dot(time + (T/2.0), *(states + step)+T*(x_k1/2.0), *(y + step)+T*(y_k1/2.0), *(z + step)+T*(z_k1/2.0), *(w + step)+T*(w_k1/2.0));
        rk[2][state] = x_dot(time + (T/2.0), *(x + step)+T*(x_k2/2.0), *(y + step)+T*(y_k2/2.0), *(z + step)+T*(z_k2/2.0), *(w + step)+T*(w_k2/2.0));
        rk[3][state] = w_dot(time + T, *(x + step)+T*x_k3, *(y + step)+T*y_k3, *(z + step)+T*z_k3, *(w + step)+T*w_k3);

        for(int state=0; state < states_number; ++state){

            // x_k1 = x_dot(time, *(x + step), *(y + step), *(z + step), *(w + step));
            // y_k1 = y_dot(time, *(x + step), *(y + step), *(z + step), *(w + step));
            // z_k1 = z_dot(time, *(x + step), *(y + step), *(z + step), *(w + step));
            // w_k1 = w_dot(time, *(x + step), *(y + step), *(z + step), *(w + step));

            // x_k2 = x_dot(time + (T/2.0), *(x + step)+T*(x_k1/2.0), *(y + step)+T*(y_k1/2.0), *(z + step)+T*(z_k1/2.0), *(w + step)+T*(w_k1/2.0));
            // y_k2 = y_dot(time + (T/2.0), *(x + step)+T*(x_k1/2.0), *(y + step)+T*(y_k1/2.0), *(z + step)+T*(z_k1/2.0), *(w + step)+T*(w_k1/2.0));
            // z_k2 = z_dot(time + (T/2.0), *(x + step)+T*(x_k1/2.0), *(y + step)+T*(y_k1/2.0), *(z + step)+T*(z_k1/2.0), *(w + step)+T*(w_k1/2.0));
            // w_k2 = w_dot(time + (T/2.0), *(x + step)+T*(x_k1/2.0), *(y + step)+T*(y_k1/2.0), *(z + step)+T*(z_k1/2.0), *(w + step)+T*(w_k1/2.0));

            // x_k3 = x_dot(time + (T/2.0), *(x + step)+T*(x_k2/2.0), *(y + step)+T*(y_k2/2.0), *(z + step)+T*(z_k2/2.0), *(w + step)+T*(w_k2/2.0));
            // y_k3 = y_dot(time + (T/2.0), *(x + step)+T*(x_k2/2.0), *(y + step)+T*(y_k2/2.0), *(z + step)+T*(z_k2/2.0), *(w + step)+T*(w_k2/2.0));
            // z_k3 = z_dot(time + (T/2.0), *(x + step)+T*(x_k2/2.0), *(y + step)+T*(y_k2/2.0), *(z + step)+T*(z_k2/2.0), *(w + step)+T*(w_k2/2.0));
            // w_k3 = w_dot(time + (T/2.0), *(x + step)+T*(x_k2/2.0), *(y + step)+T*(y_k2/2.0), *(z + step)+T*(z_k2/2.0), *(w + step)+T*(w_k2/2.0));

            // x_k4 = x_dot(time + T, *(x + step)+T*x_k3, *(y + step)+T*y_k3, *(z + step)+T*z_k3, *(w + step)+T*w_k3);
            // y_k4 = y_dot(time + T, *(x + step)+T*x_k3, *(y + step)+T*y_k3, *(z + step)+T*z_k3, *(w + step)+T*w_k3);
            // z_k4 = z_dot(time + T, *(x + step)+T*x_k3, *(y + step)+T*y_k3, *(z + step)+T*z_k3, *(w + step)+T*w_k3);
            // w_k4 = w_dot(time + T, *(x + step)+T*x_k3, *(y + step)+T*y_k3, *(z + step)+T*z_k3, *(w + step)+T*w_k3);

            *(*(states + step + 1) + state) = *(*(states + step) + state) + (1.0/6.0) * T * (rk[0][state] + 2.0 * rk[1][state] + 2.0 * rk[2][state] + rk[3][state]);

            // *(x + step + 1) = *(x + step) + (1.0/6.0) * T * (x_k1 + 2.0 * x_k2 + 2.0 * x_k3 + x_k4);
            // *(y + step + 1) = *(y + step) + (1.0/6.0) * T * (y_k1 + 2.0 * y_k2 + 2.0 * y_k3 + y_k4);
            // *(z + step + 1) = *(z + step) + (1.0/6.0) * T * (z_k1 + 2.0 * z_k2 + 2.0 * z_k3 + z_k4);
            // *(w + step + 1) = *(w + step) + (1.0/6.0) * T * (w_k1 + 2.0 * w_k2 + 2.0 * w_k3 + w_k4);
        }

    }
}

int main(){
    const double T = 1e-4;
    const double sim_time = 10.0;

    const unsigned long int steps = (unsigned long int)(sim_time/T);
    const int states_number = 4;

    double states[states_number][steps] = {};
    double initial_value[states_number] = {};

    for(int state = 0; state < states_number; ++state){
        states[state][0] = initial_value[state];
    }

    solver(states, T, steps, states_number);

    // for(unsigned long int step=0; step < steps; ++step){
    //     std::cout << z[0] << std::endl;
    // }

    std::vector<double> v(x, x+steps);
    std::vector<double> v1(z, z+steps);

    // plt::figure();
    plt::plot(v);
    plt::plot(v1);
    plt::grid(true);
    plt::show();

    return 0;
}

