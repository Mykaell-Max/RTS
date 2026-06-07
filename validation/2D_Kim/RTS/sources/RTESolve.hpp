#pragma once
#include <cuda_runtime.h>
#include <vector>

#define PI 3.14159265358979323846
#define PI4 12.56637061435917295385
#define SIGMA 5.6703744191843561648e-8

class FAM{
    private:
        double* I, Ib, Sm, kappa, sigma_s, temp, epsilon;
        double* Ax, Ay, Az, volom, alpha_r, dOmega, phi;
        double res, tol;
        int nxi, nyi, nzi, nt, np;
        void RHS_SM(double* Sm, double* I, double* Ib);
        double angular_loop(double* Sm, double* I);
        void FAMbound_in(double* I, double* Ib);
        void FAMbound_out(double* I);
    
    public:
        FAM(int nxi, int nyi, int nzi, int nt, int np, double tol, double alpha_r,
            double* I, double* Ib, double *kappa, double *sigma_s, double *temp, double *epsilon,
            double* dOmega, double* phi,
            double* Sm,  double* Ax, double* Ay, double* Az, double* volom);
        ~FAM();
        void solve();
        double *get_Srad();
        double *get_G();
        double *get_Qradw();
    
    private:
        void allocateMemory();
        void freeMemory();
        void copyToDevice();
        void copyFromDevice();
};