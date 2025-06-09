#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>

double f(double x, double y) { return y + x*x; }

void euler(double* res, double interval, int iter) {
    double slope;
    *(res++) = 0;
    for(int i=1; i<iter; i++, res++) {
        slope = f(interval * (i-1), *(res-1));
        *res = *(res-1) + slope*interval;
    }
}

void imeuler(double* res, double interval, int iter) {
    double slope;
    *(res++) = 0;
    for(int i=1; i<iter; i++, res++) {
        slope = f(interval * (i-1), *(res-1));
        *res = *(res-1) + slope*interval;
        slope = (slope + f(interval * i, *res)) / 2;
        *res = *(res-1) + slope*interval;
    }
}

void RungeKutta(double* res, double interval, int iter) {
    double k1, k2, k3, k4;
    double x, y;
    *(res++) = 0;
    for(int i=1; i<iter; i++, res++) {
        x = interval * (i-1); y = *(res-1);
        k1 = interval * f(x, y);
        k2 = interval * f(x + interval/2, y + k1/2);
        k3 = interval * f(x + interval/2, y + k2/2);
        k4 = interval * f(x + interval, y + k3);
        *res = y + (k1 + 2*k2 + 2*k3 + k4)/6;
    }
}

void precise(double* res, double interval, int iter) {
    double x;
    *(res++) = 0;
    for(int i=1; i<iter; i++, res++) {
        x = i * interval;
        *res = -(x*x) - 2*x - 2 + (2*pow(M_E, x));
    }
}

void printRes(double* arr, int num) {
    for(int i=0; i<num; i++, arr++) { printf("%10f", *arr); }
    printf("\n"); return;
}

int main(void) {
    double res[4][11] = { 0, };
    
    euler(res[0], 0.1f, 11);
    imeuler(res[1], 0.1f, 11);
    RungeKutta(res[2], 0.1f, 11);
    precise(res[3], 0.1f, 11);
    
    printf("Euler\t");
    printRes(res[0], 11);
    printf("ImEuler\t");
    printRes(res[1], 11);
    printf("R-K\t");
    printRes(res[2], 11);
    printf("Real\t");
    printRes(res[3], 11);
}