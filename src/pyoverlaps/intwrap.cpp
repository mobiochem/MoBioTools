#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>
#include"constants.h"

#define DELLEXPORT extern "C"

using namespace std;


double getNormConst(long* n, double e){
    double N = 0.;
    double f1, f2, f3_1, f3_2, norm;
    for(long i=0; i<3; i++){
        N += n[i];
    }

    f1 = pow(2.*e/pi, 0.75);
    f2 = pow(4.*e, N/2.);
    f3_1 = dfact[2*n[0]] * dfact[2*n[1]] * dfact[2*n[2]];
    f3_2 = pow(f3_1, -0.5);
    norm = f1*f2*f3_2;

    return norm;
}
// Ovlp integrals involving two primitives.
// Returns vectors gx, gy and gz in a single arr g.
void calcOvlp(double* R1, double alpha1, double* R2, double alpha2,
        double* g, long size_g, long nmax, long lj, long dj){
    double gamma = alpha1+alpha2;
    double *gx = g;
    double *gy = g + size_g;
    double *gz = g + 2*size_g;
    double r1r2[3], r1r12[3];
    long ptr;

    // Apply rhs recursion    
    r1r2[0]  = R2[0] - R1[0];
    r1r2[1]  = R2[1] - R1[1];
    r1r2[2]  = R2[2] - R1[2];
    r1r12[0] = (alpha1*R1[0] + alpha2*R2[0])/gamma - R2[0];
    r1r12[1] = (alpha1*R1[1] + alpha2*R2[1])/gamma - R2[1];
    r1r12[2] = (alpha1*R1[2] + alpha2*R2[2])/gamma - R2[2];

    gx[0] = 1;
    gy[0] = 1;
    gz[0] = 1;

    // Define S(0,1)
    if(nmax > 0){
        gx[1] = r1r12[0] * gx[0];
        gy[1] = r1r12[1] * gy[0];
        gz[1] = r1r12[2] * gz[0];
    
    }

    // Index recursion
    for(long i = 1; i < nmax; i++){
        gx[i + 1] = r1r12[0] * gx[i] + i/(2*gamma) * gx[i - 1];
        gy[i + 1] = r1r12[1] * gy[i] + i/(2*gamma) * gy[i - 1];
        gz[i + 1] = r1r12[2] * gz[i] + i/(2*gamma) * gz[i - 1];
    }

    // transfer recursion
    for(long j = 1; j<=lj; j++){
        ptr = dj * j;
        for(long i = ptr; i<= ptr + nmax -j; i++){
            gx[i] = gx[i+1-dj] + r1r2[0] * gx[i-dj];
            gy[i] = gy[i+1-dj] + r1r2[1] * gy[i-dj];
            gz[i] = gz[i+1-dj] + r1r2[2] * gz[i-dj];
        
        }
    
    }


}


DELLEXPORT void intOvlpCrt(long* bas, long* atm, double* env,
                  double* S, long nbas, long nshl){
    // Compute overlap matrix between two sets of shells
    long size1 = nbas, size2 = nbas;
    long n1[3], n2[3]; // Angular momentum vectors
    double R1[3], R2[3]; // Coordinate vectors
//    double C1, C2, alpha1, alpha2, N1, N2, r2, gamma, pi_gamma = 0.0; 
    double C1, C2, alpha1, alpha2, r2, gamma, pi_gamma = 0.0; 
    long at1, at2, l1, l2, np1, np2, lj, dj, nmax, size_g, 
         size_g_red, nbas1, nbas2;
    long idx, idy, idz;
    long idc1, idc2; //Index contractions 1 and 2
    double prefac = 0.0;

    // Iterate over first shell
    idc1 = 0;
    for(long ish = 0; ish<nshl; ish++){
        at1   = bas[ish*8]; // atom1 index
        R1[0] = env[atm[at1*6 + 1]];
        R1[1] = env[atm[at1*6 + 1] + 1];
        R1[2] = env[atm[at1*6 + 1] + 2];
        l1    = bas[ish*8 + 1]; // ang_mom 1
        np1   = bas[ish*8 + 2]; // n. primitives
        nbas1 = OFFCRT[l1 + 1] - OFFCRT[l1]; // n. contr.

        // Iterate over second shell
//        idc2 = idc1;
        idc2 = 0;
//        for(long jsh = ish; jsh<nshl; jsh++){
        for(long jsh = 0; jsh<nshl; jsh++){
            at2   = bas[jsh*8]; // atom1 index
            R2[0] = env[atm[at2*6 + 1]];
            R2[1] = env[atm[at2*6 + 1] + 1];
            R2[2] = env[atm[at2*6 + 1] + 2];
            l2    = bas[jsh*8 + 1]; // ang_mom 2
            np2   = bas[jsh*8 + 2]; // n. primitives
            nbas2 = OFFCRT[l2 + 1] - OFFCRT[l2]; // n. contr.
            
            // squared distance
            r2    = 0.0; 
            for(long i = 0; i<3; i++){
                r2 += pow(R2[i] - R1[i],2);
            }

            // Variables dependent on l1, l2
            nmax = l1 + l2;
            dj   = nmax + 1;
            if(l1>=l2){
                lj = l1;
            }
            else{
                lj = l2;
            }
//            lj   = l2;
            size_g_red = dj * (lj + 1); // length of g_i, i=x,y,z
            size_g = 3*size_g_red;      // length of g
            double* g = new double[size_g];
            for(long i = 0; i<size_g; i++){
                g[i] = 0.0;
            }


            // Iterate over primitives of ish
            for(long i = 0; i<np1; i++){
                alpha1 = env[bas[ish*8 + 5] + i];
                C1     = env[bas[ish*8 + 6] + i];
                // Iterate over primitives of jsh
                for(long j = 0; j<np2; j++){
                    alpha2 = env[bas[jsh*8 + 5] + j];
                    C2     = env[bas[jsh*8 + 6] + j];

                    // Compute prefactor
                    gamma  = alpha1 + alpha2;
                    pi_gamma = pow(pi/gamma,1.5);
                    prefac = pi_gamma * pow(e, -alpha1*alpha2*(r2)/gamma);
                
                    // Compute set of integrals corresponding to the
                    // current pair of primitives
                    calcOvlp(R1, alpha1, R2, alpha2,
                       g, size_g_red, nmax, lj, dj);

                    for(long k = 0; k<nbas1; k++){
                        n1[0] = ANGCRT[OFFCRT[l1] + k][0];
                        n1[1] = ANGCRT[OFFCRT[l1] + k][1];
                        n1[2] = ANGCRT[OFFCRT[l1] + k][2];
//                        N1    = getNormConst(n1, alpha1);
                        for(long m = 0; m<nbas2; m++){
                            n2[0] = ANGCRT[OFFCRT[l2] + m][0];
                            n2[1] = ANGCRT[OFFCRT[l2] + m][1];
                            n2[2] = ANGCRT[OFFCRT[l2] + m][2];
//                            N2    = getNormConst(n2, alpha2);
                            idx = n1[0]*dj + n2[0];
                            idy = size_g_red + n1[1]*dj + n2[1];
                            idz = 2*size_g_red + n1[2]*dj + n2[2];

                            S[(idc1 + k)*nbas + (idc2 + m)] += prefac * C1 * C2 * g[idx]*g[idy]*g[idz];
                        }
                    }
                }
            }
            delete[] g;
            idc2 += nbas2; // Update index of the first contr. of 2
        }
        idc1 += nbas1;    // Update index of the first contr. of 1;
    }
}

// Overlap between two sets of shells, with different basis sets, but the
// same basis family. The same number of basis functions is assumed
DELLEXPORT void intOvlpCrtMix(long* bas1, long* atm1, double* env1,
                              long* bas2, long* atm2, double* env2,
                              double* S, long nbasis1, long nshl1,
                              long nbasis2, long nshl2){
    // Compute overlap matrix between two sets of shells
    long size1 = nbasis1, size2 = nbasis2;
    long n1[3], n2[3]; // Angular momentum vectors
    double R1[3], R2[3]; // Coordinate vectors
    double C1, C2, alpha1, alpha2, r2, gamma, pi_gamma = 0.0; 
    long at1, at2, l1, l2, np1, np2, lj, dj, nmax, size_g, 
         size_g_red, nbas1, nbas2;
    long idx, idy, idz;
    long idc1, idc2; //Index contractions 1 and 2
    double prefac = 0.0;

    // Iterate over first shell
    idc1 = 0;
    for(long ish = 0; ish<nshl1; ish++){
        at1   = bas1[ish*8]; // atom1 index
        R1[0] = env1[atm1[at1*6 + 1]];
        R1[1] = env1[atm1[at1*6 + 1] + 1];
        R1[2] = env1[atm1[at1*6 + 1] + 2];
        l1    = bas1[ish*8 + 1]; // ang_mom 1
        np1   = bas1[ish*8 + 2]; // n. primitives
        nbas1 = OFFCRT[l1 + 1] - OFFCRT[l1]; // n. contr.

        // Iterate over second shell
//        idc2 = idc1;
        idc2 = 0;
//        for(long jsh = ish; jsh<nshl; jsh++){
        for(long jsh = 0; jsh<nshl2; jsh++){
            at2   = bas2[jsh*8]; // atom1 index
            R2[0] = env2[atm2[at2*6 + 1]];
            R2[1] = env2[atm2[at2*6 + 1] + 1];
            R2[2] = env2[atm2[at2*6 + 1] + 2];
            l2    = bas2[jsh*8 + 1]; // ang_mom 2
            np2   = bas2[jsh*8 + 2]; // n. primitives
            nbas2 = OFFCRT[l2 + 1] - OFFCRT[l2]; // n. contr.
            
            // squared distance
            r2    = 0.0; 
            for(long i = 0; i<3; i++){
                r2 += pow(R2[i] - R1[i],2);
            }

            // Variables dependent on l1, l2
            nmax = l1 + l2;
            dj   = nmax + 1;
            if(l1>=l2){
                lj = l1;
            }
            else{
                lj = l2;
            }
//            lj   = l2;
            size_g_red = dj * (lj + 1); // length of g_i, i=x,y,z
            size_g = 3*size_g_red;      // length of g
            double* g = new double[size_g];
            for(long i = 0; i<size_g; i++){
                g[i] = 0.0;
            }


            // Iterate over primitives of ish
            for(long i = 0; i<np1; i++){
                alpha1 = env1[bas1[ish*8 + 5] + i];
                C1     = env1[bas1[ish*8 + 6] + i];
                // Iterate over primitives of jsh
                for(long j = 0; j<np2; j++){
                    alpha2 = env2[bas2[jsh*8 + 5] + j];
                    C2     = env2[bas2[jsh*8 + 6] + j];

                    // Compute prefactor
                    gamma  = alpha1 + alpha2;
                    pi_gamma = pow(pi/gamma,1.5);
                    prefac = pi_gamma * pow(e, -alpha1*alpha2*(r2)/gamma);
                
                    // Compute set of integrals corresponding to the
                    // current pair of primitives
                    calcOvlp(R1, alpha1, R2, alpha2,
                       g, size_g_red, nmax, lj, dj);

                    for(long k = 0; k<nbas1; k++){
                        n1[0] = ANGCRT[OFFCRT[l1] + k][0];
                        n1[1] = ANGCRT[OFFCRT[l1] + k][1];
                        n1[2] = ANGCRT[OFFCRT[l1] + k][2];
//                        N1    = getNormConst(n1, alpha1);
                        for(long m = 0; m<nbas2; m++){
                            n2[0] = ANGCRT[OFFCRT[l2] + m][0];
                            n2[1] = ANGCRT[OFFCRT[l2] + m][1];
                            n2[2] = ANGCRT[OFFCRT[l2] + m][2];
//                            N2    = getNormConst(n2, alpha2);
                            idx = n1[0]*dj + n2[0];
                            idy = size_g_red + n1[1]*dj + n2[1];
                            idz = 2*size_g_red + n1[2]*dj + n2[2];

                            S[(idc1 + k)*nbasis2 + (idc2 + m)] += prefac * C1 * C2 * g[idx]*g[idy]*g[idz];
//                            cout<<"S["<<idc1<<"]["<<idc2<<"] = "<< S[(idc1 + k)*nbasis2 + (idc2 + m)]<<endl;
                        }
                    }
                }
            }
            delete[] g;
            idc2 += nbas2; // Update index of the first contr. of 2
        }
        idc1 += nbas1;    // Update index of the first contr. of 1;
    }
}
