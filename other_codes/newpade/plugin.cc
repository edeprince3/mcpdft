/*
 * @BEGIN LICENSE
 *
 * newpade by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/libparallel/parallel.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libpsio/psio.hpp"
//#include "psi4/libmints/mints.h"
#include "psi4/libpsio/psio.hpp"
#include<stdio.h>
#include<stdlib.h>
#include "psi4/libqt/qt.h"
#include </edfs/users/daniel/deprince-group/plugins/rttd2rdm/blas.h>
#include "psi4/libciomr/libciomr.h"


namespace psi{ namespace newpade {

extern "C"
int read_options(std::string name, Options& options)
{
    if (name == "NEWPADE"|| options.read_globals()) {
        /*- The amount of information printed to the output file -*/
        options.add_int("PRINT", 1);
                /*- File name for data -*/
        options.add_str("FILENAME", "0.00005.CO.J");
        /*- Use windowing function, default no window -*/
        options.add_int("WINDOW", 1);
        /*- Maximum quociente, default 1e5 -*/
        options.add_double("MIN_QUOCIENT", 1e5);
        /*- Damping rate -*/
        options.add_double("DAMPING", 0.005);
        /*- Initial frequency -*/
        options.add_double("INITIAL_FREQ", 0.0);
        /*- Final frequency -*/
        options.add_double("FINAL_FREQ", 40.0);
        /*- Frequency interval -*/
        options.add_double("DELTA_FREQ", 0.005);

    }

    return true;
}

extern "C"
SharedWavefunction newpade(SharedWavefunction ref_wfn, Options& options)
{
    int print = options.get_int("PRINT");

    int window = options.get_int("WINDOW");
    double min_quocient = options.get_double("MIN_QUOCIENT");
//using namespace fnocc;


    FILE * fp = fopen(options.get_str("FILENAME").c_str(),"r");
    if ( fp == NULL ) {
        throw PsiException("file does not exist.",__FILE__,__LINE__);
    }

    long int N = 0;
    double dum;
//    FILE * fp = fopen(argv[1],"r");
    char * line = (char*)malloc(1000*sizeof(char));
    do {
        fscanf(fp,"%s %le  %le  %le \n",line,&dum,&dum,&dum);
        N++;
    }while(!feof(fp));
    rewind(fp);

    double * td_time      = (double*)malloc(N*sizeof(double));
    double * corr_func_r0 = (double*)malloc(N*sizeof(double));
    double * corr_func_r1 = (double*)malloc(N*sizeof(double));

    int n = 0;
//    double t = 0.0;
//    double rate =  0.0; //0.005;

    // Choose windowing function:
    // 1) Lorentzian;
    // 2) Hanning;

    //int window = 1;
    //Set windowing parameters;

    double t0 = 0.0; //1000.0;   //My value: Lorentz (0.0), Hanning (total_time)
    double S = options.get_double("DAMPING");
    //double S  = 0.0001; //2000.0;  //My value: Lorentz (5/total_time), Hanning (2*total_time)
    do {
        double t,valrr,valri,vallr,valli;
        fscanf(fp,"%s %le  %le  %le \n",line,&t,&valrr,&valri);
        td_time[n] = t;

        if (window == 1){
            corr_func_r0[n] = valrr*exp(-S*t);
            corr_func_r1[n] = valri*exp(-S*t);
        }
        else if (window == 2){
            corr_func_r0[n] = valrr*sin((3.1415/S)*(t-t0))*sin((3.1415/S)*(t-t0));
            corr_func_r1[n] = valri*sin((3.1415/S)*(t-t0))*sin((3.1415/S)*(t-t0));
        }
        else {
            corr_func_r0[n] = valrr;
            corr_func_r1[n] = valri;
        }

        n++;
    }while(!feof(fp));

//  Begin Pade algorithm
    long int K = (N-1)/2;

    // For a real signal
    double * C = (double*)malloc(K*K*sizeof(double));
    double * b_re = (double*)malloc((K+1)*sizeof(double));
    double * temp = (double*)malloc(K*sizeof(double));
    for (int p = 0; p < K; p++){
        temp[p] = -corr_func_r0[p+K+1];
        for (int q = 0; q < K; q++){
            C[q*K+p] = corr_func_r0[K+p-q];
        }
    }

    // solving for B
    double * a_re = (double*)malloc((K+1)*sizeof(double));
    int * ipiv = (int*)malloc(K*sizeof(int));
    for (int i = 0; i < K+1; i++){
        a_re[i] = 0.0;
        b_re[i] = 0.0;
    }
    for (int i = 0; i < K; i++){
        ipiv[i] = 0.0;
    }

    int info = 0;
    info = C_DGESV(K,1,C,K,ipiv,temp,K);

    if (info == 0) {printf("Matrix inversion complete!\n");}
    else {printf("Matrix inversion failed! INFO = %i \n",info);}

//    // Normalization
    b_re[0] = 1.0;
    a_re[0] = corr_func_r0[0];

    for (int i = 1; i < K+1; i++){
        b_re[i]  = temp[i-1];
    }

    for (int p = 1; p < K+1; p++){
        for (int q = 0; q < p+1; q++){
            a_re[p] += corr_func_r0[(p-q)]*b_re[q];
        }
    }
    double start  = options.get_double("INITIAL_FREQ");
    double end    = options.get_double("FINAL_FREQ");
    double change = options.get_double("DELTA_FREQ");
    double ts = td_time[1]-td_time[0];

    for (double k = start/27.21138; k < end/27.21138; k += change/27.21138) {

        double A_re = 0.0;
        double A_im = 0.0;
        double B_re = 0.0;
        double B_im = 0.0;
        double x_re = 0.0;
        double x_im = 0.0;
        double val_re  = 0.0;
        double val_im  = 0.0;
        double val  = 0.0;

        // Check whether or not k is a root:
        for (int j = 0; j < K+1; j++) {
            x_re =  cos(k*j*ts);
            x_im =  -sin(k*j*ts);

            A_re += a_re[j]*x_re;
            A_im += a_re[j]*x_im;

            B_re += b_re[j]*x_re;
            B_im += b_re[j]*x_im;

        }
        val_re  = ((A_re*B_re + A_im*B_im)/(B_re*B_re+B_im*B_im))/(S*M_PI*(K+1));
        val_im  = ((A_im*B_re - A_re*B_im)/(B_re*B_re+B_im*B_im))/(S*M_PI*(K+1));
        val = (2./3.)*k*val_re;
        printf("%20.12lf %20.12lf %20.12lf \n",k*27.21138,val_re,val_im);
    }

    /* Your code goes here */

    // Typically you would build a new wavefunction and populate it with data
    return ref_wfn;
}

}} // End namespaces

