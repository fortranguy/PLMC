#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include "nfft3.h"
#include "data_reci.h"

static nfft_plan structure[DIM];
static nfft_plan potential[DIM];
static double Epot_reci_tab[2*Nx][2*Ny][2*Nz];

void Epot_reci_nfft_init(const int Ncol);
void Epot_reci_nfft_finalize();
void Epot_reci_init(const double Lsize[DIM], const double alpha);
double Epot_reci_move(const int lCol, const double deltaX[DIM], const double Vol);
void Epot_reci_updateX(const int lCol, const double xNew[DIM]);
double Epot_reci(double X[][DIM], double D[][DIM], const int Ncol, const double Vol);
void snapShot(const int Ncol);

// Initialisation : phi and S

void Epot_reci_nfft_init(const int Ncol){
    
    for(int iComp=0; iComp<DIM; iComp++){
        
        nfft_init_3d(&potential[iComp], 2*Nx, 2*Ny, 2*Nz, Ncol);
        nfft_init_3d(&structure[iComp], 2*Nx, 2*Ny, 2*Nz, Ncol);
        
    }
    
    return;
    
}


// Finalisation

void Epot_reci_nfft_finalize(){
 
    for (int iComp=0; iComp<DIM; iComp++){
        nfft_finalize(&structure[iComp]);
        nfft_finalize(&potential[iComp]);
    }
    
    return;
    
}

void Epot_reci_init(const double Lsize[DIM], const double alpha){

    int nb_k;
    double kxOverL, kyOverL, kzOverL;
    double kOverLsqr, alphaSqr;
    
    nb_k = 0;
    alphaSqr = alpha*alpha;
    
    for (int kx=-Nx; kx<Nx; kx++){
    
        kxOverL = (double)kx/Lsize[0];     
           
        for (int ky=-Ny; ky<Ny; ky++){
        
            kyOverL = (double)ky/Lsize[1];  
                          
            for (int kz=-Nz; kz<Nz; kz++){
                
                kzOverL = (double)kz/Lsize[2];
                
                if (kx != 0 || ky != 0 || kz != 0){

                    kOverLsqr = kxOverL*kxOverL + kyOverL*kyOverL + kzOverL*kzOverL;
                    Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz] = exp(-PI*PI/alphaSqr*kOverLsqr);
                    Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz] /= kOverLsqr;
                    nb_k++;
                
                }
                
            }
            
        }
        
    }
    
    Epot_reci_tab[Nx][Ny][Nz] = 0.;
    
    printf(" Ewald Summation : Number of wave vectors : %d\n", nb_k);
    
    return;
    
}

double Epot_reci_move(const int lCol, const double deltaX[DIM], const double Vol){
 
    // delta M
    
    // Doing the transform : structure
    
    for (int iComp=0; iComp<DIM; iComp++){
        nfft_adjoint(&structure[iComp]);
    }
    
    // delta U
    
    int ikx, iky, ik;
    double complex k_dot_structure, k_dot_structure_excess;
    double complex k_dot_xCol;
    double complex k_dot_deltaXcol;
    double complex k_dot_coefficient;
    double complex factor;
    double Epot;
    
    Epot = 0.;
    
    for (int kx=-Nx; kx<Nx; kx++){
        
        ikx = (kx + Nx)*(2*Ny * 2*Nz);
        
        for (int ky=-Ny; ky<Ny; ky++){
            
            iky = (ky + Ny)*2*Nz;
            
            for (int kz=-Nz; kz<Nz; kz++){
                
                ik = ikx + iky + (kz + Nz);
                
                k_dot_xCol = (double)kx * potential[0].x[DIM*lCol+0] +
                             (double)ky * potential[0].x[DIM*lCol+1] +
                             (double)kz * potential[0].x[DIM*lCol+2];
                             
                k_dot_deltaXcol = (double)kx * deltaX[0] +
                                  (double)ky * deltaX[1] +
                                  (double)kz * deltaX[2];
                
                k_dot_structure  = (double)kx * structure[0].f_hat[ik] +
                                   (double)ky * structure[1].f_hat[ik] +
                                   (double)kz * structure[2].f_hat[ik];
                k_dot_structure_excess = (double)kx * structure[0].f[lCol] +
                                         (double)ky * structure[1].f[lCol] +
                                         (double)kz * structure[2].f[lCol];
                k_dot_structure_excess*= exp(I*2.*PI*k_dot_xCol);
                k_dot_structure = k_dot_structure - k_dot_structure_excess;

                factor = exp(-I*2.*PI*k_dot_xCol) * (exp(-I*2.*PI*k_dot_deltaXcol) - 1.);
                
                k_dot_coefficient = (double)kx * structure[0].f[lCol] +
                                    (double)ky * structure[1].f[lCol] +
                                    (double)kz * structure[2].f[lCol];
                k_dot_coefficient*= factor;
               
                Epot = Epot + 2.*creal(k_dot_coefficient*k_dot_structure) *
                              Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz];
            
            }            
        }        
    }
    
    return 2.*PI/Vol * Epot;
    
}

void Epot_reci_updateX(const int lCol, const double xNew[DIM]){
 
    for(int iComp=0; iComp<DIM; iComp++){
        for (int iDim=0; iDim<DIM; iDim++){        
            
            potential[iComp].x[DIM*lCol+iDim] = xNew[iDim];            
            structure[iComp].x[DIM*lCol+iDim] = xNew[iDim];            
            
        }
    }
    
    // Precompute $\psi$
    /*
    for(int iComp=0; iComp<DIM; iComp++){
        
        if(potential[iComp].nfft_flags & PRE_ONE_PSI)
            nfft_precompute_one_psi(&potential[iComp]);
        if(structure[iComp].nfft_flags & PRE_ONE_PSI)
            nfft_precompute_one_psi(&structure[iComp]);
        
    }
    //*/
    return;
    
}

double Epot_reci(double X[][DIM], double D[][DIM], const int Ncol, const double Vol){
    
    // Setting the nodes
    
    for (int iCol=0; iCol<Ncol; iCol++){
        for(int iComp=0; iComp<DIM; iComp++){
            for (int iDim=0; iDim<DIM; iDim++){            
                
                potential[iComp].x[DIM*iCol+iDim] = X[iCol][iDim];            
                structure[iComp].x[DIM*iCol+iDim] = X[iCol][iDim];            
                
            }
        }
    }
    
    // Precompute $\psi$
    
    for(int iComp=0; iComp<DIM; iComp++){
        
        if(potential[iComp].nfft_flags & PRE_ONE_PSI)
            nfft_precompute_one_psi(&potential[iComp]);
        if(structure[iComp].nfft_flags & PRE_ONE_PSI)
            nfft_precompute_one_psi(&structure[iComp]);
        
    }
    
     // Setting the function : structure
    
    for (int iCol=0; iCol<Ncol; iCol++){
        for (int iComp=0; iComp<DIM; iComp++){
            structure[iComp].f[iCol] = D[iCol][iComp];
        }
    }
    
    // Doing the transform : structure
    
    for (int iComp=0; iComp<DIM; iComp++){
        nfft_adjoint_direct(&structure[iComp]);
    }
    
    // Setting the function potential : potential
    
    int ikx, iky, ik;
    double complex scalarProductCplx;
    double Epot_tabulated;
    
    for (int kx=-Nx; kx<Nx; kx++){
        
        ikx = (kx + Nx)*(2*Ny * 2*Nz);
        
        for (int ky=-Ny; ky<Ny; ky++){
            
            iky = (ky + Ny)*2*Nz;
            
            for (int kz=-Nz; kz<Nz; kz++){
                
                ik = ikx + iky + (kz + Nz);
                
                scalarProductCplx  = (double)kx * structure[0].f_hat[ik] +
                                     (double)ky * structure[1].f_hat[ik] +
                                     (double)kz * structure[2].f_hat[ik];
                
                for (int iComp=0; iComp<DIM; iComp++){
                    potential[iComp].f_hat[ik] = scalarProductCplx;
                }
                
                Epot_tabulated = Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz];
                
                potential[0].f_hat[ik] *= (double)kx * Epot_tabulated;
                potential[1].f_hat[ik] *= (double)ky * Epot_tabulated;
                potential[2].f_hat[ik] *= (double)kz * Epot_tabulated;
            
            }            
        }        
    }
    
    // Doing the transform :potential
    
    for (int iComp=0; iComp<DIM; iComp++){
        nfft_trafo_direct(&potential[iComp]);
    }
    
    // Result
    
    double Epot = 0.;
    
    for (int jCol=0; jCol<Ncol; jCol++){
    
        scalarProductCplx = 0.;
        for (int iComp=0; iComp<DIM; iComp++){
            scalarProductCplx += D[jCol][iComp]*potential[iComp].f[jCol];
        }
    
        Epot += creal(scalarProductCplx);
    
    }
    
    return 2.*PI/Vol * Epot;
    
}

void snapShot(const int Ncol){
    
    FILE *fp0 = fopen("C_snapX_0.out", "w");
    FILE *fp1 = fopen("C_snapX_1.out", "w");
    FILE *fp2 = fopen("C_snapX_2.out", "w");
 
    for (int iCol=0; iCol<Ncol; iCol++){
        for (int iDim=0; iDim<DIM; iDim++){
            fprintf(fp0, "%g ", structure[0].x[DIM*iCol+iDim]);
            fprintf(fp1, "%g ", structure[1].x[DIM*iCol+iDim]);
            fprintf(fp2, "%g ", structure[2].x[DIM*iCol+iDim]);
        }
        fprintf(fp0, "\n");
        fprintf(fp1, "\n");
        fprintf(fp2, "\n");
    }
    
    fclose(fp0);
    fclose(fp1);
    fclose(fp2);
    
    return;
}
