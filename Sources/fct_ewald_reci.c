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
double Epot_reci(double X[][DIM], double D[][DIM], const int Ncol, const double Vol);

// Initialisation : phi and S

void Epot_reci_nfft_init(const int Ncol){
    
    for(int iComp=0; iComp<DIM; iComp++){
        
        nfft_init_3d(&potential[iComp], 2*Nx, 2*Ny, 2*Nz, Ncol);
        nfft_init_3d(&structure[iComp], 2*Nx, 2*Ny, 2*Nz, Ncol);
        
    }
    
}


// Finalisation

void Epot_reci_nfft_finalize(){
 
    for (int iComp=0; iComp<DIM; iComp++){
        nfft_finalize(&structure[iComp]);
        nfft_finalize(&potential[iComp]);
    }
    
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

                    kOverLsqr = kxOverL*kxOverL + kyOverL*kyOverL + 
                        kzOverL*kzOverL;
                    Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz] = 
                        exp(-PI*PI/alphaSqr*kOverLsqr);
                    Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz] /= kOverLsqr;
                    nb_k++;
                
                }
                
            }
            
        }
        
    }
    
    Epot_reci_tab[Nx][Ny][Nz] = 0.;
    
    printf(" Ewald Summation : Number of wave vectors : %d\n", nb_k);    
    
}

double Epot_reci(double X[][DIM], double D[][DIM], const int Ncol, const double Vol){
    
    // Setting the nodes
    
    for (int iCol=0; iCol<Ncol; iCol++){    
        for (int iDim=0; iDim<DIM; iDim++){
            for(int iComp=0; iComp<DIM; iComp++){     
                
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
    double Epot_reci_tabulated;
    
    for (int kx=-Nx; kx<Nx; kx++){
        
        ikx = (kx + Nx)*(2*Ny * 2*Nz);
        
        for (int ky=-Ny; ky<Ny; ky++){
            
            iky = (ky + Ny)*2*Nz;
            
            for (int kz=-Nz; kz<Nz; kz++){
                
                ik = ikx + iky + (kz + Nz);
                
                scalarProductCplx  = (double)kx * structure[0].f_hat[ik];
                scalarProductCplx += (double)ky * structure[1].f_hat[ik];
                scalarProductCplx += (double)kz * structure[2].f_hat[ik];
                
                for (int iComp=0; iComp<DIM; iComp++){
                    potential[iComp].f_hat[ik] = scalarProductCplx;
                }
                
                Epot_reci_tabulated = Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz];
                
                potential[0].f_hat[ik] *= (double)kx * Epot_reci_tabulated;
                potential[1].f_hat[ik] *= (double)ky * Epot_reci_tabulated;
                potential[2].f_hat[ik] *= (double)kz * Epot_reci_tabulated;
            
            }            
        }        
    }
    
    // Doing the transform :potential
    
    for (int iComp=0; iComp<DIM; iComp++){
        nfft_trafo_direct(&potential[iComp]);
    }
    
    // Result
    
    double ePot_reci = 0.;
    
    for (int jCol=0; jCol<Ncol; jCol++){
    
        scalarProductCplx = 0.;
        for (int iComp=0; iComp<DIM; iComp++){
            scalarProductCplx += D[jCol][iComp]*potential[iComp].f[jCol];
        }
    
        ePot_reci += creal(scalarProductCplx);
    
    }
    
    ePot_reci *= 2.*PI/Vol;
    
    return ePot_reci;
    
}
