/* \brief Non-uniform Fast Fourier Transform + delta of energy
 * 
 *  Uses NFFT3 library based on FFTW3.
 */

#include <stdio.h>
#include <stdlib.h>
#include <tgmath.h>
#include "nfft3.h"
#include "data_reci.h"

static nfft_plan structure[DIM]; //!< ``Structure factor'' \f[ \vec{S} \f]
static nfft_plan potential[DIM]; //!< ``Potential'' \f[ \vec{\phi} \f]
static double Epot_reci_tab[2*Nx][2*Ny][2*Nz]; /*!< \f[ f(\alpha, \vec{k}) = \frac{e^{-\frac{\pi^2}
{\alpha^2} \sum_d \frac{k_d^2}{L_d}}}{\sum_d \frac{k_d^2}{L_d}} \f]
*/

void Epot_reci_nfft_init(const int Ncol);
void Epot_reci_nfft_finalize();
void Epot_reci_init(const double Lsize[DIM], const double alpha);

static void Epot_reci_fourier_new(const double xNew[DIM]);
static void Epot_reci_fourier_old(const double xOld[DIM]);

double Epot_reci_move(const int lCol, const double xNew[DIM], const double Vol);
void Epot_reci_updateX(const int lCol, const double xNew[DIM]);

double Epot_reci_rotate(const int lCol, const double mNew[DIM], const double Vol);
void Epot_reci_updateM(const int lCol, const double mNew[DIM]);

double Epot_reci(double X[][DIM], double D[][DIM], const int Ncol, const double Vol);
void snapShot(const int Ncol);

//!> Initialisation : ``Potential'' and ``Structure factor''

void Epot_reci_nfft_init(const int Ncol){
    
    for(int iComp=0; iComp<DIM; iComp++){
        
        nfft_init_3d(&potential[iComp], 2*Nx, 2*Ny, 2*Nz, Ncol);
        nfft_init_3d(&structure[iComp], 2*Nx, 2*Ny, 2*Nz, Ncol);
        
    }
    
    return;
    
}


//!> Finalisation

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
    
    FILE *fpX = fopen("C_report.txt", "w");
        fprintf(fpX, " Ewald Summation : Number of wave vectors : %d\n", nb_k);
        fprintf(fpX, "     Nx = %d\n", Nx);
        fprintf(fpX, "     Ny = %d\n", Nx);
        fprintf(fpX, "     Nz = %d\n", Nx);
    fclose(fpX);
    
    return;
    
}

void Epot_reci_fourier(const double xNew[DIM], complex exp_IkxCol_x[2*Nx], 
                                               complex exp_IkxCol_y[2*Ny],
                                               complex exp_IkxCol_z[2*Nz]){

    // x
    
    exp_IkxCol_x[Nx+0] = 1.;
    exp_IkxCol_x[Nx+1] = exp(I*2.*PI * 1.*xNew[0]);
    exp_IkxCol_x[Nx-1] = conj(exp_IkxCol_x[Nx+1]);
    
    for (int kx=2; kx<Nx; kx++){
        
        exp_IkxCol_x[Nx+kx] = exp_IkxCol_x[Nx+kx-1] * exp_IkxCol_x[Nx+1];
        exp_IkxCol_x[Nx-kx] = conj(exp_IkxCol_x[Nx+kx]);
        
    }    
    
    exp_IkxCol_x[Nx-Nx] = exp_IkxCol_x[Nx-(Nx-1)] * exp_IkxCol_x[Nx-1];
    
    // y
    
    exp_IkxCol_y[Ny+0] = 1.;
    exp_IkxCol_y[Ny+1] = exp(I*2.*PI * 1.*xNew[1]);
    exp_IkxCol_y[Ny-1] = conj(exp_IkxCol_y[Ny+1]);
    
    for (int ky=2; ky<Ny; ky++){
        
        exp_IkxCol_y[Ny+ky] = exp_IkxCol_y[Ny+ky-1] * exp_IkxCol_y[Ny+1];
        exp_IkxCol_y[Ny-ky] = conj(exp_IkxCol_y[Ny+ky]);
        
    }
    exp_IkxCol_y[Ny-Ny] = exp_IkxCol_y[Ny-(Ny-1)] * exp_IkxCol_y[Ny-1];
            
    // z
    
    exp_IkxCol_z[Nz+0] = 1.;
    exp_IkxCol_z[Nz+1] = exp(I*2.*PI * 1.*xNew[2]);
    exp_IkxCol_z[Nz-1] = conj(exp_IkxCol_z[Nz+1]);
    
    for (int kz=2; kz<Nz; kz++){
        
        exp_IkxCol_z[Nz+kz] = exp_IkxCol_z[Nz+kz-1] * exp_IkxCol_z[Nz+1];
        exp_IkxCol_z[Nz-kz] = conj(exp_IkxCol_z[Nz+kz]);
        
    }
    
    exp_IkxCol_z[Nz-Nz] = exp_IkxCol_z[Nz-(Nz-1)] * exp_IkxCol_z[Nz-1];

}

// Move --------------------------------------------------------------------------------------------

/*!> Difference of Energy \f[ \Delta U = \frac{2\pi}{V} \sum_{\vec{k} \neq 0} \Delta M^2 
 * f(\alpha, \vec{k}) \f]
 * \f[
 *  \Delta M^2 = 2\Re[
 *      (\vec{\mu}_l\cdot\vec{k}) (e^{i(\vec{k}\cdot\vec{x}^\prime_l)} -e^{i(\vec{k}\cdot\vec{x}_l)}
 *      (\vec{k}\cdot\vec{S}_l)
 *  ]
 * \f]
 * \f[  \vec{S}_l = \sum_{i \neq l} \vec{\mu}_i e^{i(\vec{k}\cdot\vec{x}_i)}
 * \f]
 * Implementation :
 * \f[ 
 *  \Delta M^2 = 2(\vec{\mu_l}\cdot\vec{k}) 
 *  [ \cos(\vec{k}\cdot\vec{x}^\prime_l) - \cos(\vec{k}\cdot\vec{x})]
 *  [\Re{(\vec{k}\cdot\vec{S})} - (\vec{k}\cdot\vec{\mu}_l) cos(\vec{k}\cdot\vec{x}_l)] -
 *  [-\sin(\vec{k}\cdot\vec{x}^\prime_l) + \sin(\vec{k}\cdot\vec{x})]
 *  [\Im{(\vec{k}\cdot\vec{S})} - (\vec{k}\cdot\vec{\mu}_l) sin(\vec{k}\cdot\vec{x}_l)]
 * \f]
 */

double Epot_reci_move(const int lCol, const double xNew[DIM], const double Vol){

    double Epot;
    int ik;
    double k_dot_mOld;
    double complex k_dot_structure;
    double complex exp_IkxNew, exp_IkxOld;
    double realPart1, realPart2;
    
    double cos_kxNew, sin_kxNew;
    double cos_kxOld, sin_kxOld;
    
    double xOld[DIM];
    
    for (int iDim=0; iDim<DIM; iDim++){
        xOld[iDim] = structure[0].x[DIM*lCol+iDim];
    }
    
    double complex exp_IkxNew_x[2*Nx];
    double complex exp_IkxNew_y[2*Ny];
    double complex exp_IkxNew_z[2*Nz];
    
    double complex exp_IkxOld_x[2*Nx];
    double complex exp_IkxOld_y[2*Ny];
    double complex exp_IkxOld_z[2*Nz];
    
    Epot_reci_fourier(xNew, exp_IkxNew_x, exp_IkxNew_y, exp_IkxNew_z);
    Epot_reci_fourier(xOld, exp_IkxOld_x, exp_IkxOld_y, exp_IkxOld_z);
        
    Epot = 0.;
    ik = 0;
    
    for (int kx=-Nx; kx<Nx; kx++){
        
        for (int ky=-Ny; ky<Ny; ky++){
            
            for (int kz=-Nz; kz<Nz; kz++){
                
                k_dot_mOld = (double)kx * creal(structure[0].f[lCol]) +
                             (double)ky * creal(structure[1].f[lCol]) +
                             (double)kz * creal(structure[2].f[lCol]);
                             
                k_dot_structure = (double)kx * structure[0].f_hat[ik] +
                                  (double)ky * structure[1].f_hat[ik] +
                                  (double)kz * structure[2].f_hat[ik];
                ik++;
                
                exp_IkxNew = exp_IkxNew_x[kx+Nx] * exp_IkxNew_y[ky+Ny] * exp_IkxNew_z[kz+Nz];
                cos_kxNew = creal(exp_IkxNew);
                sin_kxNew = cimag(exp_IkxNew);
                
                exp_IkxOld = exp_IkxOld_x[kx+Nx] * exp_IkxOld_y[ky+Ny] * exp_IkxOld_z[kz+Nz];
                cos_kxOld = creal(exp_IkxOld);
                sin_kxOld = cimag(exp_IkxOld);
                
                realPart1 = cos_kxNew - cos_kxOld;                
                realPart1*= creal(k_dot_structure) - k_dot_mOld * cos_kxOld;
                
                realPart2 =-sin_kxNew + sin_kxOld;
                realPart2*= cimag(k_dot_structure) - k_dot_mOld * sin_kxOld;
                
                Epot += 2.*k_dot_mOld * (realPart1 - realPart2) * Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz];
            
            }            
        }        
    }
    
    Epot *= 2.*PI/Vol;
    
    return Epot;
    
}

/*!> Update position -> update the ``structure factor''
 *  \f[
 *      \Delta \vec{S} = \vec{\mu}_l
 *      (e^{i(\vec{k}\cdot\vec{x}^\prime_l)} -e^{i(\vec{k}\cdot\vec{x}_l)}) 
 *  \f]
 */

void Epot_reci_updateX(const int lCol, const double xNew[DIM]){
    
    double xOld[DIM];
 
    for (int iDim=0; iDim<DIM; iDim++){
        
        xOld[iDim] = structure[0].x[DIM*lCol+iDim];
        
        for(int iComp=0; iComp<DIM; iComp++){
            
            structure[iComp].x[DIM*lCol+iDim] = xNew[iDim];
            
        }
        
    }
    
    int ik;    
    double complex exp_IkxNew;
    double complex exp_IkxOld;
    
    double complex exp_IkxNew_x[2*Nx];
    double complex exp_IkxNew_y[2*Ny];
    double complex exp_IkxNew_z[2*Nz];
    
    double complex exp_IkxOld_x[2*Nx];
    double complex exp_IkxOld_y[2*Ny];
    double complex exp_IkxOld_z[2*Nz];
    
    Epot_reci_fourier(xNew, exp_IkxNew_x, exp_IkxNew_y, exp_IkxNew_z);
    Epot_reci_fourier(xOld, exp_IkxOld_x, exp_IkxOld_y, exp_IkxOld_z);
    
    ik = 0;
    
    for (int kx=-Nx; kx<Nx; kx++){
        
        for (int ky=-Ny; ky<Ny; ky++){
            
            for (int kz=-Nz; kz<Nz; kz++){
    
                exp_IkxNew = exp_IkxNew_x[kx+Nx] * exp_IkxNew_y[ky+Ny] * exp_IkxNew_z[kz+Nz];                                
                exp_IkxOld = exp_IkxOld_x[kx+Nx] * exp_IkxOld_y[ky+Ny] * exp_IkxOld_z[kz+Nz];

                for(int iComp=0; iComp<DIM; iComp++){
                    
                    structure[iComp].f_hat[ik] += creal(structure[iComp].f[lCol]) *
                                                  (exp_IkxNew - exp_IkxOld);
                                                
                }
                ik++;
    
            }
            
        }
        
    }
    
    return;
    
}

// -------------------------------------------------------------------------------------------------

// Rotate ------------------------------------------------------------------------------------------

double Epot_reci_rotate(const int lCol, const double mNew[DIM], const double Vol){

    double Epot, Epot_k;
    int ik;
    double k_dot_mNew, k_dot_mOld;
    double complex k_dot_structure;
    double complex exp_IkxOld;
    double realPart;
    double cos_kxOld, sin_kxOld;
    
    double xOld[DIM];
    
    for (int iDim=0; iDim<DIM; iDim++){
        xOld[iDim] = structure[0].x[DIM*lCol+iDim];
    }
    
    double complex exp_IkxOld_x[2*Nx];
    double complex exp_IkxOld_y[2*Ny];
    double complex exp_IkxOld_z[2*Nz];

    Epot_reci_fourier(xOld, exp_IkxOld_x, exp_IkxOld_y, exp_IkxOld_z);
    
    Epot = 0.;
    ik = 0;
    
    for (int kx=-Nx; kx<Nx; kx++){
        
        for (int ky=-Ny; ky<Ny; ky++){
            
            for (int kz=-Nz; kz<Nz; kz++){              

                
                k_dot_mNew = (double)kx * mNew[0] +
                             (double)ky * mNew[1] +
                             (double)kz * mNew[2];               
                
                k_dot_mOld = (double)kx * creal(structure[0].f[lCol]) +
                             (double)ky * creal(structure[1].f[lCol]) +
                             (double)kz * creal(structure[2].f[lCol]);                             
                             
                k_dot_structure = (double)kx * structure[0].f_hat[ik] +
                                  (double)ky * structure[1].f_hat[ik] +
                                  (double)kz * structure[2].f_hat[ik];
                ik++;
                
                exp_IkxOld = exp_IkxOld_x[kx+Nx] * exp_IkxOld_y[ky+Ny] * exp_IkxOld_z[kz+Nz];
                cos_kxOld = creal(exp_IkxOld);
                sin_kxOld = cimag(exp_IkxOld);
                
                realPart = cos_kxOld * (creal(k_dot_structure) - k_dot_mOld * cos_kxOld);
                realPart+= sin_kxOld * (cimag(k_dot_structure) - k_dot_mOld * sin_kxOld);

                Epot_k = k_dot_mNew*k_dot_mNew - k_dot_mOld*k_dot_mOld;
                Epot_k+= 2.*(k_dot_mNew - k_dot_mOld) * realPart;
                Epot_k*= Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz];
                Epot += Epot_k;
                
            }            
        }        
    }
    
    Epot *= 2.*PI/Vol;
    
    return Epot;
    
}

void Epot_reci_updateM(const int lCol, const double mNew[DIM]){

    double mOld[DIM];
        
    for(int iDim=0; iDim<DIM; iDim++){
        
        mOld[iDim] = creal(structure[iDim].f[lCol]);
        structure[iDim].f[lCol] = mNew[iDim];
        
    }

    int ik;
    double complex exp_IkxOld;
    double xOld[DIM];
    
    for (int iDim=0; iDim<DIM; iDim++){
        xOld[iDim] = structure[0].x[DIM*lCol+iDim];
    }
    
    double complex exp_IkxOld_x[2*Nx];
    double complex exp_IkxOld_y[2*Ny];
    double complex exp_IkxOld_z[2*Nz];

    Epot_reci_fourier(xOld, exp_IkxOld_x, exp_IkxOld_y, exp_IkxOld_z);
    
    ik = 0;
    
    for (int kx=-Nx; kx<Nx; kx++){
        
        for (int ky=-Ny; ky<Ny; ky++){
            
            for (int kz=-Nz; kz<Nz; kz++){
                
                exp_IkxOld = exp_IkxOld_x[kx+Nx] * exp_IkxOld_y[ky+Ny] * exp_IkxOld_z[kz+Nz];

                for(int iComp=0; iComp<DIM; iComp++){
                    
                    structure[iComp].f_hat[ik] += (mNew[iComp] - mOld[iComp]) * exp_IkxOld;
                                                
                }
                ik++;
    
            }
            
        }
        
    }
    
    return;

}

// -------------------------------------------------------------------------------------------------

// Test particle -----------------------------------------------------------------------------------

double Epot_reci_test(const double xTest[DIM], const double mTest[DIM], const double Vol){

    double Epot;
    int ik;
    double k_dot_mTest;
    double complex k_dot_structure;
    double complex exp_IkxTest;
    double realPart;
    double cos_kxTest, sin_kxTest;
    
    double complex exp_IkxTest_x[2*Nx];
    double complex exp_IkxTest_y[2*Ny];
    double complex exp_IkxTest_z[2*Nz];

    Epot_reci_fourier(xTest, exp_IkxTest_x, exp_IkxTest_y, exp_IkxTest_z);
        
    Epot = 0.;
    ik = 0;
    
    for (int kx=-Nx; kx<Nx; kx++){
        
        for (int ky=-Ny; ky<Ny; ky++){
            
            for (int kz=-Nz; kz<Nz; kz++){
            
                k_dot_mTest = (double)kx * mTest[0] +
                              (double)ky * mTest[1] +
                              (double)kz * mTest[2];  
            
                k_dot_structure = (double)kx * structure[0].f_hat[ik] +
                                  (double)ky * structure[1].f_hat[ik] +
                                  (double)kz * structure[2].f_hat[ik];
                ik++;
                
                exp_IkxTest = exp_IkxTest_x[kx+Nx] * exp_IkxTest_y[ky+Ny] * exp_IkxTest_z[kz+Nz];
                cos_kxTest = creal(exp_IkxTest);
                sin_kxTest =-cimag(exp_IkxTest);
                
                realPart = creal(k_dot_structure) * cos_kxTest;
                realPart+= cimag(k_dot_structure) * sin_kxTest;
                
                Epot += k_dot_mTest * (k_dot_mTest + 2.*realPart) * 
                        Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz];
            
            }
            
        }
        
    }
    
    Epot *= 2.*PI/Vol;
    
    return Epot;

}

// -------------------------------------------------------------------------------------------------

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
        nfft_adjoint(&structure[iComp]);
    }
    
    // Setting the function potential : potential
    
    int ik;
    double complex k_dot_structure;
    double Epot_tabulated;
    
    ik = 0;
    
    for (int kx=-Nx; kx<Nx; kx++){
        
        for (int ky=-Ny; ky<Ny; ky++){
            
            for (int kz=-Nz; kz<Nz; kz++){
                
                k_dot_structure  = (double)kx * structure[0].f_hat[ik] +
                                   (double)ky * structure[1].f_hat[ik] +
                                   (double)kz * structure[2].f_hat[ik];
                
                for (int iComp=0; iComp<DIM; iComp++){
                    potential[iComp].f_hat[ik] = k_dot_structure;
                }
                
                Epot_tabulated = Epot_reci_tab[kx+Nx][ky+Ny][kz+Nz];
                
                potential[0].f_hat[ik] *= (double)kx * Epot_tabulated;
                potential[1].f_hat[ik] *= (double)ky * Epot_tabulated;
                potential[2].f_hat[ik] *= (double)kz * Epot_tabulated;
                
                ik ++;
            
            }            
        }        
    }
    
    // Doing the transform :potential
    
    for (int iComp=0; iComp<DIM; iComp++){
        nfft_trafo(&potential[iComp]);
    }
    
    // Result
    
    double Epot = 0.;
    double complex mCol_dot_potential;
    
    for (int jCol=0; jCol<Ncol; jCol++){
    
        mCol_dot_potential = 0.;
        for (int iComp=0; iComp<DIM; iComp++){
            mCol_dot_potential += D[jCol][iComp]*potential[iComp].f[jCol];
        }
    
        Epot += creal(mCol_dot_potential);
    
    }
    
    Epot*= 2.*PI/Vol;
    
    return Epot;
    
}

void snapShot(const int Ncol){
    
    FILE *fpX = fopen("C_snapX.out", "w");
    FILE *fpM = fopen("C_snapM.out", "w");
 
    for (int iCol=0; iCol<Ncol; iCol++){
        for (int iDim=0; iDim<DIM; iDim++){
            fprintf(fpX, "%g ", structure[0].x[DIM*iCol+iDim]);
            fprintf(fpM, "%g ", creal(structure[iDim].f[iCol]));
        }
        fprintf(fpX, "\n");
        fprintf(fpM, "\n");
    }
    
    fclose(fpX);
    fclose(fpM);
    
    return;
}
