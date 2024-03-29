//
//  main.cpp
//  transportBuddy for 3D cuprates
//  by Simon Verret, Universite de Sherbrooke
//
//  code partially adapted from afmCond.c
//  by Maxime Charlebois, Universite de Sherbrooke
//  under MIT license (september 2017)
//

#define _USE_MATH_DEFINES
#include <math.h>
#include <complex.h>
#include <sys/stat.h>
#include <stdio.h>

#include <stdlib.h>
#include <string.h>



// FUNCTIONS FOR READING FILES

void readDouble(FILE * file, char * name,  double * value) {
    rewind(file);
    char tempbuff[200];
    while(!feof(file))
    {
        if (fgets(tempbuff,200,file))
        {
            char tmpstr1[50];
            char tmpstr2[50];
            sscanf(tempbuff, "%49s %49s\n", tmpstr1, tmpstr2);
            if (strcmp(tmpstr1,name)==0) { *value = atof(tmpstr2); return;}
        }
    }
    printf("\nWARNING: cannot find the %s parameter in 'model.dat'", name);
}

void readInt(FILE * file, char * name,  int * value) {
    rewind(file);
    char tempbuff[200];
    while(!feof(file))
    {
        if (fgets(tempbuff,200,file))
        {
            char tmpstr1[50];
            char tmpstr2[50];
            sscanf(tempbuff, "%49s %49s\n", tmpstr1, tmpstr2);
            if (strcmp(tmpstr1,name)==0) { *value = atoi(tmpstr2); return;}
        }
    }
    printf("\nWARNING: cannot find the %s parameter in 'model.dat'", name);
}

void removeSubstring(char *s,const char *toremove)
{
  while( s == strstr(s,toremove) )
    memmove(s,s+strlen(toremove),1+strlen(s+strlen(toremove)));
}



//// MAIN

int main(int argc, const char * argv[]) {
    
    
    //// SETTING PARAMETERS
    double t;
    double tp;
    double tpp;
    double tppp;
    double tz;
    
    double ETA=1e-5;
    double ETAell=0.0;
    double ETAdos=0.0;
    double ETAk=0.0;
    double ETAaFL=0.0;
    double ETAbFL=0.0;
    double ETAaPL=0.0;
    double ETAw=0.0;
    double ETAkw=0.0;
    int nOmega;
    int nK;
    int nKz;
    int nMu; double muMin; double muMax;
    int nT; double Tmin; double Tmax;
    int logT=0; //bool
    double amplitudeCutoff;
    
    double etaNorm=0;
    
    //// READING PARAMETERS FROM "model.dat" FILE
    
    char model[800] = "";
    if (argc >= 2) { strcpy(model,argv[1]);}
    else { strcpy(model,"model.dat");}
    FILE * file = fopen(model,"rt");
    if(file == NULL) {printf("file %s not found", model); exit(1);}
    printf("reading parameters from %s\n", model) ;

    removeSubstring(model,"model");
    
    readDouble(file, "t",               &t);
    readDouble(file, "tp",              &tp);
    readDouble(file, "tpp",             &tpp);
    readDouble(file, "tppp",            &tppp);
    readDouble(file, "tz",              &tz);
    
    readDouble(file, "muMin",           &muMin);
    readDouble(file, "muMax",           &muMax);
    readDouble(file, "Tmin",            &Tmin);
    readDouble(file, "Tmax",            &Tmax);
    readDouble(file, "ETA",             &ETA);
    etaNorm=ETA;
    readDouble(file, "ETAell",          &ETAell);
    if (ETAell>etaNorm) etaNorm=ETAell;
    readDouble(file, "ETAdos",          &ETAdos);
    if (ETAdos>etaNorm) etaNorm=ETAdos;
    readDouble(file, "ETAk",            &ETAk);
    if (  ETAk>etaNorm) etaNorm=ETAk;
    readDouble(file, "ETAaFL",          &ETAaFL);
    readDouble(file, "ETAbFL",          &ETAbFL);
    readDouble(file, "ETAaPL",          &ETAaPL);
    readDouble(file, "ETAw",          &ETAw);
    readDouble(file, "ETAkw",          &ETAkw);
    readDouble(file, "amplitudeCutoff", &amplitudeCutoff);
    
    readInt(file, "nMu",    &nMu);
    readInt(file, "logT",   &logT); //to act as a bool
    readInt(file, "nT",     &nT);
    readInt(file, "nK",     &nK);
    readInt(file, "nKz",    &nKz);
    readInt(file, "nOmega", &nOmega);
    
    fclose(file);
    

    //// COMPUTING

    printf("transportBuddy starting\n\n");

    char nameFileT[800] = "transport_vs_T";
    strcat(nameFileT, model);
    FILE *fileOut = fopen(nameFileT,"w");

    //// PRE-CALCULATE TEMPERATURE AND FERMI DISTRIBUTIONS
    double T[nT];
    double beta[nT];
    double omega[nT][nOmega];
    double energyCutoff[nT];
    double dfermiDirac_dw[nT][nOmega];
    double dfermiDirac_dT[nT][nOmega];
    
    // loop on temperature for precalculation only
    int nn=0; for(nn=0; nn<nT; nn++) {
        T[nn] = Tmin;
        if ((Tmin != Tmax) && (nT != 1)) {
            if(logT==1) T[nn] *= exp(nn*log(Tmax/Tmin)/(nT-1));
            else T[nn] += nn*(Tmax-Tmin)/(nT-1);
        }
        beta[nn] = 1./T[nn];

        energyCutoff[nn] = 2.*2.*acosh(0.25*sqrt(beta[nn]/amplitudeCutoff)) /beta[nn];
        int n=0; for(n=0; n<nOmega; n++)
        {
            omega[nn][n]= -energyCutoff[nn] + 2.*energyCutoff[nn]*n/(nOmega-1);
            double expBw = exp(beta[nn]*omega[nn][n]);
            dfermiDirac_dw[nn][n]= -beta[nn]*expBw/((expBw+1.)*(expBw+1.));
            dfermiDirac_dT[nn][n]= -beta[nn]*omega[nn][n] * dfermiDirac_dw[nn][n];
        }
    }
    
    //// PRE-CALCULATE SINES AND COSINES
    double sink[nK];
    double sin2k[nK];
    double sink_2[nK];
    double cosk[nK];
    double cos2k[nK];
    double cosk_2[nK];

    // loop on k for precalculation only
    // only on 1/4th of the Brillouin zone (assumes fourfold symmetry)
    int ii=0; for(ii=0; ii<nK; ii++){
        double k = M_PI*(ii*1.0/nK);
        sink[ii] = sin(k); sin2k[ii] = sin(2.*k); sink_2[ii] = sin(k/2.);
        cosk[ii] = cos(k); cos2k[ii] = cos(2.*k); cosk_2[ii] = cos(k/2.);
    }
    //
    double coskz[nKz]; double coskz_2[nKz]; double sinkz_2[nKz];
    int kk=0; for(kk=0; kk<nKz; kk++){
        double kz = M_PI*(kk*2.0/(nKz)); // 0 to 2Pi assumes parity (periodicity in z is 4Pi because of the tetragonal structure)
        coskz[kk] = cos(kz); coskz_2[kk] = cos(kz/2.); sinkz_2[kk] = sin(kz/2.);
    }
    
    //// RESULTS TO FILL
    double dos[nMu];
    double density_0[nMu];
    double sigmaxx_0[nMu];
    double sigmazz_0[nMu];
    double sigmaxy_0[nMu];
    double Cv[nMu][nT];
    double density[nMu][nT];
    double sigma_xx[nMu][nT];
    double sigma_zz[nMu][nT];
    double sigma_xy[nMu][nT];
    double alpha_xx[nMu][nT];
    double alpha_zz[nMu][nT];
    double alpha_xy[nMu][nT];
    double beta_xx[nMu][nT];
    double beta_zz[nMu][nT];
    double beta_xy[nMu][nT];
    
    //// LOOP ON mu (chemical potential)
    for(int iMu=0; iMu<nMu; iMu++)
    {
        double mu = muMin;
        if (muMin!=muMax) mu += iMu*(muMax-muMin)/(nMu-1);
        
        //// INITIALIZE RESULTS TO 0
        dos[iMu] = 0.;
        density_0[iMu] = 0.;
        sigmaxx_0[iMu] = 0.;
        sigmazz_0[iMu] = 0.;
        sigmaxy_0[iMu] = 0.;

        /// INITIALIZE RESULTS vs TEMPERATURE TO 0
        int nn=0; for(nn=0; nn<nT; nn++)
        {
            Cv[iMu][nn]=0.;
            density[iMu][nn]=0.;
            sigma_xx[iMu][nn]=0.;
            sigma_zz[iMu][nn]=0.;
            sigma_xy[iMu][nn]=0.;
            alpha_xx[iMu][nn]=0.;
            alpha_zz[iMu][nn]=0.;
            alpha_xy[iMu][nn]=0.;
            beta_xx[iMu][nn]=0.;
            beta_zz[iMu][nn]=0.;
            beta_xy[iMu][nn]=0.;
        }
        
        //// INTEGRAL ON k
        for(int ii=0; ii<nK; ii++) {
            printf("muu= %i / %i, k = %i / %i \n",iMu,nMu,ii,nK);//fflush(stdout);
            for(int jj=0; jj<nK; jj++) {
                
                double epsilon_k   = -mu;
                epsilon_k         += -2.*t    * (cosk[ii] + cosk[jj]);
                epsilon_k         += -4.*tp   *  cosk[ii]*cosk[jj];
                epsilon_k         += -2.*tpp  * (cos2k[ii] + cos2k[jj]);
                epsilon_k         += -4.*tppp * (cos2k[ii]*cosk[jj] + cosk[ii]*cos2k[jj]);
                double depsilon_dkx = 0;
                depsilon_dkx      +=  2.*t    *  sink[ii];
                depsilon_dkx      +=  4.*tp   *  sink[ii]*cosk[jj];
                depsilon_dkx      +=  4.*tpp  *  sin2k[ii];
                depsilon_dkx      +=     tppp * (8.*sin2k[ii]*cosk[jj] + 4.*sink[ii]*cos2k[jj]);
                double depsilon_dkx_dky = 0;
                depsilon_dkx_dky  += -4.*tp   *  sink[ii]*sink[jj];
                depsilon_dkx_dky  += -8.*tppp * (sin2k[ii]*sink[jj] + sink[ii]*sin2k[jj]);
                double depsilon_dky = 0;
                depsilon_dky      +=  2.*t    *  sink[jj];
                depsilon_dky      +=  4.*tp   *  sink[jj]*cosk[ii];
                depsilon_dky      +=  4.*tpp  *  sin2k[jj];
                depsilon_dky      +=     tppp * (8.*sin2k[jj]*cosk[ii] + 4.*sink[jj]*cos2k[ii]);
                double depsilon_dky_dky = 0;
                depsilon_dky_dky  +=  2.*t    *  cosk[jj];
                depsilon_dky_dky  +=  4.*tp   *  cosk[jj]*cosk[ii];
                depsilon_dky_dky  +=  8.*tpp  *  cos2k[jj];
                depsilon_dky_dky  +=     tppp * (16.*cos2k[jj]*cosk[ii] + 8.*cosk[jj]*cos2k[ii]);
                
                // This will not work outside the first quadrant of the Brillouin zone
                double cos_phi12 = 0.0;
                if (ii > 0 & jj > 0 & (ETAk > 0 || ETAkw > 0)) {
                    double cos_phi = cos(atan(2*jj/ii));
                    cos_phi12 = pow(cos_phi,12);
                }

                for(int kk=0; kk<nKz; kk++) {
                    
                    // See Markiewicz 2005 for the 3D tight-binding model // derivatives from mathematica
                    double epsilonz_k         = -2*tz*coskz_2[kk]*(cosk[ii]-cosk[jj])*(cosk[ii]-cosk[jj])*cosk_2[ii]*cosk_2[jj];
                    double depsilonz_dkx      = tz*cosk_2[jj]*(-4-5*cosk[ii]+cosk[jj])*(-cosk[ii]+cosk[jj])*coskz_2[kk]*sink_2[ii];
                    double depsilonz_dky      = tz*cosk_2[ii]*(-4+cosk[ii]-5*cosk[jj])*(cosk[ii]-cosk[jj])*coskz_2[kk]*sink_2[jj];
                    double depsilonz_dkx_dky  = tz*(1/4.)*(6-5*cos2k[ii]+16*cosk[jj]+4*cosk[ii]*(4+9*cosk[jj])-5*cos2k[jj])*coskz_2[kk]*sink_2[ii]*sink_2[jj];
                    double depsilonz_dky_dky  = tz*(1/4.)*cosk_2[ii]*cosk_2[jj]*(10+16*cosk[ii]+cos2k[ii]-4*(4+9*cosk[ii])*cosk[jj]+25*cos2k[jj])*coskz_2[kk];
                    double depsilonz_dkz      = tz*cosk_2[ii]*cosk_2[jj]*(cosk[ii]-cosk[jj])*(cosk[ii]-cosk[jj])*sinkz_2[kk];
                    
                    double ep_k        =  epsilon_k   + epsilonz_k;
                    double dep_dkx     = depsilon_dkx + depsilonz_dkx;
                    double dep_dky     = depsilon_dky + depsilonz_dky;
                    double dep_dkx_dky = depsilon_dkx_dky + depsilonz_dkx_dky;
                    double dep_dky_dky = depsilon_dky_dky + depsilonz_dky_dky;
                    double dep_dkz     = depsilonz_dkz;
                    
                    double Gamma = ETA;
                    double normV_k = sqrt(dep_dkx*dep_dkx + dep_dky*dep_dky + dep_dkz*dep_dkz);
                    Gamma += ETAell*normV_k;
                    Gamma += ETAdos/(normV_k+10e-10);
                    Gamma += ETAk*cos_phi12;
                    double Ak_0   = -(1./M_PI)*cimag(1.0/ (I*Gamma - ep_k));

                    //// T-INDEPENDENT
                    double kernel_xx  = dep_dkx*dep_dkx;
                    double kernel_zz  = dep_dkz*dep_dkz;
                    double kernel_xy  = -(2./3.)*(dep_dkx * (dep_dkx * dep_dky_dky - dep_dky * dep_dkx_dky));
                    
                    dos[iMu]       += Ak_0;
                    density_0[iMu] += 1.0/(1.0+exp(1000*ep_k));
                    sigmaxx_0[iMu] += kernel_xx * Ak_0*Ak_0;
                    sigmaxy_0[iMu] += kernel_xy * Ak_0*Ak_0*Ak_0;
                    sigmazz_0[iMu] += kernel_zz * Ak_0*Ak_0;
                    
                    //// LOOP ON TEMPERATURES
                    int nn=0; for(nn=0; nn<nT; nn++)
                    {
                        density[iMu][nn]    += 1.0/(1.0+exp(beta[nn]*ep_k));
                        double GammaT = Gamma + ETAaPL*T[nn];
                               GammaT += M_PI*M_PI*ETAaFL*T[nn]*T[nn];
                               GammaT += M_PI*M_PI*ETAbFL*T[nn]*T[nn]/(normV_k+10e-10);

                        //// INTEGRAL IN ENERGY
                        int n=0; for(n=0; n<nOmega; n++)
                        {
                            double GammaOmega = GammaT;
                                   GammaOmega += ETAaFL*omega[nn][n]*omega[nn][n];
                                   GammaOmega += ETAbFL*omega[nn][n]*omega[nn][n]/(normV_k+10e-10);
                                   GammaOmega += ETAw*omega[nn][n];
                                   GammaOmega += ETAkw*omega[nn][n]*cos_phi12;
                            
                            double complex z = omega[nn][n] + GammaOmega * I;
                            double A_k = -(1./M_PI)*cimag(1.0/ (z-ep_k) );

                            double CvKernel = omega[nn][n]*dfermiDirac_dT[nn][n]*A_k;
                            Cv[iMu][nn]       += CvKernel;

                            double frequencyKernel_xx = -dfermiDirac_dw[nn][n]*kernel_xx*A_k*A_k;
                            double frequencyKernel_zz = -dfermiDirac_dw[nn][n]*kernel_zz*A_k*A_k;
                            double frequencyKernel_xy = -dfermiDirac_dw[nn][n]*kernel_xy*A_k*A_k*A_k;
                            
                            sigma_xx[iMu][nn] += frequencyKernel_xx;
                            sigma_zz[iMu][nn] += frequencyKernel_zz;
                            sigma_xy[iMu][nn] += frequencyKernel_xy;
                            double betaOmega = beta[nn]*omega[nn][n];
                            alpha_xx[iMu][nn] += betaOmega * frequencyKernel_xx;
                            alpha_zz[iMu][nn] += betaOmega * frequencyKernel_zz;
                            alpha_xy[iMu][nn] += betaOmega * frequencyKernel_xy;
                            double betaOmega2 = beta[nn]*beta[nn] * omega[nn][n]*omega[nn][n];
                            beta_xx[iMu][nn]  += betaOmega2 * frequencyKernel_xx;
                            beta_zz[iMu][nn]  += betaOmega2 * frequencyKernel_zz;
                            beta_xy[iMu][nn]  += betaOmega2 * frequencyKernel_xy;
                        }
                    }
                }
            }
        }
        
        double f_0 = 2.0/nK/nK/nKz;
        printf("% 4.8f % 4.8f % 4.8f ", mu, 1.0-f_0*density_0[iMu], f_0*dos[iMu]);
        printf("% 4.8f % 4.8f % 4.8f ", f_0*sigmaxx_0[iMu], f_0*sigmaxy_0[iMu], f_0*sigmazz_0[iMu]);
        printf("\n");
        fflush(stdout);

        
        //// FILE GROUPED BY T
        fprintf(fileOut, "mu\tn0\tp0\tdos\teta\t");
        fprintf(fileOut, "sigmaxx0\tsigmaxy0\tsigmazz0\t");
        fprintf(fileOut, "T\tn\tp\tCv\t");
        fprintf(fileOut, "sigmaxx\tsigmaxy\tsigmazz\t");
        fprintf(fileOut, "alphaxx\talphaxy\talphazz\t");
        fprintf(fileOut, "betaxx\tbetaxy\tbetazz\t");
        fprintf(fileOut, "\n");

        //// LOOP ON TEMPERATURES
        nn=0; for(nn=0; nn<nT; nn++)
        {
            double f = (2.*energyCutoff[nn]) /(nOmega)*f_0;
            fprintf(fileOut,"%e\t%e\t%e\t%e\t%e\t", mu, f_0*density_0[iMu], 1.0-f_0*density_0[iMu], f_0*dos[iMu], etaNorm);
            fprintf(fileOut,"%e\t%e\t%e\t", f_0*sigmaxx_0[iMu], f_0*sigmaxy_0[iMu], f_0*sigmazz_0[iMu]);
            fprintf(fileOut,"%e\t%e\t%e\t%e\t", T[nn],f_0*density[iMu][nn], 1.0-f_0*density[iMu][nn], f*Cv[iMu][nn]);
            fprintf(fileOut,"%e\t%e\t%e\t", f*sigma_xx[iMu][nn], f*sigma_xy[iMu][nn], f*sigma_zz[iMu][nn] );
            fprintf(fileOut,"%e\t%e\t%e\t", f*alpha_xx[iMu][nn], f*alpha_xy[iMu][nn], f*alpha_zz[iMu][nn] );
            fprintf(fileOut,"%e\t%e\t%e\t", f*beta_xx[iMu][nn],  f*beta_xy[iMu][nn],  f*beta_zz[iMu][nn] );

            fprintf(fileOut, "\n");
        }
        fprintf(fileOut, "\n\n");
        fflush(fileOut);
        
        
    }
    fclose(fileOut);
    
    
    //// FILE GROUPED BY mu
    char nameFileMu[800] = "transport_vs_Mu";
    strcat(nameFileMu, model);
    FILE *fileOutMu = fopen(nameFileMu,"w");

    nn=0; for(nn=0; nn<nT; nn++)
    {
        
        double f_0 = 2.0/nK/nK/nKz;
        double f = (2.*energyCutoff[nn])/(nOmega)*f_0;
        
        fprintf(fileOutMu, "mu\tn0\tp0\tdos\teta\t");
        fprintf(fileOutMu, "sigmaxx0\tsigmaxy0\tsigmazz0\t");
        fprintf(fileOutMu, "T\tn\tp\tCv\t");
        fprintf(fileOutMu, "sigmaxx\tsigmaxy\tsigmazz\t");
        fprintf(fileOutMu, "alphaxx\talphaxy\talphazz\t");
        fprintf(fileOutMu, "betaxx\tbetaxy\tbetazz\t");
        fprintf(fileOutMu, "\n");

        for(int iMu=0; iMu<nMu; iMu++)
        {
            double mu = muMin + iMu*(muMax-muMin)/(nMu-1);
        
            fprintf(fileOutMu,"%e\t%e\t%e\t%e\t%e\t", mu, f_0*density_0[iMu], 1.0-f_0*density_0[iMu], f_0*dos[iMu], etaNorm);
            fprintf(fileOutMu,"%e\t%e\t%e\t", f_0*sigmaxx_0[iMu], f_0*sigmaxy_0[iMu], f_0*sigmazz_0[iMu]);
            fprintf(fileOutMu,"%e\t%e\t%e\t%e\t", T[nn],f_0*density[iMu][nn], 1.0-f_0*density[iMu][nn], f*Cv[iMu][nn]);
            fprintf(fileOutMu,"%e\t%e\t%e\t", f*sigma_xx[iMu][nn], f*sigma_xy[iMu][nn], f*sigma_zz[iMu][nn] );
            fprintf(fileOutMu,"%e\t%e\t%e\t", f*alpha_xx[iMu][nn], f*alpha_xy[iMu][nn], f*alpha_zz[iMu][nn] );
            fprintf(fileOutMu,"%e\t%e\t%e\t", f*beta_xx[iMu][nn],  f*beta_xy[iMu][nn],  f*beta_zz[iMu][nn] );
            fprintf(fileOutMu, "\n");
        
        }
        fprintf(fileOutMu, "\n\n");
    }
    
    printf("\n transportBuddy over.\n");
    return 0;
}

