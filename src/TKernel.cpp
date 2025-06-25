#include <iostream>
#include <math.h>
#include <sstream>
#include <gsl/gsl_integration.h>
#include "TKernel.h"
#include "TFastReso_formulas.h"

using namespace std;
// struct for eta and phi numerical integrals
struct dK_dphideta_params{
  double phi,pT,ur;
  int index;
  int m;
  TKernel * kernel_class;
};
struct dK_dphi_params{
  gsl_function  *Fdphideta;
  gsl_integration_workspace *weta ;
  double etamax;
};
double dK_dphideta(double eta, void * p) {
  double &phi          = ((struct dK_dphideta_params *) p)->phi;
  double &pT           = ((struct dK_dphideta_params *) p)->pT;
  double &ur           = ((struct dK_dphideta_params *) p)->ur;
  int &index           = ((struct dK_dphideta_params *) p)->index;
  int &m           = ((struct dK_dphideta_params *) p)->m;
  TKernel * kernel_class  = ((struct dK_dphideta_params *) p)->kernel_class;
  //double dK_mu_detadphi[8];
  //kernel_class->get_K_mu_detadphi(eta, phi, pT, ur, dK_mu_detadphi);
  double dK_mu_detadphi = kernel_class->get_K_mu_detadphi_i(eta, phi, pT, ur,  index, m);
  return 2*dK_mu_detadphi;
}

double dK_dphi(double phi, void * p) {
  gsl_function *Fdphideta          = ((struct dK_dphi_params *) p)->Fdphideta;
  gsl_integration_workspace *weta  = ((struct dK_dphi_params *) p)->weta;
  double &etamax                   = ((struct dK_dphi_params *) p)->etamax;
  ((struct dK_dphideta_params *) Fdphideta->params)->phi = phi;
  double result, error;
  using namespace qag_params;
int status =  gsl_integration_qag (Fdphideta, 0.0, etamax,  fEpsAbs,fEpsRel,fLimit,fKey,weta, &result, &error);
  // Check status of integration and reduce accuracy if failing
  if (status) {
status =  gsl_integration_qag (Fdphideta, 0.0, etamax,  fEpsAbs,1e-4,fLimit,GSL_INTEG_GAUSS15,weta, &result, &error);
    if (status) {
status =  gsl_integration_qag (Fdphideta, 0.0, etamax,  1e-4,1e-4,fLimit,GSL_INTEG_GAUSS15,weta, &result, &error);
      if (status) {      std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << result << "+_" << error << std::endl;
        exit(EXIT_FAILURE);
      }else {
        std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << result << "+-" << error <<std::endl;
    }
  }


  return 2*result;
}

TKernel::TKernel(TParticle *particle): fParticle(particle)  {
  initialize();
}
TKernel::~TKernel(){
  gsl_integration_workspace_free (fWeta);
  gsl_integration_workspace_free (fWphi);
}
void TKernel::initialize(){
  fWphi = gsl_integration_workspace_alloc (qag_params::fLimit);
  fWeta = gsl_integration_workspace_alloc (qag_params::fLimit);
  }

// pT = pT ur =sinh(chi)
double TKernel::get_K_mu_detadphi_i(double eta, double phi,double pT, double ur, int index, int m) {

  double fMb = fParticle->getM();
  int fNu = fParticle->getNu();
  double Ebar=sqrt(fMb*fMb+pT*pT)*sqrt(1+ur*ur)*cosh(eta)-pT*ur*cos(phi);

  double pbar = sqrt(Ebar*Ebar-fMb*fMb);
  switch (index) {
    case 0:
      {
        double feq1 = fNu*fParticle->get_Fj(0, Ebar)/pbar;
        double feq2 = fNu*fParticle->get_Fj(1, Ebar)/pbar;
        return feq1*cosh(eta)*sqrt(fMb*fMb+pT*pT) + (feq2-feq1)*Ebar*sqrt(1+ur*ur); //component 0 (tau)
      }
      break;
    case 1:
      {
        double feq1 = fNu*fParticle->get_Fj(0, Ebar)/pbar;
        double feq2 = fNu*fParticle->get_Fj(1, Ebar)/pbar;
        return feq1*cos(phi)*pT + (feq2-feq1)*Ebar*ur; //component 1 (r)
      }
      break;
    case 2: //mu = 0, taking eta component
      { 
        double fshear1 = fNu*fParticle->get_Fj(2, Ebar)/pbar/pbar/pbar;
        double fshear2 = fNu*fParticle->get_Fj(3, Ebar)/pbar/pbar/pbar;
        double fshear3 = fNu*fParticle->get_Fj(4, Ebar)/pbar/pbar/pbar;
        return (fshear1*cosh(eta)*sqrt(fMb*fMb+pT*pT) + (fshear3-fshear1)*Ebar*sqrt(1+ur*ur))*
          ((fMb*fMb+pT*pT)*sinh(eta)*sinh(eta)-pow(sqrt(fMb*fMb+pT*pT)*ur*cosh(eta)-pT*sqrt(1+ur*ur)*cos(phi), 2)  )
          +(fshear2-fshear1)*2.0/5.0*(Ebar*Ebar-fMb*fMb)*ur*(sqrt(fMb*fMb+pT*pT)*ur*cosh(eta)-pT*sqrt(1+ur*ur)*cos(phi)) ;
      }
      break;
    case 3: 
      {
        double fshear1 = fNu*fParticle->get_Fj(2, Ebar)/pbar/pbar/pbar;
        double fshear2 = fNu*fParticle->get_Fj(3, Ebar)/pbar/pbar/pbar;
        double fshear3 = fNu*fParticle->get_Fj(4, Ebar)/pbar/pbar/pbar;
        return (fshear1*cos(phi)*pT + (fshear3-fshear1)*Ebar*ur)*
          ((fMb*fMb+pT*pT)*sinh(eta)*sinh(eta)-pow(sqrt(fMb*fMb+pT*pT)*ur*cosh(eta)-pT*sqrt(1+ur*ur)*cos(phi), 2)  )
          +(fshear2-fshear1)*2.0/5.0*(Ebar*Ebar-fMb*fMb)*sqrt(1+ur*ur)*(sqrt(fMb*fMb+pT*pT)*ur*cosh(eta)-pT*sqrt(1+ur*ur)*cos(phi));
      }
      break;
    case 4:
      {
        double fshear1 = fNu*fParticle->get_Fj(2, Ebar)/pbar/pbar/pbar;
        double fshear2 = fNu*fParticle->get_Fj(3, Ebar)/pbar/pbar/pbar;
        double fshear3 = fNu*fParticle->get_Fj(4, Ebar)/pbar/pbar/pbar;
        return (fshear1*cosh(eta)*sqrt(fMb*fMb+pT*pT)+ (fshear3-fshear1)*Ebar*sqrt(1+ur*ur))*
          (pT*pT*sin(phi)*sin(phi)-pow(sqrt(fMb*fMb+pT*pT)*ur*cosh(eta)-pT*sqrt(1+ur*ur)*cos(phi), 2)  )
          +(fshear2-fshear1)*2.0/5.0*(Ebar*Ebar-fMb*fMb)*ur*(sqrt(fMb*fMb+pT*pT)*ur*cosh(eta)-pT*sqrt(1+ur*ur)*cos(phi));
      }
      break;
    case 5:
      {
        double fshear1 = fNu*fParticle->get_Fj(2, Ebar)/pbar/pbar/pbar;
        double fshear2 = fNu*fParticle->get_Fj(3, Ebar)/pbar/pbar/pbar;
        double fshear3 = fNu*fParticle->get_Fj(4, Ebar)/pbar/pbar/pbar;
        return (fshear1*cos(phi)*pT + (fshear3-fshear1)*Ebar*ur)*
          (pT*pT*sin(phi)*sin(phi)-pow(sqrt(fMb*fMb+pT*pT)*ur*cosh(eta)-pT*sqrt(1+ur*ur)*cos(phi), 2)  )
          +(fshear2-fshear1)*2.0/5.0*(Ebar*Ebar-fMb*fMb)*sqrt(1+ur*ur)*(sqrt(fMb*fMb+pT*pT)*ur*cosh(eta)-pT*sqrt(1+ur*ur)*cos(phi));
      }
      break;
    case 6:
      {
        double fbulk1 = fNu*fParticle->get_Fj(5, Ebar)/pbar;
        double fbulk2 = fNu*fParticle->get_Fj(6, Ebar)/pbar;
        return fbulk1*cosh(eta)*sqrt(fMb*fMb+pT*pT) + (fbulk2-fbulk1)*Ebar*sqrt(1+ur*ur); //component 0 (tau) 
      }
      break;
    case 7:
      {
        double fbulk1 = fNu*fParticle->get_Fj(5, Ebar)/pbar;
        double fbulk2 = fNu*fParticle->get_Fj(6, Ebar)/pbar;
        return fbulk1*cos(phi)*pT + (fbulk2-fbulk1)*Ebar*ur;
      }
      break;
    case 8:
      {
        double ftemp1 = fNu*fParticle->get_Fj(7, Ebar)/pbar;
        double ftemp2 = fNu*fParticle->get_Fj(8, Ebar)/pbar;
        return (ftemp1*cosh(eta)*sqrt(fMb*fMb+pT*pT) + (ftemp2-ftemp1)*Ebar*sqrt(1+ur*ur) )*cos(m*phi);
      }
      break;
    case 9:
      {
        double ftemp1 = fNu*fParticle->get_Fj(7, Ebar)/pbar;
        double ftemp2 = fNu*fParticle->get_Fj(8, Ebar)/pbar;
        return (ftemp1*cos(phi)*pT + (ftemp2-ftemp1)*Ebar*ur)*cos(m*phi);
      }
      break;
    case 10:
      {
        double fvel1 = fNu*fParticle->get_Fj(9, Ebar)/pbar/pbar;
        double fvel2 = fNu*fParticle->get_Fj(10, Ebar)/pbar/pbar;
        double fvel3 = fNu*fParticle->get_Fj(11, Ebar)/pbar/pbar;
        return ((fvel1*cosh(eta)*sqrt(fMb*fMb+pT*pT) + (fvel3-fvel1)*Ebar*sqrt(1+ur*ur))*
          (-(fMb*fMb+pT*pT)*ur*cosh(eta)+ pT*sqrt(1+ur*ur)*cos(phi)  )
          +(fvel2-fvel1)*1.0/3.0*(Ebar*Ebar-fMb*fMb)*ur )*cos(m*phi) ;
      }
      break;
    case 11:
      {
        double fvel1 = fNu*fParticle->get_Fj(9, Ebar)/pbar/pbar;
        double fvel2 = fNu*fParticle->get_Fj(10, Ebar)/pbar/pbar;
        double fvel3 = fNu*fParticle->get_Fj(11, Ebar)/pbar/pbar;
        return ((fvel1*cos(phi)*pT + (fvel3-fvel1)*Ebar*ur)*
          (-(fMb*fMb+pT*pT)*ur*cosh(eta)+ pT*sqrt(1+ur*ur)*cos(phi)  )
          +(fvel2-fvel1)*1.0/3.0*(Ebar*Ebar-fMb*fMb)*sqrt(1+ur*ur) )*cos(m*phi) ;
      }
      break;
    case 12:
      {
        double fvel1 = fNu*fParticle->get_Fj(9, Ebar)/pbar/pbar;
        //double fvel2 = fNu*fParticle->get_Fj(10, Ebar)/pbar/pbar;
        double fvel3 = fNu*fParticle->get_Fj(11, Ebar)/pbar/pbar;
        return ((fvel1*cosh(eta)*sqrt(fMb*fMb+pT*pT) + (fvel3-fvel1)*Ebar*sqrt(1+ur*ur))* (pT*sin(phi)) )*sin(m*phi) ; //phi component
      }
      break;
    case 13:
      {
        double fvel1 = fNu*fParticle->get_Fj(9, Ebar)/pbar/pbar;
        //double fvel2 = fNu*fParticle->get_Fj(10, Ebar)/pbar/pbar;
        double fvel3 = fNu*fParticle->get_Fj(11, Ebar)/pbar/pbar;
        return ((fvel1*cos(phi)*pT + (fvel3-fvel1)*Ebar*ur)* (pT*sin(phi)) )*sin(m*phi) ; //eta component 
      }
      break;
    case 14:
      {
        double Lambda=3.14*0.155;
        double v0 = 0.5;
        double m0 = 0.1;
        double y = pbar/Lambda;
        double fac = 1/(1+y*y/2+pow(y,4));
        fMb = sqrt(fMb*fMb*(1-fac)+m0*m0*fac);
        double v = sqrt(1-fac+v0*v0*fac);
        Ebar=sqrt(fMb*fMb+v*v*pbar*pbar);
        double QMu=0.0;
        double Tfo = 0.155;
        double feq1 = fNu*get_thermal_F(Ebar, Tfo, QMu, fParticle->getType()); //fNu*fParticle->get_Fj(0, Ebar)/pbar;
        return feq1*cosh(eta)*sqrt(fMb*fMb+pT*pT) ;// + (feq2-feq1)*Ebar*sqrt(1+ur*ur);
      }
      break;
    case 15:
      {
        double Lambda=3.14*0.155;
        double v0 = 0.5;
        double m0 = 0.1;
        double y = pbar/Lambda;
        double fac = 1/(1+y*y/2+pow(y,4));
        fMb = sqrt(fMb*fMb*(1-fac)+m0*m0*fac);
        double v = sqrt(1-fac+v0*v0*fac);
        Ebar=sqrt(fMb*fMb+v*v*pbar*pbar);
        double QMu=0.0;
        double Tfo=0.155;
        double feq1 = fNu*get_thermal_F(Ebar, Tfo, QMu, fParticle->getType()); //fNu*fParticle->get_Fj(0, Ebar)/pbar;
        return feq1*cos(phi)*pT ; //+ (feq2-feq1)*Ebar*ur;
      }
      break;
    case 16:
      {
        double fdiff1 = fNu*fParticle->get_Fj(12, Ebar)/pbar/pbar; //divide and multiply by pbar (in FastReso_formulas.h) for numerical stability
        double fdiff2 = fNu*fParticle->get_Fj(13, Ebar)/pbar/pbar;
        double fdiff3 = fNu*fParticle->get_Fj(14, Ebar)/pbar/pbar;
        //return -ur*sqrt(fMb*fMb+pT*pT)*cosh(eta)*(Ebar*fdiff3+2*ur*pT*cos(phi)*fdiff1)+ur*ur*ur/sqrt(1+ur*ur)*(fMb*fMb+pT*pT)*cosh(eta)*cosh(eta)*fdiff1+
        //1.0/3.0*(Ebar*Ebar-fMb*fMb)*ur/sqrt(1+ur*ur)*(fdiff2-fdiff1)
        //+sqrt(1+ur*ur)*pT*Ebar*cos(phi)*fdiff3+sqrt(1+ur*ur)*ur*pT*pT*cos(phi)*cos(phi)*fdiff1;

        return ((fdiff1*cosh(eta)*sqrt(fMb*fMb+pT*pT) + (fdiff3-fdiff1)*Ebar*sqrt(1+ur*ur))*
          (-(fMb*fMb+pT*pT)*ur*cosh(eta)+ pT*sqrt(1+ur*ur)*cos(phi)  )
          +(fdiff2-fdiff1)*1.0/3.0*(Ebar*Ebar-fMb*fMb)*ur ) ;
      }
      break;
    case 17:
      {
        double fdiff1 = fNu*fParticle->get_Fj(12, Ebar)/pbar/pbar;
        double fdiff2 = fNu*fParticle->get_Fj(13, Ebar)/pbar/pbar;
        double fdiff3 = fNu*fParticle->get_Fj(14, Ebar)/pbar/pbar;
        
        //return -ur*sqrt(fMb*fMb+pT*pT)*cosh(eta)*(2*sqrt(1+ur*ur)*pT*cos(phi)*fdiff1+ur/sqrt(1+ur*ur)*Ebar*fdiff3)+ur*ur*cosh(eta)*cosh(eta)*(fMb*fMb+pT*pT)*fdiff1+
        //(1+ur*ur)*pT*pT*cos(phi)*cos(phi)*fdiff1
        //+ur*pT*Ebar*cos(phi)*fdiff3+1.0/3.0*(Ebar*Ebar-fMb*fMb)*(fdiff2-fdiff1);


        return ((fdiff1*cos(phi)*pT + (fdiff3-fdiff1)*Ebar*ur)*
          (-(fMb*fMb+pT*pT)*ur*cosh(eta)+ pT*sqrt(1+ur*ur)*cos(phi)  )
          +(fdiff2-fdiff1)*1.0/3.0*(Ebar*Ebar-fMb*fMb)*sqrt(1+ur*ur) ) ;
      }
      break;


    default:
      cerr  << "\033[1mTKernel.cpp::get_K_mu_detadphi_i\033[0m : \033[1;31merror\033[0m :  index out of bounds i =  "  << index << endl;
      exit(EXIT_FAILURE);

  }

}

double TKernel::get_K_mu_qag(double pT, double ur, int index, int m) {
  using namespace qag_params;
  double K_mu =0;
  gsl_function Fdphi;
  Fdphi.function = dK_dphi;
  gsl_function Fdphideta;
  Fdphideta.function = dK_dphideta;
  dK_dphideta_params dphideta_params = {0.0, pT,ur,index, m, this};
  Fdphideta.params = &dphideta_params;
  dK_dphi_params dphi_params = {&Fdphideta,fWeta,fEtaMax};
  Fdphi.params = &dphi_params;
  double error;
  int status;

  status = gsl_integration_qag (&Fdphi, 0.0, M_PI,fEpsAbs,fEpsRel,fLimit,fKey ,fWphi, &K_mu, &error);

  if (status) {
  status = gsl_integration_qag (&Fdphi, 0.0, M_PI,fEpsAbs,1e-4,fLimit,GSL_INTEG_GAUSS15 ,fWphi, &K_mu, &error);
    if (status) {
  status = gsl_integration_qag (&Fdphi, 0.0, M_PI,1e-4,1e-4,fLimit,GSL_INTEG_GAUSS15 ,fWphi, &K_mu, &error);
      if (status) {      std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31merror\033[0m : integration failure! status = " << status << ". " << K_mu << "+_" << error << std::endl;
        exit(EXIT_FAILURE);
      }else {
        std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : reduced relative and absolute error (1e-4) integration. Result = " << K_mu<< "+-" << error <<std::endl;
      }
    } else {
      std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : reduced relative error (1e-4) integration. Result = " << K_mu << "+-" << error <<std::endl;
    }

  }
return K_mu;
}

void TKernel::print(string tag, string columns){

    stringstream sstr(columns);
    string col;
    bool IsKeq=false;
    bool IsKshear=false;
    bool IsKbulk=false;
    bool IsKtemp=false;
    bool IsKvel=false;
    bool IsKcrit=false;
    bool IsKdiff=false;
    int m = 0;
    while (getline(sstr, col, ' ')){
      if (col=="Keq")         { IsKeq=true; }
      else if (col=="Kshear") { IsKshear=true; }
      else if (col=="Kbulk")  { IsKbulk=true; }
      else if (col=="Ktemp")  { IsKtemp=true; }
      else if (col=="Kvel")   { IsKvel=true; }
      else if (col=="Kcrit")   { IsKcrit=true; }
      else if (col=="Kdiff")   { IsKdiff=true; }
      else if (col[0]=='m')   { m=atoi((col.substr(2)).c_str()) ;
      }
      else { std::cerr << "\033[1mTKernel.cpp\033[0m : \033[1;31mwarning\033[0m : unkwnown label " << col  <<std::endl; }

    }

    FILE * pFile;
    string fname = tag+"_Kj.out";
    pFile = fopen (fname.c_str(),"w");
    fprintf(pFile,"#%15s\t%s", "1:pT [GeV]","2:Ur");
    int index = 3;
    if (IsKeq) {
      fprintf(pFile,"\t%d:%s\t%d:%s", index, "Keq 1", index+1, "Keq 2");
      index +=2;
    }
   if (IsKshear) {
      fprintf(pFile,"\t%d:%s\t%d:%s\t%d:%s\t%d:%s",
          index, "Kshear 1", index+1, "Kshear 2" , index+2,  "Kshear 3", index+3, "Kshear 4");
      index +=4;
    }
   if (IsKbulk) {
      fprintf(pFile,"\t%d:%s\t%d:%s", index,  "Kbulk 1", index+1, "Kbulk 2");
      index +=2;
    }
    if (IsKtemp) {
      for (int im=0; im<m; im++){
        fprintf(pFile,"\t%d:%s,%d\t%d:%s,%d" , index, "Ktemp 1",im, index+1, "Ktemp 2", im);
        index +=4;
      }
    }
    if (IsKcrit) {
      fprintf(pFile,"\t%d:%s\t%d:%s", index, "Kcrit 1", index+1, "Kcrit 2");
        index +=2;
    }
    if (IsKdiff) {
      fprintf(pFile,"\t%d:%s\t%d:%s", index, "Kdiff 1", index+1, "Kdiff 2");
      index +=2;
    }

  fprintf(pFile,"\n");
  for(int j = 0; j <grid_params::fNur; j++){
      double ur = tan(atan(grid_params::fUrMax)*(j)/(grid_params::fNur-1));
    for(int i = 0; i <grid_params::fNpT; i++){
    double pT = grid_params::fMref*tan(atan(grid_params::fPTMax/grid_params::fMref)*(i+0.5)/(grid_params::fNpT-1));
      fprintf(pFile,"%15.5e\t%15.5e", pT, ur);

        if (IsKeq) {
      double Keq1 = get_K_mu_qag(pT,  ur, 0) ;
      double Keq2 = get_K_mu_qag(pT,  ur, 1) ;
        fprintf(pFile,"\t%15.5e\t%15.5e", Keq1, Keq2);
        }
       if (IsKshear) {
      double Kshear1 = get_K_mu_qag(pT,  ur, 2) ;
      double Kshear2 = get_K_mu_qag(pT,  ur, 3) ;
      double Kshear3 = get_K_mu_qag(pT,  ur, 4) ;
      double Kshear4 = get_K_mu_qag(pT,  ur, 5) ;
        fprintf(pFile,"\t%15.5e\t%15.5e\t%15.5e\t%15.5e", Kshear1,Kshear2,Kshear3,Kshear4);
        }
      if (IsKbulk) {
      double Kbulk1 = get_K_mu_qag(pT,  ur, 6) ;
      double Kbulk2 = get_K_mu_qag(pT,  ur, 7) ;
        fprintf(pFile,"\t%15.5e\t%15.5e", Kbulk1, Kbulk2);
        }
       if (IsKtemp) {
          for (int im=0; im<m; im++){
      double Ktemp1 = get_K_mu_qag(pT,  ur, 8, im) ;
      double Ktemp2 = get_K_mu_qag(pT,  ur, 9, im) ;
        fprintf(pFile,"\t%15.5e\t%15.5e", Ktemp1, Ktemp2);
          }
        }
      if (IsKvel) {
          for (int im=0; im<m; im++){
      double Kvel1 = get_K_mu_qag(pT,  ur, 10, im) ;
      double Kvel2 = get_K_mu_qag(pT,  ur, 11, im) ;
      double Kvel3 = get_K_mu_qag(pT,  ur, 12, im) ;
      double Kvel4 = get_K_mu_qag(pT,  ur, 13, im) ;
        fprintf(pFile,"\t%15.5e\t%15.5e\t%15.5e\t%15.5e", Kvel1, Kvel2, Kvel3, Kvel4);
          }
        }
      if (IsKcrit) {
      double Kcrit1 = get_K_mu_qag(pT,  ur, 14) ;
      double Kcrit2 = get_K_mu_qag(pT,  ur, 15) ;
        fprintf(pFile,"\t%15.5e\t%15.5e", Kcrit1, Kcrit2);
        }
      if (IsKdiff) {
      double Kdiff1 = get_K_mu_qag(pT,  ur, 16) ;
      double Kdiff2 = get_K_mu_qag(pT,  ur, 17) ;
        fprintf(pFile,"\t%15.5e\t%15.5e", Kdiff1, Kdiff2);
        }

      fprintf(pFile,"\n");
    }
    fprintf(pFile,"\n");
  }
  fclose(pFile);
}
