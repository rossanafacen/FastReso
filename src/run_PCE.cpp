/*
 * Copyright (c) 2018-2021 Aleksas Mazeliauskas, Stefan Floerchinger, 
 *                    Eduardo Grossi, and Derek Teaney
 * All rights reserved.
 *
 * FastReso is distributed under MIT license;
 * see the LICENSE file that should be present in the root
 * of the source distribution, or alternately available at:
 * https://github.com/amazeliauskas/FastReso/
 */

#include <iostream>
#include <cmath>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <TFastReso_AZYHYDRO.h>
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_multiroots.h"
#include "qag_params.h"
#include "eos_func.h"
#include <format>
#include <omp.h>

#include <nlohmann/json.hpp>
#include <fstream>

using namespace std ;
using json = nlohmann::json;

struct rparams
{
  double **partial_charges;
  double *masses;
  int *sign;
  int *deg;
  double *ci;
  double Tkin;
  int NPart;
  int NCharge;
};
  void
print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = % .3f % .3f "
      "f(x) = % .3e % .3e\n",
      (int) iter,
      gsl_vector_get (s->x, 0),
      gsl_vector_get (s->x, 10),
      gsl_vector_get (s->f, 0),
      gsl_vector_get (s->f, 10));
}

int partial_f (const gsl_vector * x, void *params,
    gsl_vector * f)
{
  double **partial_charges = ((struct rparams *) params)->partial_charges;
  double *masses = ((struct rparams *) params)->masses;
  int *sign = ((struct rparams *) params)->sign;
  int *deg = ((struct rparams *) params)->deg;
  double *ci = ((struct rparams *) params)->ci;
  double &Tkin = ((struct rparams *) params)->Tkin;
  int &Npart = ((struct rparams *) params)->NPart;
  int &Ncharge = ((struct rparams *) params)->NCharge;

  double * fi ;
  fi = new double[Ncharge];

  for (int i=0; i<Ncharge; i++){ fi[i]=0; }

  double s=0;
  for (int i=0; i<Npart; i++){
    if (masses[i] < 0.01) continue;
    double QMu=0; for (int pj=0; pj<Ncharge; pj++){ QMu +=partial_charges[i][pj]*gsl_vector_get (x, pj); }
    s += deg[i]*sfunc( masses[i]/Tkin,QMu/Tkin, sign[i]);
    for (int pi=0; pi<Ncharge; pi++){ fi[pi] += partial_charges[i][pi]*deg[i]*nfunc( masses[i]/Tkin,QMu/Tkin, sign[i]); }
  }

  for (int pj=0; pj<Ncharge; pj++){ gsl_vector_set (f, pj, fi[pj]/s-ci[pj]); }
  delete [] fi;
  return GSL_SUCCESS;
}
int main(int argc, const char * argv[])
{
  bool verbose = false;
  //bool four_decays_verbose = false;

  if (argc<6) {
    cerr << argv[0] << " not enough input arguments " << endl;
    cerr << argv[0] << " reverse_decays_PDG2016_massordered.dat  decays_PDG2016_massordered.dat ./outputfolder/ /json_file/ job_id" << endl;
    exit(EXIT_FAILURE) ;
  }
/////////////////// DANGER ///////////////////////////////////////////////////
  // If integration fails, you can turn the handler off. Integration functions
  // in TFastReso_AZYHYDRO.cpp and TFastReso_formulas.h will try reducing accuracy
  // for the failing cases.  Use with caution!
  gsl_set_error_handler_off ();
/////////////////// DANGER ///////////////////////////////////////////////////

  // list of all particles
  string pdata(argv[1]); // = "reverse_decays_PDG2016_massordered.dat";
  string inputname(argv[1]);
  string ddata(argv[2]);// = "decays_without_weak.data";
  string tag(argv[3]);
  string config(argv[4]); 
  int jobposition = std::stoi(argv[5]);

  
  ifstream file(config);
  json j;
  file >> j;
  
  int ns = j["physics"]["pce"]["n_metastable"];
  int pdg[ns];
  for(int i = 0; i < ns; i++) {
        pdg[i] = j["physics"]["pce"]["metastable_particles"][i];
  } 
 
  
  int ns2 = j["physics"]["n_particles"];
  int outpdg[ns2];
  for(int i = 0; i < ns2; i++) {
        outpdg[i] = j["physics"]["particles"][i];
  } 

  
  
  double Tchem = j["physics"]["Tchem"][jobposition];
  const double Tmin = j["physics"]["Tkin"][0];
  const double Tmax = j["physics"]["Tkin"][1];
  double dT = j["physics"]["dT"];
  int NT = round((Tmax - Tmin) / dT) + 1;

  bool dothermal = j["physics"]["thermal"];
  
  clog << "T chem " << Tchem << " GeV" << endl; 


  TFastReso_AZYHYDRO fastreso;

  // read all particles
  fastreso.read_particles_data(pdata);
  int fParticleNumber = fastreso.fParticleNumber;
  double **partial_charges = new double *[fParticleNumber];

  for (int i=0; i<fParticleNumber; i++){
    partial_charges[i] = new double [ns];
    for (int pi=0; pi<ns; pi++){
      partial_charges[i][pi] =0.0;
    }
  }

  for (int pi=0; pi<ns; pi++){
    partial_charges[fastreso.getParticleIndexByPDG(pdg[pi])][pi] +=1;
  }
  string line;

  // Open the input file
  ifstream fInputFile(inputname);
  if (fInputFile.fail()) {
    cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m :  input file does not exist "  << inputname<< endl;
    exit(EXIT_FAILURE);
  }
  if (verbose){
    // Print the header
    printf("#%9s\t%9s\t%9s\t%9s\t%9s\n",
        "Father", "Child1", "Child2","Child3","BranchingRatio");
  }

  int Father, Child1, Child2, Child3, Child4, Child5;
  int ndecayproducts;
  string Name;
  double BranchingRatio;

  int count_decays=0;
  int count_2decays=0;
  int count_3decays=0;

  while (getline(fInputFile, line)){
    bool FathernotQuasiStable=true;
    //! check that there is the correct number of elements per line (5 for 2-body decay, 6 for three body decay), continue otherwise.
    stringstream sstrc(line); int count=0; while (sstrc >> Name){count++;} if (not (count==8)) continue;
    stringstream sstr(line);
    sstr >> Name; 
    // check that the first non-whitespace character is not #, i.e. commented line
    if(Name.front()!='#'){

      Father = atoi(Name.c_str());
      sstr>> ndecayproducts >> BranchingRatio >>  Child1 >> Child2 >> Child3 >> Child4 >> Child5; 
      ndecayproducts = abs(ndecayproducts);
      if (ndecayproducts==1) { continue;}
      else if (ndecayproducts==2) {         if (verbose) {printf("%9d\t%9d\t%9d\t%f\n", Father, Child1, Child2, BranchingRatio);}
      }
      else if (ndecayproducts==3) {
        if (verbose) { printf("%9d\t%9d\t%9d\t%9d\t%f\n", Father, Child1, Child2, Child3, BranchingRatio);}
      }
      else if (ndecayproducts==4) {
        //cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : skipping 4 particle decay! " << endl;
        //printf("%9d\t%9d\t%9d\t%9d\t%9d\t%f\n", Father, Child1, Child2, Child3, Child4, BranchingRatio);
        continue;
      }
      else { 
        cerr << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : unrecognised number of decays! " << ndecayproducts << endl;
        cout << line << endl;
        exit(EXIT_FAILURE);}
      for (int pi=0; pi<ns; pi++){
        if (Father==pdg[pi]){
          FathernotQuasiStable=false;
        }
      }
      if (ndecayproducts==2)
      {

        double Ma = fastreso.getParticleByPDG(Father)->getM();
        double Mb = fastreso.getParticleByPDG(Child1)->getM();
        double Mc = fastreso.getParticleByPDG(Child2)->getM();
        // Skip decays violating mass inequality (only possible by including resonance widths)
        if (Ma < Mb+Mc) {
          cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc << endl;
          cout << fastreso.getParticleByPDG(Father)->getName() << "\t";
          cout << fastreso.getParticleByPDG(Child1)->getName() << "\t";
          cout << fastreso.getParticleByPDG(Child2)->getName() << endl;
        continue;
        }

      count_decays++;
count_2decays++;
        //// Skip decays whose parents have decay width less than 10KeV 
        //if (GammaA <  10e-6) {
        //  cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle width smaller than 10 KeV ! Gamma = " << GammaA << " GeV" << endl;
        //  continue;
        //}

        double QS = fastreso.getParticleByPDG(Father)->getQS();
        double QS1 = fastreso.getParticleByPDG(Child1)->getQS();
        double QS2 = fastreso.getParticleByPDG(Child2)->getQS();
        if (QS != QS1+QS2){
          cout << fastreso.getParticleByPDG(Father)->getName() << endl;
          exit(EXIT_FAILURE);
        }

        if (FathernotQuasiStable){
          for (int pi=0; pi<ns; pi++){
            partial_charges[ fastreso.getParticleIndexByPDG(Father) ][pi] += BranchingRatio*partial_charges[ fastreso.getParticleIndexByPDG(Child1) ][pi];
            partial_charges[ fastreso.getParticleIndexByPDG(Father) ][pi] += BranchingRatio*partial_charges[ fastreso.getParticleIndexByPDG(Child2) ][pi];

          }
        }


      }
      else if (ndecayproducts==3)
      {

        double Ma = fastreso.getParticleByPDG(Father)->getM();
        double Mb = fastreso.getParticleByPDG(Child1)->getM();
        double Mc = fastreso.getParticleByPDG(Child2)->getM();
        double Md = fastreso.getParticleByPDG(Child3)->getM();
        if (Ma < Mb+Mc+Md) {
          cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc+Md << endl;
          continue;
        }

      count_decays++;
count_3decays++;
        double QS = fastreso.getParticleByPDG(Father)->getQS();
        double QS1 = fastreso.getParticleByPDG(Child1)->getQS();
        double QS2 = fastreso.getParticleByPDG(Child2)->getQS();
        double QS3 = fastreso.getParticleByPDG(Child3)->getQS();
        if (QS != QS1+QS2+QS3){
          cout << fastreso.getParticleByPDG(Father)->getName() << endl;
          exit(EXIT_FAILURE);
        }

        if (FathernotQuasiStable){
          for (int pi=0; pi<ns; pi++){
            partial_charges[ fastreso.getParticleIndexByPDG(Father) ][pi] += BranchingRatio*partial_charges[ fastreso.getParticleIndexByPDG(Child1) ][pi];
            partial_charges[ fastreso.getParticleIndexByPDG(Father) ][pi] += BranchingRatio*partial_charges[ fastreso.getParticleIndexByPDG(Child2) ][pi];
            partial_charges[ fastreso.getParticleIndexByPDG(Father) ][pi] += BranchingRatio*partial_charges[ fastreso.getParticleIndexByPDG(Child3) ][pi];
          }
        }
      }
      else if (ndecayproducts==4)
      {
        double Ma = fastreso.getParticleByPDG(Father)->getM();
        //double QBa = fastreso.getParticleByPDG(Father)->getQB();
        double Mb = fastreso.getParticleByPDG(Child1)->getM();
        double Mc = fastreso.getParticleByPDG(Child2)->getM();
        double Md = fastreso.getParticleByPDG(Child3)->getM();
        double Me = fastreso.getParticleByPDG(Child4)->getM();
        if (Ma < Mb+Mc+Md+Me) {
          cerr  << "\033[1mTFastReso_AZYHYDRO.cpp::do_decays\033[0m : \033[1;31merror\033[0m : father particle lighter than childern! " << Ma << " < " << Mb+Mc+Md << endl;
          continue;
        }

        double QS = fastreso.getParticleByPDG(Father)->getQS();
        double QS1 = fastreso.getParticleByPDG(Child1)->getQS();
        double QS2 = fastreso.getParticleByPDG(Child2)->getQS();
        double QS3 = fastreso.getParticleByPDG(Child3)->getQS();
        double QS4 = fastreso.getParticleByPDG(Child4)->getQS();
        if (QS != QS1+QS2+QS3+QS4){
          cout << fastreso.getParticleByPDG(Father)->getName() << endl;
          exit(EXIT_FAILURE);
        }

        if (FathernotQuasiStable){
          for (int pi=0; pi<ns; pi++){
            partial_charges[ fastreso.getParticleIndexByPDG(Father) ][pi] += BranchingRatio*partial_charges[ fastreso.getParticleIndexByPDG(Child1) ][pi];
            partial_charges[ fastreso.getParticleIndexByPDG(Father) ][pi] += BranchingRatio*partial_charges[ fastreso.getParticleIndexByPDG(Child2) ][pi];
            partial_charges[ fastreso.getParticleIndexByPDG(Father) ][pi] += BranchingRatio*partial_charges[ fastreso.getParticleIndexByPDG(Child3) ][pi];
            partial_charges[ fastreso.getParticleIndexByPDG(Father) ][pi] += BranchingRatio*partial_charges[ fastreso.getParticleIndexByPDG(Child4) ][pi];
          }
        }

      }


      }
    }


  fInputFile.close();
    if (verbose) {
      clog << "Total number of decays read " << count_decays << endl;
      clog << "Number of 2-body decays read " << count_2decays << endl;
      clog << "Number of 3-body decays read " << count_3decays << endl;
    }


  // Write the partial charge table
  const int Npart = fParticleNumber;
  FILE *pFile = fopen ((tag+"partial_charges.data").c_str(),"w");
  fprintf(pFile,"# Name\t");
  fprintf(pFile,"total\t");
  for (int pi=0; pi<ns; pi++){

    fprintf(pFile,"%d\t", pdg[pi]);
  }
  fprintf(pFile,"\n");
  for (int i=0; i<Npart; i++){
    fprintf(pFile, "%s\t",(fastreso.getParticle(i)->getName()).c_str());
    double tot=0;
    for (int pi=0; pi<ns; pi++){
      tot +=partial_charges[i][pi];
    }

    fprintf(pFile, "%3.2f\t",tot);
    for (int pi=0; pi<ns; pi++){
      if  ( partial_charges[i][pi]>0){
        fprintf(pFile, "%3.2f\t",partial_charges[i][pi]);
      } else {
        fprintf(pFile,"   \t");
      }
    }
    fprintf(pFile,"\n");
  }

  fclose(pFile);
  // write out ordinary HRG equation of state
  pFile = fopen ((tag+"eos.dat").c_str(),"w");
  fprintf(pFile,"# 1:T [GeV] \t");
  fprintf(pFile," 2:e/T^4\t");
  fprintf(pFile," 3:p/T^4\t");
  fprintf(pFile," 4:n/T^3\t");
  fprintf(pFile," 5:s/T^3\n");


 
  for (int j=0; j<NT; j++)
  {
    double   T = Tmin+dT*j;
    double e=0;
    double s=0;
    double n=0;
    double p=0;
    for (int i=0; i<Npart; i++){
      if (fastreso.getParticle(i)->getM() < 0.01) continue;
      int sign = ( fastreso.getParticle(i)->getType()==EParticleType::kBoson ? 1 : -1  );
      e += fastreso.getParticle(i)->getNu()*efunc( fastreso.getParticle(i)->getM()/T, 0.0, sign);
      s += fastreso.getParticle(i)->getNu()*sfunc( fastreso.getParticle(i)->getM()/T, 0.0, sign);
      p += fastreso.getParticle(i)->getNu()*pfunc( fastreso.getParticle(i)->getM()/T, 0.0, sign);
      n += fastreso.getParticle(i)->getNu()*nfunc( fastreso.getParticle(i)->getM()/T, 0.0, sign);

    }
    fprintf(pFile,"%8f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\n", T, e, p, n,s);
  }

  fclose(pFile);





  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
  // Calculate partial chemical potentials. that is solve for
  //! Ni(Tchem)/s(Tchem) = ni(Tchem,mu)/s(Tkin, mu)
  //! where Ni = sumj b_{j->i} nj, is sum of particles decying to i with branching ratios.
  //! s is total entropy
  //! for finite chemical potential, the "charge" is calculated as a sum of branching ratios
  //! So QMu = sumi m_i b{j->i}
  //! One ten solves for mu_i.


  // Initialize solver
  double *ci = new double[ns];
  double *fi= new double[ns];
  double *mui= new double[ns];
  for (int pi=0; pi<ns; pi++){ ci[pi]=0; fi[pi]=0; mui[pi]=0.0; }

  // Find Ni/s at Tchem and initialize ci
  double s=0;
  for (int i=0; i<Npart; i++){
    if (fastreso.getParticle(i)->getM() < 0.01) continue;
    int sign = ( fastreso.getParticle(i)->getType()==EParticleType::kBoson ? 1 : -1  );
    s += fastreso.getParticle(i)->getNu()*sfunc( fastreso.getParticle(i)->getM()/Tchem, 0.0, sign);
  }
  for (int pi=0; pi<ns; pi++){
    for (int i=0; i<Npart; i++){
      if (fastreso.getParticle(i)->getM() < 0.01) continue;
      int sign = ( fastreso.getParticle(i)->getType()==EParticleType::kBoson ? 1 : -1  );
      ci[pi] += partial_charges[i][pi]*fastreso.getParticle(i)->getNu()*nfunc( fastreso.getParticle(i)->getM()/Tchem, 0.0, sign);
    }
  }
  for (int pi=0; pi<ns; pi++){ ci[pi]/=s; }

  // initialize mass, sign and dof tables
  double * masses = new double[Npart]; 
  int *sign = new int[Npart]; 
  int *deg = new int[Npart]; 
  for (int i=0; i<Npart; i++){
    masses[i] = fastreso.getParticle(i)->getM();
    sign[i] = (fastreso.getParticle(i)->getType()==EParticleType::kBoson ? 1 : -1 );
    deg[i] = fastreso.getParticle(i)->getNu();
  }

  const gsl_multiroot_fsolver_type *Type;
  gsl_multiroot_fsolver *sol;
  gsl_vector *x = gsl_vector_alloc (ns);


  // open equation of state and chemical potential files for writing: (but only if not doing decays)
  FILE * pFile2;
    pFile = fopen((tag+std::format("peos_Tchem_{:.3f}.dat", Tchem)).c_str(),"w");
    pFile2 = fopen((tag + std::format("mu_Tchem_{:.3f}.dat", Tchem)).c_str(), "w");

    
    fprintf(pFile,"# 1:T [GeV] \t");
    fprintf(pFile," 2:e/T^4\t");
    fprintf(pFile," 3:p/T^4\t");
    fprintf(pFile," 4:n/T^3\t");
    fprintf(pFile," 5:s/T^3\n");

    fprintf(pFile2,"# 1:T [GeV]");
    for (int pi=0; pi<ns; pi++){
      fprintf(pFile2," \t%d:%d", pi+2, pdg[ns-1-pi]);
    }
    fprintf(pFile2,"\n");

  // Table for partial chemical potential for each particle
  double **QMu = new double *[fParticleNumber];
  for (int i=0; i<fParticleNumber; i++){
    QMu[i] = new double [NT];
    for (int pi=0; pi<NT; pi++){
      QMu[i][pi] =0.0;
    }
  }
  // loop over temperatures and solve for mui and compute thermodynamics

  for (int j=0; j<NT; j++){
    double Tkin= Tmax-dT*j;
    int Nl;

    for (int pi=0; pi<ns; pi++){
      gsl_vector_set (x, pi, 0.0);
    }
    Nl = 10; //round(abs(Tchem-Tkin)/0.001)+1;
    double dTkin = (Tkin-Tchem)/(Nl);

    int status;
    size_t iter = 0;
    for (int l=1; l<Nl+1; l++)
    {
      double Tkini=Tchem + dTkin*l;
      // initialize solver function
      struct rparams pr = {partial_charges, masses, sign, deg, ci, Tkini, Npart, ns};
      gsl_multiroot_function f = {&partial_f, ns, &pr};

      // re-use previous solution of mui to speed up calculation
      for (int pi=0; pi<ns; pi++){
        if (abs(Tkini-Tchem)<0.005){
          gsl_vector_set (x, pi, 0.0);
        } else {

          gsl_vector_set (x, pi, mui[pi]);
        }
      }
      // set up solver itself
      Type = gsl_multiroot_fsolver_hybrids;
      sol = gsl_multiroot_fsolver_alloc (Type, ns);
      gsl_multiroot_fsolver_set (sol, &f, x);

      // itterated until close to solution
      do
      {
        iter++;
        status = gsl_multiroot_fsolver_iterate (sol);


        if (status)   /* check if solver is stuck */
          break;
        status = gsl_multiroot_test_residual (sol->f, 1e-10);
      }
      while (status == GSL_CONTINUE && iter < 1000);


      // extract solution for chemical potential
      for (int pi=0; pi<ns; pi++){ mui[pi] = gsl_vector_get (sol->x, pi);}

    } 
    cerr << "Tkin  " << Tkin << " iter = " << iter << " " << gsl_strerror (status) << endl;

      fprintf(pFile2,"%f", Tkin );
      for (int pi=0; pi<ns; pi++){
        fprintf(pFile2," \t%g",mui[ns-1-pi]);
      }
      fprintf(pFile2,"\n");
      //calculate equation of state
      double e=0;
      double s=0;
      double n=0;
      double p=0;
      for (int i=0; i<Npart; i++){

        if (fastreso.getParticle(i)->getM() < 0.01) continue;
        double QMu=0;
        for (int pj=0; pj<ns; pj++){
          QMu +=partial_charges[i][pj]*mui[pj];
        }
        int sign = ( fastreso.getParticle(i)->getType()==EParticleType::kBoson ? 1 : -1  );
        e += fastreso.getParticle(i)->getNu()*efunc( fastreso.getParticle(i)->getM()/Tkin, QMu/Tkin, sign);
        s += fastreso.getParticle(i)->getNu()*sfunc( fastreso.getParticle(i)->getM()/Tkin, QMu/Tkin, sign);
        p += fastreso.getParticle(i)->getNu()*pfunc( fastreso.getParticle(i)->getM()/Tkin, QMu/Tkin, sign);
        n += fastreso.getParticle(i)->getNu()*nfunc( fastreso.getParticle(i)->getM()/Tkin, QMu/Tkin, sign);

      }
      fprintf(pFile,"%8f\t%15.9e\t%15.9e\t%15.9e\t%15.9e\n", Tkin, e, p, n,s);

    for (int i=0; i<Npart; i++){
      if (fastreso.getParticle(i)->getM() < 0.01) continue;
      for (int pj=0; pj<ns; pj++){
        QMu[i][j] +=partial_charges[i][pj]*mui[pj];
      }
    }
  }

  
  // close files and clean up memory
  gsl_vector_free(x);
  gsl_multiroot_fsolver_free(sol);

  delete [] ci;
  delete [] fi;
  delete [] mui;
  delete [] masses;
  delete [] sign;
  delete [] deg;

  for (int i=0; i<fParticleNumber; i++){
    delete [] partial_charges[i];
  }
  delete [] partial_charges;

  fclose(pFile);
  fclose(pFile2);

  
  cout << NT << "\n";
  #pragma omp parallel for
  for (int j=0; j<NT; j++){
    double Tkin = Tmin + dT*j;
    TFastReso_AZYHYDRO fastreso2;
    std::cout << "Thread " << omp_get_thread_num() << " processing j = " << j << std::endl;
    // read all particles
    fastreso2.read_particles_data(pdata);
    // make tag for spectra
    char buffer [50];
    sprintf (buffer, "%.3f", Tkin);
    string Ttag(buffer);
    sprintf (buffer, "%.3f", Tchem);
    string Tctag(buffer);
    if (dothermal){
      cout << " Do thermal T kin = " << Tkin  <<" GeV" <<endl;
      // initialize spectra with finite chemical potential
      for (int i=0; i<Npart; i++){
        if (fastreso2.getParticle(i)->getM() < 0.01) continue;
        fastreso2.do_thermal(i, Tkin, QMu[i][j]);
      }
      // print out thermal irreducible functions of some particles
      for (int k =0; k<ns2; k++){
        sprintf (buffer, "PDGid_%d", outpdg[k]);
        string name(buffer);
        fastreso2.getParticleByPDG(outpdg[k])->print(tag+name+"_thermal_Tchem_"+Tctag+"_Tkin_"+Ttag); }
      }
    
    
    
    cout << " Do decays T kin = " <<Tkin  <<" GeV" <<endl;
    fastreso2.do_decays(ddata);
    cout << " Finished decays T kin = " << Tkin  <<" GeV" <<endl;
    //// print out final ireducible functions
    for (int k =0; k<ns2; k++){
      sprintf (buffer, "PDGid_%d", outpdg[k]);
      string name(buffer);
      fastreso2.getParticleByPDG(outpdg[k])->print(tag+name+"_total_Tchem_"+Tctag+"_Tkin_"+Ttag); 
    }
  } // end of parallel for loop

  for (int i=0; i<fParticleNumber; i++){
  delete [] QMu[i];
}
delete [] QMu;

}







