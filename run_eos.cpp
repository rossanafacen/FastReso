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
#include <memory>
#include <string>
#include <fstream>
#include <cmath>
#include <sstream>
#include <TFastReso_AZYHYDRO.h>
#include <TKernel.h>
#include "eos_func.h"

using namespace std ;


//! This code computes equation of state for given content of resonances.
int main(int argc, const char * argv[])
{

  if (argc<3) {
    cerr << argv[0] << " not enough input arguments " << endl;
    cerr << argv[0] << " decays_PDG2016_massordered.dat ./outputfolder/" << endl;
    exit(EXIT_FAILURE) ;
  }

  // freeze-out temperature
  string pdata(argv[1]);// = "decays_PDG2016_massordered.dat";
  string ddata(argv[1]);// = "decays_PDG2016_massordered.dat";
  string tag(argv[2]);


  TFastReso_AZYHYDRO fastreso;
  // read all particles
  fastreso.read_particles_data(pdata);
  //   fastreso.print_partial_charges();
  int fParticleNumber = fastreso.fParticleNumber;
  double Tmax=0.2;
  double Tmin=0.0;
  double dT = 0.001;
  int NT= ((int) round ((Tmax-Tmin)/dT ))+1;
//  double dT=(Tmax-Tmin)/NT;

    FILE * pFile;
    pFile = fopen ((tag+"eos.dat").c_str(),"w");
    fprintf(pFile,"# 1:T [GeV] \t");
    fprintf(pFile," 2:e/T^4\t");
    fprintf(pFile," 3:p/T^4\t");
    fprintf(pFile," 4:n/T^3\t");
    fprintf(pFile," 5:s/T^3\n");



  for (int j=0; j<NT; j++){
    double   T= Tmax-dT*j;
    double e=0;
    double s=0;
    double n=0;
    double p=0;
    for (int i=0; i<fParticleNumber; i++){
      if (fastreso.getParticle(i)->getM() < 0.01) continue;
      int sign = ( fastreso.getParticle(i)->getType()==EParticleType::kBoson ? 1 : -1  );
      e += fastreso.getParticle(i)->getNu()*efunc( fastreso.getParticle(i)->getM()/T, 0.0, sign);
      s += fastreso.getParticle(i)->getNu()*sfunc( fastreso.getParticle(i)->getM()/T, 0.0, sign);
      p += fastreso.getParticle(i)->getNu()*pfunc( fastreso.getParticle(i)->getM()/T, 0.0, sign);
      n += fastreso.getParticle(i)->getNu()*nfunc( fastreso.getParticle(i)->getM()/T, 0.0, sign);

    }

    fprintf(pFile,"%8f\t%15.9f\t%15.9f\t%15.9f\t%15.9f\n", T, e, p, n,s);
//    cout <<  T << " " << e  << " " << p << " " << n << " " << s << endl;
  }

    fclose(pFile);
}



