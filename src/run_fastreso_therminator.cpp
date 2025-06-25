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
#include <TFastReso_THERMINATOR.h>
#include <gsl/gsl_errno.h>

#include <nlohmann/json.hpp>
#include <fstream>

using namespace std ;
using json = nlohmann::json;

int main(int argc, const char * argv[])
{

  if (argc<6) {
    cerr << argv[0] << " not enough input arguments " << endl;
    cerr << argv[0] << " particles.data decays.data ./outputfolder/ /json_file/ job_id" << endl;
    exit(EXIT_FAILURE) ;
  }
/////////////////// DANGER ///////////////////////////////////////////////////
  // If integration fails, you can turn the handler off. Integration functions
  // in TFastReso_THERMINATOR.cpp and TFastReso_formulas.h will try reducing accuracy
  // for the failing cases.  Use with caution!
  gsl_set_error_handler_off ();
/////////////////// DANGER ///////////////////////////////////////////////////


  string pdata(argv[1]);// = "particles.data";
  string ddata(argv[2]);// = "decays_without_weak.data";
  string tag(argv[3]);  
  string config(argv[4]); 
  int jobposition = std::stoi(argv[5]);


  ifstream file(config);
  json j;
  file >> j;

  int ns = j["physics"]["n_particles"];
  int pdg[ns];
  for(int i = 0; i < ns; i++) {
       pdg[i] = j["physics"]["particles"][i];
  } 
  double Tkin = j["physics"]["Tkin"][jobposition];
  bool dothermal = j["physics"]["thermal"];
  bool dodiff = j["physics"]["diffusion"];
  cout << "Job position: " << jobposition << " with temperature "<< Tkin << endl;
  
  bool verbose = true ; // print what is being done
       
   {

        double Tfo=Tkin; //GeV freeze-out temperature
        char buffer [50];
        sprintf (buffer, "%.4f", Tfo);
        string Ttag(buffer);
        string inarg;
        if (dodiff) {
          inarg= "Feq Fshear Fbulk Fdiff"; // which components to print out
        } else {
          inarg= "Feq Fshear Fbulk"; // which components to print out
        }
        // create the TFastReso class
        TFastReso_THERMINATOR fastreso;
        // read all particles
        fastreso.read_particles_data(pdata,inarg, verbose);

        if (dothermal) {
          // initialize components with thermal distributions at freeze-out temperature
          cout << " Do thermal T = " <<Tfo  <<" GeV" <<endl;
          fastreso.do_thermal(Tfo);

          // print out particles
          for (int k =0; k<ns; k++){ 
          sprintf (buffer, "PDGid_%d", pdg[k]);
          string name(buffer);
            fastreso.getParticleByPDG(pdg[k])->print(tag+name+"_thermal_T"+Ttag); }
        }

        // perform decays
        cout << " Do decays T = " <<Tfo  <<" GeV" <<endl;
        fastreso.do_decays(ddata, verbose);

        // print out components
        for (int k =0; k<ns; k++){ 
        sprintf (buffer, "PDGid_%d", pdg[k]);
        string name(buffer);
          fastreso.getParticleByPDG(pdg[k])->print(tag+name+"_total_T"+Ttag); }
    }
}
