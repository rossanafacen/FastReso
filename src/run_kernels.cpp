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
#include <string>
#include <TKernel.h>
#include <TParticle_AZYHYDRO.h>
#include <nlohmann/json.hpp>
#include <fstream>
#include <gsl/gsl_errno.h>
using namespace std ;
using json = nlohmann::json;

int main(int argc, const char * argv[])
{
  if (argc<5) {
    cerr << argv[0] << " not enough input arguments " << endl;
    cerr << argv[0] << " ./inputfolder/ ./outputfolder/ /json_file/ job_id" << endl;
    exit(EXIT_FAILURE) ;
  }

  string tag1(argv[1]);
  string tag2(argv[2]);
  string config(argv[3]); 
  int jobposition = std::stoi(argv[4]);

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
  bool diff_flag = j["physics"]["diffusion"];
   
  /////////////////// DANGER ///////////////////////////////////////////////////
  // If integration fails, you can turn the handler off. Integration functions
  // in TFastReso_AZYHYDRO.cpp and TFastReso_formulas.h will try reducing accuracy
  // for the failing cases.  Use with caution!
  gsl_set_error_handler_off ();
  /////////////////// DANGER ///////////////////////////////////////////////////

   

  double Tfo = Tkin; //GeV
  char buffer [50];
  sprintf (buffer, "%.3f", Tfo);
  string Ttag(buffer);    
  cout << " Do kernels = " << Ttag  <<" GeV" <<endl;
  
  if (dothermal){
    for (int k =0; k<ns; k++){
      string tag;  
      sprintf (buffer, "PDGid_%d", pdg[k]);
      string name(buffer);

      cout << " Do particle  " <<tag1+name+"_thermal_Tkin_"+Ttag  <<endl;
      tag = tag1+name+"_thermal_Tkin_"+Ttag ; 
      TParticle_AZYHYDRO particle(tag);    
      TKernel kernel(&particle);
      if (diff_flag) {
        kernel.print(tag2+name+"_thermal_Tkin_"+Ttag,"Keq Kshear Kbulk Kdiff");
      } else {
        kernel.print(tag2+name+"_thermal_Tkin_"+Ttag,"Keq Kshear Kbulk Ktemp Kvel");
      }       
      }
    }

  for (int k =0; k<ns; k++){
    sprintf (buffer, "PDGid_%d", pdg[k]);
    string name(buffer);
    string tag;
      
    cout << " Do particle  " <<tag1+name+"_total_Tkin_"+Ttag  <<endl;
    tag = tag1+name+"_total_Tkin_"+Ttag ; 
    TParticle_AZYHYDRO particle(tag);
    TKernel kernel(&particle);
    // if (diff_flag) {
    //       kernel.print(tag2+name+"_total_Tkin_"+Ttag,"Keq Kshear Kbulk Kdiff" );
    //     } else {
        kernel.print(tag2+name+"_total_Tkin_"+Ttag,"Keq Kshear Kbulk Ktemp Kvel" );
    //}
    }
}



