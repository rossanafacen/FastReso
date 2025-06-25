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
#include <TFastReso_AZYHYDRO.h>
#include <gsl/gsl_errno.h>
#include <nlohmann/json.hpp>
#include <fstream>

using namespace std ;
using json = nlohmann::json;
int main(int argc, const char * argv[])
{
    if (argc<5) {
        cerr << argv[0] << " not enough input arguments " << endl;
        cerr << argv[0] << " decays_PDG2016_massordered.dat ./outputfolder/ /json_file/ job_id" << endl;
        exit(EXIT_FAILURE) ;
    }
/////////////////// DANGER ///////////////////////////////////////////////////
  // If integration fails, you can turn the handler off. Integration functions
  // in TFastReso_THERMINATOR.cpp and TFastReso_formulas.h will try reducing accuracy
  // for the failing cases.  Use with caution!
  gsl_set_error_handler_off ();
/////////////////// DANGER ///////////////////////////////////////////////////



    // freeze-out temperature
    string pdata(argv[1]);// = "decays_PDG2016_massordered.dat";
    string ddata(argv[1]);// = "decays_PDG2016_massordered.dat";
    string tag(argv[2]);
    string config(argv[3]); 
    int jobposition = std::stoi(argv[4]);
    bool verbose = true ; // print what is being done
    
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
    int NT = 2;
    #pragma omp parallel for //to delete
    for (int j=0; j<NT; j++)
    {
        double Tfo=Tkin; //GeV freeze-out temperature
        char buffer [50];
        sprintf (buffer, "%.3f", Tfo);
        string Ttag(buffer);
        string inarg;
        if (dodiff) {
          inarg= "Feq Fshear Fbulk Ftemp Fvel Fdiff"; // which components to print out
        } else {
          inarg= "Feq Fshear Fbulk Ftemp Fvel"; 
        }

        // create the TFastReso class
        TFastReso_AZYHYDRO fastreso;
        // read all particles
        fastreso.read_particles_data(pdata,inarg, verbose);
        if(dothermal) {
            // initialize components with thermal distributions at freeze-out temperature
            cout << " Do thermal T = " <<Tfo  <<" GeV" <<endl;
            fastreso.do_thermal(Tfo);

            // print out particles
            for (int k =0; k<ns; k++){ 
            sprintf (buffer, "%d_PDGid_%d", j, pdg[k]);
            string name(buffer);
            fastreso.getParticleByPDG(pdg[k])->print(tag+name+"_thermal_Tkin_"+Ttag); }
        }

        // perform decays
        cout << " Do decays T = " << Tfo  <<" GeV" <<endl;
        fastreso.do_decays(ddata, verbose);

        // print out components
        for (int k =0; k<ns; k++){ 
        sprintf (buffer, "%d_PDGid_%d", j, pdg[k]);
        string name(buffer);
          fastreso.getParticleByPDG(pdg[k])->print(tag+name+"_total_Tkin_"+Ttag); }
        
    }
}
