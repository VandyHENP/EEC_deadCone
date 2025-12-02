// Copyright (C) 2025 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Authors: Dag Gillberg <dag.gillberg@cern.ch>

// Keywords: root; event display

// This is a program that uses ROOT to visualize the particles produced by
// Pythia. Particles are drawn in (y, phi)-space to depict the E/p flow.
// A pdf file is produced with multiple pages showing WH->qqbb events.

#include <sys/stat.h>

#include "Pythia8/Pythia.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "TLine.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TRandom3.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Recluster.hh"

using namespace Pythia8;
using namespace fastjet;

// Example main program to draw an event display.

//Function that takes in a particle ID from a particle object, and returns a value of 0 if it is light flavor, 1 if it is a gluon, 2 if it is a charm quark, and 3 if it is a bottom quark.
int getPartSpecies(int partID){
  if(abs(partID) <= 3 && partID != 0) return 0;
  if(abs(partID) == 21) return 1;
  if(abs(partID) == 4) return 2;
  if(abs(partID) == 5) return 3;

  return -1;
}

//Most likely hadrons that include b quarks
vector<int> bHadronList = {511, 521, 553, 5122, 5112};

//Takes a particles ID and determines if it is a hardon contain a b quark (one of the ones in the list)
bool isBHadron(int parID)
{
    for(int i=0; i<(int)bHadronList.size(); i++)
    {
        if(abs(parID) == bHadronList[i]) return true;
    }
    return false;
}

int main(int argc, char *argv[])
{
  //Establishes all variables for different collision types and different jet cuts that will happen after jets are identified.
  int ECM = 200;
  int pTHat = 35;
  int nEvents = 1000;
  double jetRadius = 0.8;
  double detEta = 1.1;
  bool decayOn = false;
  int jobIndex = -1;
  bool MPIon = true;
  std::string outDir = "";

  //Checks to make sure all arguments are valid, and then updates the variables above accordingly.
  if (argc > 14)
  {
    std::cout<<"error: please remember to use this program as: FastJetEECLight -e ECM -p pTHat -n nEvents; thanks so much, and I know you will be able to figure this out!!";
    return 0;
  }
  for(int i = 1; i < argc; i++)
  {
    if (argv[i] == std::string("-e"))
    {
      ECM = std::atoi(argv[i + 1]);
      std::cout<<"You did it!!"<< ECM;
    }
    else if (argv[i] == std::string("-p"))
    {
      pTHat = std::atoi(argv[i + 1]);
      std::cout<<"You did again with momentum this time. Nice Job!!" << pTHat;
    }
    else if (argv[i] == std::string("-n"))
    {
      nEvents = std::atoi(argv[i + 1]);
      std::cout<<"Number of Itterations Time. Lets Goooooooo!!" << nEvents;
    }
    else if (argv[i] == std::string("decay_on"))
    {
      decayOn = true;
    }
    else if (argv[i] == std::string("-f"))
    {
      jobIndex = std::atoi(argv[i + 1]);
      std::cout<<"Job index" << jobIndex;
    }
    else if (argv[i] == std::string("-d"))
    {
      outDir = std::string(argv[i + 1]);
      std::cout << "outDirectory: " << outDir;
    }
    else if (argv[i] == std::string("-r"))
    {
      jetRadius = std::atof(argv[i + 1]);
      std::cout<<"we got the Jet radius of " << jetRadius;
    }
    else if (argv[i] == std::string("MPIoff"))
    {
      MPIon = false;
    }

    std::cout<<"["<<i<<"] "<<argv[i]<<"\n";
  }


  std::cout << "outDir: " << outDir << std::endl;

  struct stat sb;
//Checks to make sure the inputted directory for file saving is real
  bool dirExists = (stat(Form("/data/rke_group/millsh1/pp/%s",outDir.c_str()),&sb) == 0);
  if(!dirExists){
    //Change /data/rke_group/millsh1/pp to wherever all of your data will be stored (do this everywhere this appears).
    std::cout << "Error, output directory " << Form("/data/rke_group/millsh1/pp/%s",outDir.c_str()) << " does not exist. Exiting" << std::endl;
    return 0;
  }
  std::cout << "directory exists! " << Form("/data/rke_group/millsh1/pp/%s",outDir.c_str()) << std::endl;
  //All ranges of pT tested based throughout the investigation, add pT mins and maxes if you would like
  int min_pT[8] = {10, 20, 30, 40, 60, 100, 150, 500};
  int max_pT[8] = {20, 30, 40, 60, 100, 150, 200, 550};

  int pTHatIndex = -1;
  //Takes the pTHat inputted and if it is valid for one of the pT ranges above, sets the index that then refers to those pT ranges to that value.
  switch(pTHat){
    case 5:
      pTHatIndex = 0;
      std::cout<<"Valid pTHat 5 inputted. Nice Job!!";
      break;
    case 15:
      pTHatIndex = 1;
      std::cout<<"Valid pTHat 15 inputted. Nice Job!!";
      break;
    case 25:
      pTHatIndex = 2;
      std::cout<<"Valid pTHat 25 inputted. Nice Job!!";
      break;
    case 35:
      pTHatIndex = 3;
      std::cout<<"Valid pTHat 35 inputted. Nice Job!!";
      break;
    case 45:
      pTHatIndex = 4;
      std::cout<<"Valid pTHat 45 inputted. Nice Job!!";
      break;
    case 80:
      pTHatIndex = 5;
      std::cout<<"Valid pTHat 80 inputted. Nice Job!!";
      break;
    case 125:
      pTHatIndex = 6;
      std::cout<<"Valid pTHat 125 inputted. Nice Job!!";
      break;
    case 450:
      pTHatIndex = 7;
      std::cout<<"Valid pTHat 450 inputted. Nice Job!!";
      break;
    default:
      break;

  }
  //Stops the code if an invalid pTHat is inputted
  if(pTHatIndex == -1){
    std::cout << "incorrect pTHat value of " << pTHat << " passed. Try a valid value. Exiting" << std::endl;
    return 2;
  }


  // Adjust ROOTs display options.
  gStyle->SetOptStat(0);
  //Allows to call with the getPartSpecies first letter of each species for easy naming.
  string partSpeciesNames[4] = {"l","g","c","b"};

  //Randomly generates a seed upon which event generation will be based.
  TRandom3 *r = new TRandom3();
  r->SetSeed(0);
  unsigned int seed = r->GetSeed();
  while(seed >= 900000000){
    seed = seed/10;
  }

  // Generator.
  Pythia pythia;

  // PYTHIA setup. Process selection. RHIC initialization.
  pythia.readString(Form("Beams:eCM = %d",ECM));
  pythia.readString("HardQCD:all = on");
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Random:setSeed = on");
  pythia.readString(Form("Random:seed = %d",seed));
  pythia.readString(Form("PhaseSpace:pTHatMin = %d", pTHat));
  //Turns on and off MPI based on user input above
  if (!MPIon)
  {
    pythia.readString("PartonLevel:MPI = off");
  }
  if (!decayOn)
  {
    //charm hadrons
    /*
    pythia.readString("411:mayDecay = off");
    pythia.readString("-411:mayDecay = off");
    pythia.readString("421:mayDecay = off");
    pythia.readString("-421:mayDecay = off");
    pythia.readString("443:mayDecay = off");
    pythia.readString("4122:mayDecay = off");
    pythia.readString("100443:mayDecay = off");
    */

    //bottom hardons
    for(int b=0; b<(int)bHadronList.size(); b++)
    {
        pythia.readString(Form("%d:mayDecay = off",bHadronList[b]));
        pythia.readString(Form("-%d:mayDecay = off",bHadronList[b]));
    }
    //pythia.readString("511:mayDecay = off");
    //pythia.readString("-511:mayDecay = off");
    //pythia.readString("521:mayDecay = off");
    //pythia.readString("-521:mayDecay = off");
    //pythia.readString("553:mayDecay = off");
    //pythia.readString("5122:mayDecay = off");
    //pythia.readString("100553:mayDecay = off");
    //pythia.readString("200553:mayDecay = off");
  }



  // If Pythia fails to initialize, exit with error.
  if (!pythia.init())
    return 1;

  //Set the parameters for the histogram to be filled by the EEC
  //Number of Bins
  const int nbins = 30;
  //Smallest bin min
  double mindr = 0.01;
  //Largest Bin max
  double maxdr = jetRadius*2;
  double drwidth = (log(maxdr)- log(mindr))/nbins;
  double bins[nbins + 1];

  //Creates an array of bins logrithmically sized
  for (int i = 0; i <= nbins; i++)
  {
    bins[i] =  mindr * exp(i*drwidth);
  }
  
//Make Light Quark Jet Histogram
  TFile *f = new TFile(Form("/data/rke_group/millsh1/pp/%s/FJEEC_noEdgeEffects_eCM%d_jetRadius%.1f_pTHat%d_decays%s_MPI%s_n%dk_%d.root",outDir.c_str(),ECM, jetRadius,pTHat,(decayOn ? "On" : "Off"),(MPIon ? "On" : "Off"),nEvents/1000,jobIndex), "RECREATE");
  //Makes histograms that will house different EECs
  TH1D *EEC_hist[4];
  TH1D *EEC_hist_BMeson[4];
  TH1D *EEC_hist_WTA[4];
  TH1D *jetSpec[4];
  TH1D *nBJets = new TH1D("nBJets","Number of B Jets",5,-0.5,4.5);
  //Creates a histogram for each flavor of quark
  for(int i=0; i<4; i++){
    EEC_hist[i] = new TH1D(Form("EEC_hist_%s",partSpeciesNames[i].c_str()),Form("%s EEC Histogram;#DeltaR;EEC",partSpeciesNames[i].c_str()), nbins, bins);
    EEC_hist_BMeson[i] = new TH1D(Form("EEC_hist_BMeson_%s",partSpeciesNames[i].c_str()),Form("%s EEC Histogram B Meson;#DeltaR;EEC",partSpeciesNames[i].c_str()), nbins, bins);
    EEC_hist_WTA[i] = new TH1D(Form("EEC_hist_WTA_%s",partSpeciesNames[i].c_str()),Form("%s EEC Histogram WTA;#DeltaR;EEC",partSpeciesNames[i].c_str()), nbins, bins);
    jetSpec[i] = new TH1D(Form("jetSpec_%s",partSpeciesNames[i].c_str()),Form("%s Jet Spectrum;p_{T} [GeV/c];Number of Jets",partSpeciesNames[i].c_str()), 100,0.0,100.0);
  }

  //TH1D *EEC_hist = new TH1D("EEC_hist", "Light EEC Histogram;#DeltaR;EEC", nbins, bins);
  //TH1D *jetSpec = new TH1D("jetSpec","Jet Spectrum;#DeltaR;EEC", nbins, bins);

//Make Gluon Jet Histogram
  //TFile *g = new TFile(Form("FJEECGluonFile_eCM%d_pTHat%d_decays%s_n%dk.root", ECM,pTHat,(decayOn ? "On" : "Off"),nEvents/1000), "RECREATE");
  //TH1D *Gluon_hist = new TH1D("Gluon_hist", "Gluon EEC Histogram;#DeltaR;EEC", nbins, bins);
  //TH1D *jetSpecOG = new TH1D("jetSpecOG", "Jet Spectrum;DeltaR;EEC", nbins, bins);

//Create Pythia Event in which collisions are run
  Event &event = pythia.event;

//Define the jets that we are looking for using antikt algorithm (highest momentum particle down)
  JetDefinition jetdef(antikt_algorithm, jetRadius);
  JetDefinition wtadef(antikt_algorithm, jetRadius, WTA_pt_scheme);
  Selector jetSel=SelectorAbsEtaMax(detEta - jetRadius);

  Recluster reclust(wtadef);

//Set up vectors to store different particle sets as we sort
  std::vector<PseudoJet> particles;
  std::vector<PseudoJet> ch_particles;
  std::vector<PseudoJet> constituents;
  std::vector<double> side_lengths;


//Loop through the number of events selcted by the user only selecting events that meet the requirement of a leading quark type (u,d,s)
  int iEvent = 0;
  int cEvent = 0;
  int bEvent = 0;
  int gEvent = 0;
  while(iEvent<nEvents || cEvent<nEvents || bEvent<nEvents || gEvent<nEvents)
  //for (int iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    //Clear vectors for comparison
    particles.clear();
    ch_particles.clear();

    // Generate event. (Skip to next if pythia.next() returns false = error.)
    if (!pythia.next())
    {
      continue;
    }


    
    //selects u, d, and s leading particles and does not increment events if this is not met.
    if( (abs(event[5].id()) > 5 && abs(event[5].id()) != 21) || (abs(event[6].id()) > 5 && abs(event[6].id()) != 21) ){
      continue;
    }
    //Tells us how many events by the thousand the program has actually run (Different than pythia output as we reject the events without u, d, s leading quarks)
    if(bEvent % 1000 == 0){
      std::cout << "working on event " << bEvent << std::endl;
    }

    bool hasBMeson = false;

    //Iterates through every particle in the selected event 
    for (int j = 0; j < event.size(); j++)
    { 
      //Checks if the particle is done doing radiations and decay    
      if (event[j].isFinal() && abs(event[j].eta()) <= detEta)
      {
        PseudoJet temp(event[j].px(), event[j].py(), event[j].pz(), event[j].e());
        //Stores the particle IDs in a seperate array.
        temp.set_user_index(event[j].id());
        particles.push_back(temp);
        //Checks if particles are charged or if it is a BHadron
        if(!hasBMeson && isBHadron(event[j].id())) hasBMeson = true;

        //if it is, it gets added to the vector of particles we will work with
        if (event[j].isCharged() || isBHadron(event[j].id()))
        {
          ch_particles.push_back(temp);
        }
      }
    }
    //Determines what flavor the initiating quarks are and then itterates the corresponding event variable
    if (abs(event[5].id()) < 4 || abs(event[6].id()) < 4)
    {
      iEvent++;
    }
    if (abs(event[5].id()) == 4 || abs(event[6].id()) == 4)
    {
      cEvent++;
    }
    if (abs(event[5].id()) == 5 || abs(event[6].id()) == 5)
    {
      bEvent++;
    }
    if (abs(event[5].id()) == 21 || abs(event[6].id()) == 21)
    {
      gEvent++;
    }
    //Creates the jets based on predefined algorithms
      ClusterSequence cs(particles, jetdef);
      std::vector<PseudoJet> jets = sorted_by_pt(jetSel(cs.inclusive_jets()));
      PseudoJet part5(event[5].px(), event[5].py(), event[5].pz(), event[5].e());
      PseudoJet part6(event[6].px(), event[6].py(), event[6].pz(), event[6].e());
      //Matches the initiating parton with a jet
      for (auto &jet : jets)
      {
        //checks if the jet falls into the specified pT range 
        if(jet.pt() > min_pT[pTHatIndex] && jet.pt() < max_pT[pTHatIndex])
        {
          int partSpecies = -1;

          if (part5.delta_R(jet) < 0.75*jetRadius)
          {
            partSpecies = getPartSpecies(event[5].id());
          }
          else if (part6.delta_R(jet) < 0.75*jetRadius)
          {
            partSpecies = getPartSpecies(event[6].id());
          }

          if(partSpecies == -1) continue;

          if(partSpecies == 3) nBJets->Fill(0);

          //If there are B mesons in the jet, it gets added to a special list

          vector<int> bMesonsInJet;
          if(hasBMeson)
          {
            for(int jc=0; jc<(int)jet.constituents().size(); jc++)
            {
                if(isBHadron(jet.constituents()[jc].user_index()))
                {
                    bMesonsInJet.push_back(jc);
                }
            }
          }
          //Checks if there is a Bmeson in the jets and adds 1 to corresponding categories if that is the case.
          if(bMesonsInJet.size() > 0) nBJets->Fill(1);
          if(bMesonsInJet.size() > 0 && partSpecies == 3) nBJets->Fill(3);

          //Creates a jet spectra for the different particle species.
          jetSpec[partSpecies]->Fill(jet.pt());
          //Conducts the EEC on all particles within the jet radius
          constituents.clear();
          for (auto &part : ch_particles)
          {
            if (part.delta_R(jet) < (jetRadius * 2))
            {
              constituents.push_back(part);
            }
          }
          //Iterates though every particle in the jet radius and fills a histogram with its distnce from every other particle weighted by the product of their momentums.
          for (unsigned int p = 0; p < constituents.size(); p++)
          {
            for (unsigned int q = p + 1; q < constituents.size(); q++)
            {
              double dr = constituents[p].delta_R(constituents[q]);
              EEC_hist[partSpecies]->Fill(dr, constituents[p].pt() * constituents[q].pt());
              if(bMesonsInJet.size() > 0)
              {
                EEC_hist_BMeson[partSpecies]->Fill(dr, constituents[p].pt() * constituents[q].pt());
              }
            }
          }
          //reculsters all Bmeson included jets for Winner Takes All
          if(bMesonsInJet.size() > 0)
          {
            //for WTA EEC
            PseudoJet wtaJet = reclust(jet);
            bool isBwta = false;
            for(int bMes=0; bMes<(int)bMesonsInJet.size(); bMes++)
            {
                if(jet.constituents()[bMesonsInJet[bMes]].delta_R(wtaJet) == 0.0)
                //if(jet.constituents()[bMesonsInJet[bMes]].eta() == wtaJet.eta() && jet.constituents()[bMesonsInJet[bMes]].phi() == wtaJet.phi())
                {
                    isBwta = true;
                    break;
                }
            }

            if(isBwta) nBJets->Fill(2);
            if(isBwta && partSpecies == 3) nBJets->Fill(4);

            //Runs seperarte EEC on WTA reclustered jets
            if(isBwta)
            {
                constituents.clear();
                for (auto &part : ch_particles)
                {
                    if (part.delta_R(wtaJet) < (jetRadius * 2))
                    {
                    constituents.push_back(part);
                    }
                }
                for (unsigned int p = 0; p < constituents.size(); p++)
                {
                    for (unsigned int q = p + 1; q < constituents.size(); q++)
                    {
                        double dr = constituents[p].delta_R(constituents[q]);
                        EEC_hist_WTA[partSpecies]->Fill(dr, constituents[p].pt() * constituents[q].pt());
                    }
                }
            }
          }


        }
      }

   
  }

  
  //TCanvas *c1 = new TCanvas();
  //c1->SetLogy();
  //c1->SetLogx();

  //Saves the markers and their colors for each EEC histogram.
  int color[4] = {418, 617, 632, 600};
  int marker[4] = {20, 33, 21, 22};
  TCanvas *e3ccanvas[4];
  for(int i=0; i<4; i++)
  {
    //c1->Clear();
    EEC_hist[i]->SetLineColor(color[i]);
    EEC_hist[i]->SetMarkerStyle(marker[i]);
    EEC_hist[i]->SetMarkerColor(color[i]);
    
    EEC_hist_BMeson[i]->SetLineColor(color[i]);
    EEC_hist_BMeson[i]->SetMarkerStyle(marker[i]);
    EEC_hist_BMeson[i]->SetMarkerColor(color[i]);

    EEC_hist_WTA[i]->SetLineColor(color[i]);
    EEC_hist_WTA[i]->SetMarkerStyle(marker[i]);
    EEC_hist_WTA[i]->SetMarkerColor(color[i]);
    //EEC_hist[i]->Draw("p");
    //c1->SaveAs(Form("EECHist_%s.png",partSpeciesNames[i].c_str()));
  }
  
  //Bin labels for diagnostic histogram for Bmeson analysis
  string binLabels[5] = {"b quark","B meson","B WTA","b quark & B meson","b quark & B WTA"};

  for(int i=0; i<5; i++)
  {
    nBJets->GetXaxis()->SetBinLabel(i+1,binLabels[i].c_str());
  }
  //Saves all histograms to a specific file location specified above 
  f->cd();
  for(int i=0; i<4; i++){
    EEC_hist[i]->Write();
    EEC_hist_BMeson[i]->Write();
    EEC_hist_WTA[i]->Write();
    jetSpec[i]->Write();
  }
  nBJets->Write();
  f->Close();

  // Done.
  return 0;
}