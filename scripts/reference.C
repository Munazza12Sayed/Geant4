
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
#include <algorithm>    // std::max
#include <stdlib.h>
#include <math.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TTree.h>
#include <TCut.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TChain.h>
#include <TChainElement.h>

bool maxorder(const std::pair<Int_t,Int_t>  &v1, const std::pair<Int_t,Int_t> &v2) { return (v1.second > v2.second ); }


void m()
{

  TFile *f = new TFile("150GeV.root");
  TTree *t = (TTree*)f->Get("MyRun");

  Int_t Event;
  Int_t nSteps;
  Double_t EDepI_perstep[200];
  Double_t  posIAlStepx[200][200];
  Double_t  posIAlStepy[200][200];
  Double_t  posIAlStepz[200][200];
  Int_t Nclust_perstep[200];//allows for effective use of architecture with wider formats
 
  std::vector<std::pair<Int_t,Int_t> > StepMaxClust;

  int nentries;
  
  t->SetBranchAddress("Event",&Event);
  t->SetBranchAddress("nSteps",&nSteps);
  t->SetBranchAddress("EDepI_perstep",&EDepI_perstep);
  t->SetBranchAddress("posIAlStepx",&posIAlStepx);
  t->SetBranchAddress("posIAlStepy",&posIAlStepy);
  t->SetBranchAddress("posIAlStepz",&posIAlStepz);
  t->SetBranchAddress("Nclust_perstep",&Nclust_perstep);
  
  nentries = t->GetEntries();

    ofstream myfile;
  myfile.open ("tripleGEM_150GeV_0.txt");

  cout << "Event"<< setw(7)<<"  "<<"N.clusters"<< setw(7)<<"  "<<"Cluster #"<< setw(7)<<"   "<<"Edepionization"<<"     "<<"E/cluster"<<"      "<<"x-pos" << setw(7)<<"    "<< "y-pos"<<  setw(7)<<"    "<< "z-pos" << endl;
   
 myfile << "Event"<< setw(7)<<"  "<<"N.clusters"<< setw(7)<<"  "<<"Cluster #"<< setw(7)<<"   "<<"Edepionization"<<"     "<<"E/cluster"<<"      "<<"x-pos" << setw(7)<<"    "<< "y-pos"<<  setw(7)<<"    "<< "z-pos" << endl;

  //std::fstream ofs("3.txt",std::ofstream::out);

  for(int i=0;i<1000;i++){ //loop for number of events
    t->GetEntry(i);
    
    std::pair<Int_t,Int_t> temp;
    
    //    cout<<" Event   "<<i<<endl;

    for(int k=0;k<nSteps;k++){
           
      //      cout<<" step number  "<<k<<" cluster per step  "<<Nclust_perstep[k]<<endl;
      temp.first = k;
      temp.second = Nclust_perstep[k];     
      StepMaxClust.push_back(temp);

    }    //  loop for steps
    
    if(StepMaxClust.size()>0){
      
      std::sort(StepMaxClust.begin(),StepMaxClust.end(),maxorder);
      //      cout<<" Step with max cluster "<<StepMaxClust[0].first<<"  # of clusters   "<<StepMaxClust[0].second<<endl;

      Int_t nclust = Nclust_perstep[StepMaxClust[0].first];
      
      //      cout<<" number of clusters to loop  "<<nclust<<endl;
      
      cout<<"       "<<endl;
      
      for(int j=0;j<nclust;j++){
	cout << Event <<setw(12)<<""<< Nclust_perstep[StepMaxClust[0].first] <<setw(15)<<""<<j<<setw(15)<<"   "<<std::setprecision(6)<<EDepI_perstep[StepMaxClust[0].first]<<setw(10)<<""<<EDepI_perstep[StepMaxClust[0].first]/Nclust_perstep[StepMaxClust[0].first]<<setw(8)<<""<<fixed<< std::setprecision(3)<< posIAlStepx[StepMaxClust[0].first][j] <<setw(7)<<""<< posIAlStepy[StepMaxClust[0].first][j]<<setw(7)<<""<<posIAlStepz[StepMaxClust[0].first][j] << endl; 


	myfile << Event <<setw(12)<<""<< Nclust_perstep[StepMaxClust[0].first] <<setw(15)<<""<<j<<setw(15)<<"   "<<std::setprecision(6)<<EDepI_perstep[StepMaxClust[0].first]<<setw(10)<<""<<EDepI_perstep[StepMaxClust[0].first]/Nclust_perstep[StepMaxClust[0].first]<<setw(8)<<""<<fixed<< std::setprecision(3)<< posIAlStepx[StepMaxClust[0].first][j] <<setw(7)<<""<< posIAlStepy[StepMaxClust[0].first][j]<<setw(7)<<""<<posIAlStepz[StepMaxClust[0].first][j] << endl;      
      }
    }

    StepMaxClust.clear();
  
  } // loop for events
  
  myfile.close();
  //ofs.close();
}


