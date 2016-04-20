/// This is a trivial script for analyzing the timing resolution for DT5720 data for 3" Hamamatsu PMT.
/// We fit the pulses for the monitor and 3" PMT and calculate the time difference between them.
/// 
/// The difficult part of this script was tweaking the parameters for the
/// Exponentially modified Gaussian for the different run conditions
/// This particular script is for PMT running at 8e6 gain, with no shapers.
///
/// T. Lindner , April 2016
/// 

double funcEMG(double* x, double* p){
  //Exponential gaussian used for fitting
  // p[0]: amplitude
  // p[1]: gaussian mu
  // p[2]: gaussian sig
  // p[3]: exponential decay constant  -  ???
  // p[4]: baseline

  double y = p[4] + (p[0]/0.3)*(p[3]/2.)*exp((p[1]+p[2]*p[2]/p[3]/2.-x[0])/(p[3]))*
    TMath::Erfc((p[1]+p[2]*p[2]/p[3] -x[0])/sqrt(2.)/p[2]);

  return y;

}



// Generalizede Normal Version 2 (Skewed Gaussian)

double funcGNV2(double* x, double* p) {
  //http://en.wikipedia.org/wiki/Generalized_normal_distribution   -- version 2
  //p[3] = baseline
  //p[4] = mag
  //p[5] = sigma
  //p[6] = start
  //p[7] = end
  

  double alpha, kappa, Xi, y, phi, range;  // scale, shape, location, --, standard normal gaussian, pulse region
  alpha = p[0];
  Xi = p[1];
  range = p[7]-p[6];
  kappa = p[2];

  // For Normal dist  (kappa = 0).
  if(kappa == 0.)
    y = (x[0] - Xi)/alpha;
  else
    y = - log(1 - kappa*(x[0] - Xi)/alpha) / kappa;  


  phi = exp(-y*y/p[5]) / sqrt(2*TMath::Pi());


  if(x[0] < (p[6]-range)  ||  x[0] > (p[7]+range)) {
    return p[3];  // return baseline; helps stabilize
  } else {
    return p[3] + p[4] * phi / (alpha - kappa*(x[0] - Xi))   *   exp(-0.5*pow((x[0]-Xi)/p[5],2)); 
    // Gaussian factor at end
  }


}


  

void findMin(TGraph *gr,double &min, double &min_bin){
  
  double x,y;
  for(int j = 0; j < gr->GetN(); j++){
    gr->GetPoint(j,x,y);
    if(y < min){
      min = y;
      min_bin = x;
    }
    
  }
  //  std::cout << "min " << min_bin << " " << min << std::endl;
}



void analyze_3inch_data(int nnnn){


  {


    TFile *f = new TFile("waveforms_r12199_1250V_dt720.root");

    TTree *waveforms = (TTree*)f->Get("waveformTree");

    if(waveforms)
      std::cout << "N entries: " << waveforms->GetEntries() << std::endl;
    TCanvas *c3 = new TCanvas("C3");
 
    TGraph *gr1 = new TGraph();
    waveforms->SetBranchAddress("g_ch1",&gr1);

    TGraph *gr2 = new TGraph();
    waveforms->SetBranchAddress("g_ch2",&gr2);


    TH1F *tdiff = new TH1F("diff","diff",200,0,100);
    TH1F *tdifflow = new TH1F("difflow","diff 1PE",200,0,100);
    TH1F *tdiffhigh = new TH1F("diffhigh","diff >1PE",200,0,100);

    TH1F *ph = new TH1F("ph","pulse height",2000,0,200);
    
    if(nnnn == 0)
      nnnn = waveforms->GetEntries();
    
    for(int i = 1; i < nnnn; i++){

      waveforms->GetEntry(i);

      
      // Loop over channels
      double times[2] = {-9999,-9999};

      double pmtph = -999;
      
      for(int ch = 0; ch < 2; ch++){

        TGraph *gr;
        if(ch == 0) gr = gr1;
        else gr = gr2;
        
        // Find the minimum bin
        double min = 999999, min_bin = -1;
        findMin(gr,min,min_bin);
      
        
      
        if(ch == 0){
          // Currently only has range set up for EMG, but EMG works well
          TF1 *fitted_ref = new TF1("EMG",funcEMG,min_bin-18,min_bin+14,5);
          fitted_ref->SetParameter(0,-150);
          fitted_ref->SetParameter(1,min_bin);
          fitted_ref->FixParameter(2,2);
          fitted_ref->FixParameter(3,0.3);
          fitted_ref->FixParameter(4,2162);
          //Exponential gaussian used for fitting
          // p[0]: amplitude
          // p[1]: gaussian mu
          // p[2]: gaussian sig
          // p[3]: exponential decay constant  -  ???
          // p[4]: baseline    
          gr->Fit("EMG","Q","",min_bin-18,min_bin+14);
          
          times[ch] =  fitted_ref->GetParameter(1);
          fitted_ref->Draw("SAME");
          //sleep(5);
          //char c;
          //cin << c;
          delete fitted_ref;
        }else{

          if(min <= 3938){ 

            gr->Draw("AP*");
            //            std::cout << "Fitting ch2" << std::endl;
            // Currently only has range set up for EMG, but EMG works well
            TF1 *fitted_ref = new TF1("EMG",funcEMG,min_bin-18,min_bin+14,5);
            fitted_ref->SetParameter(0,-20);
            fitted_ref->SetParameter(1,min_bin);
            fitted_ref->FixParameter(2,2.5);
            fitted_ref->FixParameter(3,3);
            fitted_ref->FixParameter(4,3943);
            //Exponential gaussian used for fitting
            // p[0]: amplitude
            // p[1]: gaussian mu
            // p[2]: gaussian sig
            // p[3]: exponential decay constant  -  ???
            // p[4]: baseline    
            gr->Fit("EMG","Q","",min_bin-18,min_bin+14);
            
            //            std::cout << fitted_ref->GetParameter(1) << std::endl;
            fitted_ref->Draw("SAME");
            times[ch] =  fitted_ref->GetParameter(1);
            ph->Fill(-fitted_ref->GetParameter(0));
            pmtph =  -fitted_ref->GetParameter(0);
            //sleep(5);
            //char c;
            //cin << c;
            delete fitted_ref;
          }
        }
      }

      double diff = times[1] - times[0];
      //      std::cout << "diff " << diff << std::endl;
      if(times[1] > 0){
        tdiff->Fill(diff);
        if(pmtph < 6){
          tdifflow->Fill(diff);
        }else{
          tdiffhigh->Fill(diff);
        }

      }

    }


    TCanvas *c = new TCanvas("C2");
    tdiff->Draw();

    
    TCanvas *c2 = new TCanvas("tdiff split");
    tdifflow->Draw();
    tdiffhigh->Draw("SAME");
    tdiffhigh->SetLineColor(2);

    TCanvas *cph = new TCanvas("C pulse height");
    ph->Draw();
  }


}
