void draw_E2C_High_Energy_Bjet(int ECM = 200, double jetRadius = 0.4, bool decayOn = false, bool MPIon = true, int nEvents = 1, string outDir = "")
{
  //list of different flavors of leading parton that we will be looking for
    string partSpeciesNames[4] = {"Gluon","Light Flavor","Charm","Bottom"};
  
  //Set all of the style preferences and the matricies of parameters specified within the file name and the inputted values by the user
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
    double globalmin = 0.013581;
    double globalmax = 1.95972;
    int color[4] = {617, 418, 632, 600};
    int altcolor[4] = {797, 433, 625, 4};
    int marker[4] = {33, 20, 21, 22};
    int pTHat[7] = {5, 15, 25, 35, 45, 80, 120};
    int min_pT[8] = {10, 20, 30, 40, 60, 100, 150, 500};
    int max_pT[8] = {20, 30, 40, 60, 100, 150, 200, 550};
    int pTmidrange[4] = {15, 25, 35, 50};
    //Hard coded limits on ranges for plots, change these if looking at different pT ranges or collision energies as needed.
    double slopeMinRange[2] = {-1.6, -2.5};
    double slopeMaxRange[2] = {-0.3, -0.3};
    double peakMinRange[2] = {0, 0};
    double peakMaxRange[2] = {0.55, 0.4};
    //Allows for comparisons of both 0.4 and 0.8 radius jets.
    int radiusIndex = -1;
    if (jetRadius == 0.4)
    {
      radiusIndex = 1;
    }
    else if (jetRadius == 0.8)
    {
      radiusIndex = 0;
    }

    //Set all of the fitting values

    //Light Quarks
    //Fitting minimums for the power law fit to the EECs
    double lminl[2][4] = {{0.6, 0.5, 0.3 ,0.25}, {0.325, 0.35, 0.3101, 0.3075}};
    //Fitting maximums for the power law fit to the EECs
    double lmaxl[2][4] = {{0.75, 0.75, 0.75, 0.75}, {0.425, 0.425, 0.455, 0.425}};
    //Fitting minimums for the gaussian fit to the EECs
    double lming[2][4] = {{0.1, 0.08, 0.05, 0.045}, {0.11, 0.07, 0.05, 0.04}};
    //Fitting maximums for the gaussian fit to the EECs
    double lmaxg[2][4] = {{0.275, 0.2, 0.125, 0.1}, {0.23, 0.14, 0.1, 0.09}};

    //Charm Quarks
    double cminl[2][4] = {{0.6, 0.5, 0.3, 0.25}, {0.325, 0.35, 0.3101, 0.3075}};

    double cmaxl[2][4] = {{0.75, 0.75, 0.75, 0.75}, {0.425, 0.425, 0.455, 0.425}};

    double cming[2][4] = {{0.15, 0.08, 0.05, 0.05}, {0.15, 0.08, 0.06, 0.04}};

    double cmaxg[2][4] = {{0.375, 0.2, 0.15, 0.12}, {0.31, 0.17, 0.13, 0.11}};

    //Bottom Quarks
    double bminl[2][4] = {{0.6, 0.5, 0.3, 0.25}, {0.325, 0.35, 0.3101, 0.3075}};
  
    double bmaxl[2][4] = {{0.75, 0.75, 0.75, 0.75}, {0.425, 0.425, 0.455, 0.425}};

    double bming[2][4] = {{0.2, 0.15, 0.1, 0.085}, {0.3, 0.15, 0.12, 0.09}};

    double bmaxg[2][4] = {{0.51, 0.4, 0.275, 0.2}, {0.425, 0.3, 0.22, 0.2}};

    /*
    //Gluons
    double gminl[2][4] = {{0.6, 0.5, 0.3, 0.25}, {0.325, 0.32, 0.3, 0.275}};
    
    double gmaxl[2][4] = {{0.75, 0.75, 0.75, 0.75}, {0.385, 0.385, 0.385, 0.385}};

    double gming[2][4] = {{0.2, 0.15, 0.1, 0.08}, {0.15, 0.11, 0.09, 0.07}};

    double gmaxg[2][4] = {{0.5, 0.4, 0.25, 0.2}, {0.33, 0.25, 0.2, 0.16}};
    */

    //E3C Fitting
    //double e3cminl[4][3] ={0.2}


    int pTHatIndex = -1;
    string histlist[7] = {"eec5", "eec15", "eec 25", "eec 35", "eec45", "eec80", "eec125"};
    

    //Import the files from root locations
    //Sets different bins for different pT ranges
    double bins[5] = {10, 20, 30, 40, 60};
    double badbins[4] = {20,30,40,60};
    //Creates the histograms to be filled with the slope of the power law fits and the peaks of the gaussian fit
    TH1D *slopeHist[4];
    TH1D *peakHist[4];
    TCanvas *c1_final[4];
    //TCanvas *Meson[4];
    TFile *eec[4];
    TCanvas *ratio[4];
    //TCanvas *WTA[4];
    //TCanvas *BMeson[4];

    //Creates the histograms for the slopes and peaks to be recorded.
  for (int i = 0; i<4; i++)
  {
    slopeHist[i] = new TH1D(Form("slopeHist_%s",partSpeciesNames[i].c_str()), Form("Slope Histogram %s;p_{T}^{Jet};Slope",partSpeciesNames[i].c_str()), 3, badbins);
    slopeHist[i]->SetMarkerStyle(marker[i]);
    slopeHist[i]->SetMarkerColor(color[i]);
    slopeHist[i]->SetLineColor(color[i]);
    peakHist[i] = new TH1D(Form("peakHist_%s",partSpeciesNames[i].c_str()), Form("Peak Histogram%s;p_{T}^{Jet};Peak #DeltaR", partSpeciesNames[i].c_str()), 4, bins);
    peakHist[i]->SetMarkerStyle(marker[i]);
    peakHist[i]->SetMarkerColor(color[i]);
    peakHist[i]->SetLineColor(color[i]);
  }
  //Manually sets the canvas up for holding all 4 pT ranges in a 2x2 grid.
  TCanvas *cAll4 = new TCanvas("cAll4","cAll4",1600,1600);
  cAll4->cd()->SetTopMargin(0.15);
  cAll4->Divide(2,2,0.001,0.0);
  //loops through all 4 pT ranges that I looked at
  for (int i=0; i<4; i++)
  {
    //double peaks[4] = {0};

    //Pulls the file from where the user specified
    eec[i] = new TFile(Form("/data/rke_group/millsh1/pp/%s/FJEEC_noEdgeEffects_eCM%d_jetRadius%.1f_pTHat%d_decays%s_MPI%s_n%dk_all.root",outDir.c_str(),ECM,jetRadius,pTHat[i],(decayOn ? "On" : "Off"),(MPIon ? "On" : "Off"),nEvents), "READ");
    //Meson[i] = new TCanvas(Form("Meson%d",i), Form("Meson_%d",i));
    //WTA[i] =  new TCanvas(Form("WTA%d",i), Form("WTA_%d",i));
    //BMeson[i] =  new TCanvas(Form("B Meson%d",i), Form("B Meson%d",i));

    //Creates the individual canvases for each pT range
    c1_final[i] = new TCanvas(Form("c1_final%d",i), Form("Combined Canvas %d", (i+1)));
    c1_final[i]->cd();
    gPad->SetTickx(2);
    gPad->SetTicky(2);
    //Assign the histograms of each type to lh, ch, and bh respectively
    TH1D *lh = (TH1D*)eec[i]->Get("EEC_hist_l");
    lh->SetMarkerStyle(marker[1]);
    lh->SetLineColor(color[1]);
    lh->SetMarkerColor(color[1]);
    TH1D *lhclone = (TH1D*)lh->Clone("lhclone");
    TH1D *ch = (TH1D*)eec[i]->Get("EEC_hist_c");
    ch->SetMarkerStyle(marker[2]);
    ch->SetLineColor(color[2]);
    ch->SetMarkerColor(color[2]);
    TH1D *chclone = (TH1D*)ch->Clone("chclone");
    TH1D *bh = (TH1D*)eec[i]->Get("EEC_hist_b");
    bh->SetMarkerStyle(marker[3]);
    bh->SetLineColor(color[3]);
    bh->SetMarkerColor(color[3]);
    TH1D *bhclone = (TH1D*)bh->Clone("bhclone");
    bhclone->SetMarkerStyle(marker[3]);
    bhclone->SetLineColor(color[3]);
    bhclone->SetMarkerColor(color[3]);
    /*TH1D *bhcloneWTA = (TH1D*)bh->Clone("bhcloneWTA");

    TH1D *bhMeson = (TH1D*)eec[i]->Get("EEC_hist_BMeson_b");
    bhMeson->Add((TH1D*)eec[i]->Get("EEC_hist_BMeson_c"));
    bhMeson->Add((TH1D*)eec[i]->Get("EEC_hist_BMeson_l"));
    bhMeson->Add((TH1D*)eec[i]->Get("EEC_hist_BMeson_g"));
    bhMeson->SetMarkerStyle(marker[3]);
    bhMeson->SetLineColor(color[3]);
    bhMeson->SetMarkerColor(color[3]);

    TH1D *bhWTA = (TH1D*)eec[i]->Get("EEC_hist_WTA_b");
    bhWTA->Add((TH1D*)eec[i]->Get("EEC_hist_WTA_l"));
    bhWTA->Add((TH1D*)eec[i]->Get("EEC_hist_WTA_c"));
    bhWTA->Add((TH1D*)eec[i]->Get("EEC_hist_WTA_g"));
    bhWTA->SetMarkerStyle(marker[3]);
    bhWTA->SetLineColor(color[3]);
    bhWTA->SetMarkerColor(color[3]);
    TH1D *gh = (TH1D*)eec[i]->Get("EEC_hist_g");
    gh->SetMarkerStyle(marker[0]);
    gh->SetLineColor(color[0]);
    gh->SetMarkerColor(color[0]);
    TH1D *ghclone = (TH1D*)gh->Clone("ghclone");
    */
    
    //Fits the perturbative region of the EEC with a power law fit (looks linear when plotted in log space)
    TF1 *ll = new TF1("ll", "[0]*pow(10.0,[1]*log10(x))", lminl[radiusIndex][i], lmaxl[radiusIndex][i]);
    ll->SetNpx(10000);
    ll->SetLineColor(color[1]-2);
    ll->SetLineColor(kTeal-9);
    //ll->SetLineStyle(9);
    ll->SetLineWidth(3);

    //Does the same for the charm initiated jet EEC
    TF1 *cl = new TF1("cl", "[0]*pow(10.0,[1]*log10(x))", cminl[radiusIndex][i], cmaxl[radiusIndex][i]);
    cl->SetNpx(10000);
    cl->SetLineColor(kOrange);
    //cl->SetLineColor(kCyan+2);
    //cl->SetLineStyle(9);
    cl->SetLineWidth(3);

    //Same for Bottom
    TF1 *bl = new TF1("bl", "[0]*pow(10.0,[1]*log10(x))", bminl[radiusIndex][i], bmaxl[radiusIndex][i]);
    bl->SetNpx(10000);
    bl->SetLineColor(870);
    bl->SetLineColor(kViolet-9);
    //bl->SetLineStyle(9);
    bl->SetLineWidth(3);


    //TF1 *gl = new TF1("gl", "[0]*pow(10.0,[1]*log10(x))", gminl[radiusIndex][i], gmaxl[radiusIndex][i]);
    //gl->SetLineColor(897);
    //TF1 *lgauss = new TF1("lgauss","[0]*exp(-0.5*pow(log(x-[1])/[2],2))", 0.03, 0.1);


    //Plots the gaussian fits for the EEC transisition region 
    TF1 *lgauss = new TF1("lgauss","gaus", lming[radiusIndex][i], lmaxg[radiusIndex][i]);
    lgauss->SetNpx(10000);
    lgauss->SetLineColor(color[1] + 2);
    lgauss->SetLineColor(kSpring+7);
    //lgauss->SetLineStyle(2);
    lgauss->SetLineWidth(3);
    TF1 *cgauss = new TF1("cgauss","gaus", cming[radiusIndex][i], cmaxg[radiusIndex][i]);
    cgauss->SetNpx(10000);
    cgauss->SetLineColor(color[2] + 2);
    cgauss->SetLineColor(kPink+1);
    //cgauss->SetLineStyle(2);
    cgauss->SetLineWidth(3);
    TF1 *bgauss = new TF1("bgauss","gaus", bming[radiusIndex][i], bmaxg[radiusIndex][i]);
    bgauss->SetNpx(10000);
    bgauss->SetLineColor(color[3] + 2);
    bgauss->SetLineColor(kAzure+8);
    //bgauss->SetLineStyle(2);
    bgauss->SetLineWidth(3);
    //TF1 *ggauss = new TF1("ggauss","gaus", gming[radiusIndex][i], gmaxg[radiusIndex][i]);
    //ggauss->SetLineColor(color[0] + 2);


    //ll->SetLineColor(kRed);
    //lh->Fit("ll","R","",0.2,0.8);
    //lh->Draw("P");
    //ll->Draw("Psame");

    //TF1 *cl = new TF1("cl");
    //TF1 *bl = new TF1("bl");
    //TF1 *gl = new TF1("gl");


    //Assgin the jetSpec historgrams
    TH1D *lc = (TH1D*)eec[i]->Get("jetSpec_l");
    TH1D *cc = (TH1D*)eec[i]->Get("jetSpec_c");
    TH1D *bc = (TH1D*)eec[i]->Get("jetSpec_b");
    //TH1D *gc = (TH1D*)eec[i]->Get("jetSpec_g");

    //Name the graph and center its axis
    lh->SetTitle("");
    lh->GetXaxis()->CenterTitle(true);
    lh->GetYaxis()->CenterTitle(true);

    //Begins scaling the histograms to all fit on the same plot 
    lh->Scale( 1.0, "width");
    ch->Scale( 1.0, "width");
    bh->Scale( 1.0, "width");
    //bhMeson->Scale( 1.0, "width");
    //gh->Scale( 1.0, "width");


    //Sends the labels for the individual plots into the stratosphere so that there are not extra labels when the 4x4 grid is plotted.
    if (i == 0)
    {
      lh->GetXaxis()->SetLabelOffset(999);
      lh->GetXaxis()->SetTitleOffset(999);
    }
    if (i == 1)
    {
      lh->GetYaxis()->SetLabelOffset(999);
      lh->GetYaxis()->SetTitleOffset(999);
      lh->GetXaxis()->SetLabelOffset(999);
      lh->GetXaxis()->SetTitleOffset(999);
    }
    if (i == 3)
    {
      lh->GetYaxis()->SetLabelOffset(999);
      lh->GetYaxis()->SetTitleOffset(999);
    }



    //TCanvas *c2 = new TCanvas("c2", "Fit for Light");
    //c2->SetLogx();
    //c2->SetLogy();
    //lh->Draw("P");

    //ll->Draw("same");
    //c2->cd();
    //c2->Update();
    //c2->Draw();


    //Set the scaling to line up the graphs so that they are all the same in the free hardon region

    double lNorm = lh->Integral(lh->FindBin(2e-3),lh->FindBin(2e-2));
    double cNorm = ch->Integral(ch->FindBin(2e-3),ch->FindBin(2e-2));
    double bNorm = bh->Integral(bh->FindBin(2e-3),bh->FindBin(2e-2));
    //double gNorm = gh->Integral(bh->FindBin(2e-3),bh->FindBin(2e-2));
    
    lh->Scale(1.0/lNorm);
    ch->Scale(1.0/cNorm);
    bh->Scale(1.0/bNorm);

    //gh->Scale(1.0/gNorm);
    
    //Fit the falling region
    if (i > 0)
    {
    lh->Fit(ll, "R0", "", lminl[radiusIndex][i], lmaxl[radiusIndex][i]);
    slopeHist[1]->Fill(pTmidrange[i], ll->GetParameter(1));
    slopeHist[1]->SetBinError(i, ll->GetParError(1));
    ch->Fit(cl, "R0", "", cminl[radiusIndex][i], cmaxl[radiusIndex][i]);
    slopeHist[2]->Fill(pTmidrange[i], cl->GetParameter(1));
    slopeHist[2]->SetBinError(i, cl->GetParError(1));
    bh->Fit(bl, "R0", "", bminl[radiusIndex][i], bmaxl[radiusIndex][i]);
    slopeHist[3]->Fill(pTmidrange[i], bl->GetParameter(1));
    slopeHist[3]->SetBinError(i, bl->GetParError(1));
    /*gh->Fit(gl, "R0", "", gminl[radiusIndex][i], gmaxl[radiusIndex][i]);
    slopeHist[0]->Fill(pTmidrange[i], gl->GetParameter(1));
    slopeHist[0]->SetBinError(i, gl->GetParError(1));
    */
    }

    //Add the peaks from the gaussian fits to the EEC plot
    lh->Fit(lgauss, "R0", "", lming[radiusIndex][i], lmaxg[radiusIndex][i]);
    peakHist[1]->Fill(pTmidrange[i], lgauss->GetParameter(1));
    peakHist[1]->SetBinError(i+1, lgauss->GetParError(1));
    ch->Fit(cgauss, "R0", "", cming[radiusIndex][i], cmaxg[radiusIndex][i]);
    peakHist[2]->Fill(pTmidrange[i], cgauss->GetParameter(1));
    peakHist[2]->SetBinError(i+1, cgauss->GetParError(1));
    bh->Fit(bgauss, "R0", "", bming[radiusIndex][i], bmaxg[radiusIndex][i]);
    peakHist[3]->Fill(pTmidrange[i], bgauss->GetParameter(1));
    peakHist[3]->SetBinError(i+1, bgauss->GetParError(1));
    /*gh->Fit(ggauss, "R0", "", gming[radiusIndex][i], gmaxg[radiusIndex][i]);
    peakHist[0]->Fill(pTmidrange[i], ggauss->GetParameter(1));
    peakHist[0]->SetBinError(i+1, ggauss->GetParError(1));
    */
    /*lh->Scale(1./lh->Integral());
    ch->Scale(1./ch->Integral());
    bh->Scale(1./bh->Integral());
    */

    //Automatically find the maximum for each plot (used when each of the 4 pT plots were plotted seperately).
    double max=lh->GetMaximum();
    if(ch->GetMaximum() > max)
    {
        max = ch->GetMaximum();
    }
    
    if(bh->GetMaximum() > max)
    {
        max = bh->GetMaximum();
    }
    
    /*if(gh->GetMaximum() > max)
    {
        max = gh->GetMaximum();
    }*/


    //Does the same thing that was done for the Max for the min
    double min=1e10;
    for (int i = 1; i <= lh->GetNbinsX();i++)
    {
        if (lh->GetBinContent(i)<min)
        {
            min = lh->GetBinContent(i);
        }
        if (ch->GetBinContent(i)<min)
        {
            min = ch->GetBinContent(i);
        }
        if (bh->GetBinContent(i)<min)
        {
            min = bh->GetBinContent(i);
        }
        /*if (gh->GetBinContent(i)<min)
        {
            min = gh->GetBinContent(i);
        }*/
    }
    if (min > max)
    {
        min = 1;
    }
    if (min<=0)
    {
        min = 0.001;
    }
    //Sets the user ranges to all be the same to fit when plotted 2x2
    lh->GetYaxis()->SetRangeUser(globalmin,1.25 * globalmax);
    lh->GetXaxis()->SetRangeUser(0.01,0.4);
    std::cout<<min<<"   ;)   "<<max<<"\n";

    //Create new canvas with logrithmic axis
    c1_final[i]->SetLogx();
    c1_final[i]->SetLogy();

    

    //Draw the actual histogram on the canvas
    lh->Draw("p");
    ch->Draw("psame");
    bh->Draw("psame");
    //gh->Draw("psame");
    ll->Draw("same");
    bl->Draw("same");
    cl->Draw("same");
    //gl->Draw("same");
    lgauss->Draw("same");
    cgauss->Draw("same");
    bgauss->Draw("same");
    //ggauss->Draw("same");



    //Add a legend explaining the symbols
    TLegend *leg = new TLegend(0.15, 0.15, 0.85, 0.35);
    leg->SetNColumns(3);
    leg->SetTextSize(0.0275);
    leg->AddEntry(lh, "Light","p");
    leg->AddEntry(ll, "Light Perturbative Region Fit","l");
    leg->AddEntry(lgauss, "Light Transition Region Fit","l");
    leg->AddEntry(ch, "Charm","p");
    leg->AddEntry(cl, "Charm Perturbative Region Fit","l");
    leg->AddEntry(cgauss, "Charm Transition Region Fit","l");
    leg->AddEntry(bh, "Bottom","p");
    leg->AddEntry(bl, "Bottom Perturbative Region Fit","l");
    leg->AddEntry(bgauss, "Bottom Transition Region Fit","l");
    /*leg->AddEntry(gh, "Gluon", "p");
    leg->AddEntry(gl, "Gluon Perturbative Region Fit","l");
    leg->AddEntry(ggauss, "Gluon Transition Region Fit","l");*/
    //leg->Draw();

    //Lable what we did and how we did it
    TPaveText *pT = new TPaveText(0.325, 0.15, 0.675, 0.35, "NDC");
    pT->AddText(Form("PYTHIA-8 pp #sqrt{s} = %d GeV",ECM));
    pT->AddText(Form("%d < p_{T}^{Jet} < %d GeV/c",min_pT[i],max_pT[i]));
    pT->AddText(Form("Anti-k_{T} R = %.1f", jetRadius));
    pT->SetBorderSize(0);
    pT->SetFillStyle(0);
    pT->Draw();

    //Save
    c1_final[i]->SaveAs(Form("/data/rke_group/millsh1/pp/%s/FJEECComp_noEdgeEffects_eCM%d_jetRadius%.1f_pTHat%d_decays%s_MPI%s_n%dk.pdf",outDir.c_str(),ECM,jetRadius,pTHat[i],(decayOn ? "On" : "Off"),(MPIon ? "On" : "Off"),nEvents));

    cAll4->cd(i+1);
    cAll4->cd(i+1)->SetLogy();
    cAll4->cd(i+1)->SetLogx();
    gPad->SetTickx(2);
    gPad->SetTicky(2);

    //Organizes the plots in a 2x2 grid for final presentation
    if(i==1 || i== 3) cAll4->cd(i+1)->SetRightMargin(0.01);
    //if(i==0 || i== 2) cAll4->cd(i+1)->SetLeftMargin(0.06);
    if(i==0 || i==1) cAll4->cd(i+1)->SetTopMargin(0.04);
    if(i==2 || i==3) cAll4->cd(i+1)->SetBottomMargin(0.175);

    //Creates the labels for the final 2x2 plot
    if(i==0 || i==2)
    {
      lh->GetYaxis()->SetTitleSize(0.06);
      lh->GetYaxis()->SetLabelSize(0.045);
      lh->GetYaxis()->SetTitleOffset(0.4);
    }
    if(i==2 || i==3)
    {
      lh->GetXaxis()->SetTitleSize(0.06);
      lh->GetXaxis()->SetLabelSize(0.045);
    }


    //Draws all of the EECs and their various fits with the parameters specified above.
    lh->Draw("p");
    ch->Draw("psame");
    bh->Draw("psame");
    //gh->Draw("psame");
    ll->Draw("same");
    bl->Draw("same");
    cl->Draw("same");
    //gl->Draw("same");
    lgauss->Draw("same");
    cgauss->Draw("same");
    bgauss->Draw("same");

    //TPaveText *pTAll = new TPaveText(0.125, 0.65, 0.475, 0.85, "NDC");
    TPaveText *pTAll = new TPaveText(0.015, 0.9, 0.05, 1.05);
    pTAll->AddText(Form("%d < p_{T}^{Jet} < %d GeV/c",min_pT[i],max_pT[i]));
    pTAll->SetTextSize(0.05);
    pTAll->SetBorderSize(0);
    pTAll->SetFillStyle(0);
    pTAll->Draw();

    //Adds the correct legends to the 2x2 plots
    if(i==0)
    {
      TLegend *lLeg = new TLegend(0.0175, 0.0175, 0.04, 0.07,"","br");
      lLeg->SetTextSize(0.04);
      lLeg->AddEntry(lh, "Light","p");
      lLeg->AddEntry(ch, "Charm","p");
      lLeg->AddEntry(bh, "Bottom","p");
      lLeg->Draw("Same");
    }
    else if(i==1)
    {
      TLegend *cLeg = new TLegend(0.0175, 0.0175, 0.25, 0.07,"","br");
      cLeg->SetTextSize(0.04);
      cLeg->AddEntry(ll, "Light Perturbative Region Fit","l");
      cLeg->AddEntry(cl, "Charm Perturbative Region Fit","l");
      cLeg->AddEntry(bl, "Bottom Perturbative Region Fit","l");
      cLeg->Draw("same");
    }
    else if(i==2){
      TLegend *bLeg = new TLegend(0.0175, 0.0175, 0.25, 0.07,"","br");
      bLeg->SetTextSize(0.04);
      bLeg->AddEntry(lgauss, "Light Transition Region Fit","l");
      bLeg->AddEntry(cgauss, "Charm Transition Region Fit","l");
      bLeg->AddEntry(bgauss, "Bottom Transition Region Fit","l");
      bLeg->Draw("same");
    }

    //The next section is B Meson stuff that we kind of got working that would be cool to continue researching if you feel so inclined

    /*BMeson[i]->cd();
    bhMeson->Draw("P");
    BMeson[i]->SetLogx();
    BMeson[i]->SetLogy();
    gPad->SetTickx(2);
    gPad->SetTicky(2);

    bhclone->Scale(1.0,"width");

    double bMesonNorm = bhMeson->Integral(bhMeson->FindBin(2e-3),bhMeson->FindBin(2e-2));
    bhMeson->Scale(1.0/bMesonNorm);
    double bhcloneNorm = bhclone->Integral(bhclone->FindBin(2e-3),bhclone->FindBin(2e-2));
    bhclone->Scale(1.0/bhcloneNorm);

    bhMeson->Divide(bhclone);
    bhMeson->SetTitle(";#DeltaR;#frac{B Meson}{True b#bar{b}}");
    bhMeson->GetXaxis()->CenterTitle(true);
    bhMeson->GetYaxis()->CenterTitle(true);
    Meson[i]->cd();
    gPad->SetTickx(2);
    gPad->SetTicky(2);
    Meson[i]->SetLogx();
    //Meson[i]->SetLogy();
    bhMeson->Draw("P");
    TPaveText *pTMeson = new TPaveText(0.3, 0.55, 0.65, 0.75, "NDC");
      pTMeson->AddText(Form("PYTHIA-8 pp #sqrt{s} = %d GeV",ECM));
      pTMeson->AddText(Form("%d < p_{T}^{Jet} < %d GeV/c",min_pT[i],max_pT[i]));
      pTMeson->AddText(Form("Anti-k_{T} R = %.1f", jetRadius));
      //pT->AddText(Form("MPI %s", (MPIon ? "On" : "Off")));
      pTMeson->SetBorderSize(0);
      pTMeson->SetFillStyle(0);
      pTMeson->Draw();
    Meson[i]->SaveAs(Form("/data/rke_group/millsh1/pp/%s/FJEECComp_noEdgeEffects_BMeson%d_Analysis_eCM%d_jetRadius%.1f_pTHat%d_decays%s_MPI%s_n%dk.pdf",outDir.c_str(),i,ECM,jetRadius,pTHat[i],(decayOn ? "On" : "Off"),(MPIon ? "On" : "Off"),nEvents));

    bhcloneWTA->Scale(1.0,"width");
    bhWTA->Scale(1.0,"width");
    double bWTANorm = bhWTA->Integral(bhWTA->FindBin(2e-3),bhWTA->FindBin(2e-2));
    bhWTA->Scale(1.0/bWTANorm);
    double bhcloneNormWTA = bhcloneWTA->Integral(bhcloneWTA->FindBin(2e-3),bhcloneWTA->FindBin(2e-2));
    bhcloneWTA->Scale(1.0/bhcloneNormWTA);

    bhWTA->Divide(bhcloneWTA);
    bhWTA->SetTitle(";#DeltaR;#frac{B Meson WTA}{True b#bar{b}}");
    bhWTA->GetXaxis()->CenterTitle(true);
    bhWTA->GetYaxis()->CenterTitle(true);
    WTA[i]->cd();
    gPad->SetTickx(2);
    gPad->SetTicky(2);
    WTA[i]->SetLogx();
    //WTA[i]->SetLogy();
    bhWTA->Draw("P");
    TPaveText *pTWTA = new TPaveText(0.3, 0.55, 0.65, 0.75, "NDC");
      pTWTA->AddText(Form("PYTHIA-8 pp #sqrt{s} = %d GeV",ECM));
      pTWTA->AddText(Form("%d < p_{T}^{Jet} < %d GeV/c",min_pT[i],max_pT[i]));
      pTWTA->AddText(Form("Anti-k_{T} R = %.1f", jetRadius));
      //pT->AddText(Form("MPI %s", (MPIon ? "On" : "Off")));
      pTWTA->SetBorderSize(0);
      pTWTA->SetFillStyle(0);
      pTWTA->Draw();
    WTA[i]->SaveAs(Form("/data/rke_group/millsh1/pp/%s/FJEECComp_noEdgeEffects_BMesonWTA%d_Analysis_eCM%d_jetRadius%.1f_pTHat%d_decays%s_MPI%s_n%dk.pdf",outDir.c_str(),i,ECM,jetRadius,pTHat[i],(decayOn ? "On" : "Off"),(MPIon ? "On" : "Off"),nEvents));
*/
    /*for(int j=0; j<4; j++)
    {
      slopeHist[i]->Fill(j,slopes[j]);
      slopeHist[i]->SetBinError(j,slopesErr[j]);
      slopeHist[i]->GetXaxis()->SetBinLabel(j+1,partSpeciesNames[j].c_str());
      peakHist[i]->Fill(j,peaks[j]);
      peakHist[i]->SetBinError(j,peaksErr[j]);
      peakHist[i]->GetXaxis()->SetBinLabel(j+1,partSpeciesNames[j].c_str());
    }*/

    /*lhsize[r]->Scale( 1.0, "width");
    chsize[r]->Scale( 1.0, "width");
    bhsize[r]->Scale( 1.0, "width");
    ghsize[r]->Scale( 1.0, "width");*/
    //TCanvas *c2 = new TCanvas("c2", "Fit for Light");
    //c2->SetLogx();
    //c2->SetLogy();
    //lh->Draw("P");

    //ll->Draw("same");
    //c2->cd();
    //c2->Update();
    //c2->Draw();
    //Set the scaling to line up the graphs
    //Fit the falling region
    /*lh->Fit(ll, "R0", "", lminl[i], lmaxl[i]);
    slopeHist[1]->Fill(pTmidrange[i], ll->GetParameter(1));
    slopeHist[1]->SetBinError(i, ll->GetParError(1));
    ch->Fit(cl, "R0", "", cminl[i], cmaxl[i]);
    slopeHist[2]->Fill(pTmidrange[i], cl->GetParameter(1));
    slopeHist[2]->SetBinError(i, cl->GetParError(1));
    bh->Fit(bl, "R0", "", bminl[i], bmaxl[i]);
    slopeHist[3]->Fill(pTmidrange[i], bl->GetParameter(1));
    slopeHist[3]->SetBinError(i, bl->GetParError(1));
    gh->Fit(gl, "R0", "", gminl[i], gmaxl[i]);
    slopeHist[0]->Fill(pTmidrange[i], gl->GetParameter(1));
    slopeHist[0]->SetBinError(i, gl->GetParError(1));

    //Fit the peaks
    lh->Fit(lgauss, "R0", "", lming[i], lmaxg[i]);
    peakHist[1]->Fill(pTmidrange[i], lgauss->GetParameter(1));
    peakHist[1]->SetBinError(i, lgauss->GetParError(1));
    ch->Fit(cgauss, "R0", "", cming[i], cmaxg[i]);
    peakHist[2]->Fill(pTmidrange[i], cgauss->GetParameter(1));
    peakHist[2]->SetBinError(i, cgauss->GetParError(1));
    bh->Fit(bgauss, "R0", "", bming[i], bmaxg[i]);
    peakHist[3]->Fill(pTmidrange[i], bgauss->GetParameter(1));
    peakHist[3]->SetBinError(i, bgauss->GetParError(1));
    gh->Fit(ggauss, "R0", "", gming[i], gmaxg[i]);
    peakHist[0]->Fill(pTmidrange[i], ggauss->GetParameter(1));
    peakHist[0]->SetBinError(i, ggauss->GetParError(1));
    */
    /*lh->Scale(1./lh->Integral());
    ch->Scale(1./ch->Integral());
    bh->Scale(1./bh->Integral());
    */
    


    //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
    //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

  }
  //Sets the 2x2 plot as current and adds a title (specifically for STS purposes)
  cAll4->cd();
  TPaveText *cAll4Title = new TPaveText(0.05,0.925,0.95,0.98, "NDC");
  cAll4Title->SetFillStyle(0);
  cAll4Title->SetBorderSize(0);
  cAll4Title->AddText("PYTHIA-8 pp #sqrt{s}=200 GeV Anti-k_{T} R=0.4");
  cAll4Title->Draw("SAME");

  cAll4->SaveAs(Form("/data/rke_group/millsh1/pp/%s/FJEECComp_noEdgeEffects_eCM%d_jetRadius%.1f_allTogether_decays%s_MPI%s_n%dk.pdf",outDir.c_str(),ECM,jetRadius,(decayOn ? "On" : "Off"),(MPIon ? "On" : "Off"),nEvents));
  cAll4->SaveAs(Form("/data/rke_group/millsh1/pp/%s/FJEECComp_noEdgeEffects_eCM%d_jetRadius%.1f_allTogether_decays%s_MPI%s_n%dk.png",outDir.c_str(),ECM,jetRadius,(decayOn ? "On" : "Off"),(MPIon ? "On" : "Off"),nEvents));


    //TFile *s = new TFile(Form("/data/rke_group/millsh1/pp/%s/FJEEC_noEdgeEffects_Slope_eCM%d_pTHat%d_decays%s_n%dk_all.root",outDir.c_str(),ECM,pTHat,(decayOn ? "On" : "Off"),nEvents), "READ");
    TCanvas *slope = new TCanvas("slope", "Slopes");
    gPad->SetTickx(2);
    gPad->SetTicky(2);

    //Automatically sets the mins and maxes for the Slope from the perturbative region power law fits histograms
    double max=slopeHist[1]->GetMaximum();
    if(slopeHist[2]->GetMaximum() > max)
    {
        max = slopeHist[2]->GetMaximum();
    }
    
    if(slopeHist[3]->GetMaximum() > max)
    {
        max = slopeHist[3]->GetMaximum();
    }

    double min=1e10;
    for (int i = 1; i <= slopeHist[1]->GetNbinsX();i++)
    {
        if (slopeHist[1]->GetBinContent(i)< min)
        {
            min = slopeHist[1]->GetBinContent(i);
        }
        if (slopeHist[2]->GetBinContent(i)< min)
        {
            min = slopeHist[2]->GetBinContent(i);
        }
        if (slopeHist[3]->GetBinContent(i)< min)
        {
            min = slopeHist[3]->GetBinContent(i);
        }
    }
    /*if (min > max)
    {
        min = 1;
    }*/
    /*if (min<=0)
    {
        min = 0.001;
    }*/
    slopeHist[1]->GetYaxis()->SetRangeUser(1.25*min+0.25,1.25*max+0.5);
    //slopeHist[0]->GetYaxis()->SetRangeUser(slopeMinRange[radiusIndex],slopeMaxRange[radiusIndex]);
    //slopeHist[0]->Draw("P");
    //slope->cd();
    slopeHist[1]->SetTitle("");
    slopeHist[1]->GetXaxis()->CenterTitle(true);
    slopeHist[1]->GetYaxis()->CenterTitle(true);

    for(int i=1; i<4; i++)
    {
      if(i==1) slopeHist[i]->Draw("p");
      else slopeHist[i]->Draw("psame");
    }
    TPaveText *stext = new TPaveText(0.475, 0.4, 0.875, 0.65, "NDC");
    stext->AddText(Form("PYTHIA-8 pp #sqrt{s} = %d GeV",ECM));
    stext->AddText(Form("Anti-k_{T} R = %.1f", jetRadius));
    stext->AddText("Power Law Fit to EEC");
    stext->AddText("Perturbative Region");
    stext->SetBorderSize(0);
    stext->SetFillStyle(0);
    stext->Draw();
    TLegend *slopeleg = new TLegend(0.625, 0.6625, 0.875, 0.8875);
    slopeleg->AddEntry(slopeHist[1], "Light Flavor","p");
    slopeleg->AddEntry(slopeHist[2], "Charm","p");
    slopeleg->AddEntry(slopeHist[3], "Bottom","p");
    //slopeleg->AddEntry(slopeHist[0], "Gluon","p");
    slopeleg->Draw();

    slope->SaveAs(Form("/data/rke_group/millsh1/pp/%s/FJEECComp_noEdgeEffects_Slope_eCM%d_jetRadius%.1f_pTHat_all_decays%s_MPI%s_n%dk.pdf",outDir.c_str(),ECM,jetRadius,(decayOn ? "On" : "Off"),(MPIon ? "On" : "Off"),nEvents));

    //Sets the minimum and maxiumum viewing range automatically to include all of the data points.
    max=peakHist[0]->GetMaximum();
    if(peakHist[1]->GetMaximum() > max)
    {
        max = peakHist[1]->GetMaximum();
    }
    
    if(peakHist[2]->GetMaximum() > max)
    {
        max = peakHist[2]->GetMaximum();
    }
    
    if(peakHist[3]->GetMaximum() > max)
    {
        max = peakHist[3]->GetMaximum();
    }

    min=1e10;
    for (int i = 1; i <= peakHist[0]->GetNbinsX();i++)
    {
        if (peakHist[0]->GetBinContent(i)< min)
        {
            min = peakHist[0]->GetBinContent(i);
        }
        if (peakHist[1]->GetBinContent(i)< min)
        {
            min = peakHist[1]->GetBinContent(i);
        }
        if (peakHist[2]->GetBinContent(i)< min)
        {
            min = peakHist[2]->GetBinContent(i);
        }
        if (peakHist[3]->GetBinContent(i)< min)
        {
            min = peakHist[3]->GetBinContent(i);
        }
    }
    std:cout<<"min "<<min<<"\n";
    std::cout<<"max "<<max<<"\n";
    peakHist[0]->GetYaxis()->SetRangeUser(0, 1.25*max);
    TCanvas *peak = new TCanvas("peak", "Peaks");
    gPad->SetTickx(2);
    gPad->SetTicky(2);
    //peakHist[0]->GetYaxis()->SetRangeUser(peakMinRange[radiusIndex],peakMaxRange[radiusIndex]);
    peakHist[0]->SetTitle("");
    peakHist[0]->GetXaxis()->CenterTitle(true);
    peakHist[0]->GetYaxis()->CenterTitle(true);
    //Draws all of the peak histograms for each flavor
    for(int i=0; i<4; i++)
    {
      if(i==0) peakHist[i]->Draw("p");
      else peakHist[i]->Draw("psame");
    }
    //Plots the Text for the Gaussian fits peaks histogram
    TPaveText *ptext = new TPaveText(0.325, 0.675, 0.625, 0.875, "NDC");
    ptext->AddText(Form("PYTHIA-8 pp #sqrt{s} = %d GeV",ECM));
    ptext->AddText(Form("Anti-k_{T} R = %.1f", jetRadius));
    ptext->AddText("Gaussian Fit to EEC Peak");
    ptext->SetBorderSize(0);
    ptext->SetFillStyle(0);
    ptext->Draw();
    //Plots the legend for the gaussian fit peaks histogram
    TLegend *pleg = new TLegend(0.625, 0.6625, 0.875, 0.8875);
    pleg->AddEntry(peakHist[1], "Light Flavor","p");
    pleg->AddEntry(peakHist[2], "Charm","p");
    pleg->AddEntry(peakHist[3], "Bottom","p");
    //pleg->AddEntry(peakHist[0], "Gluon", "p");
    pleg->Draw();
    peak->SaveAs(Form("/data/rke_group/millsh1/pp/%s/FJEECComp_noEdgeEffects_Peak_eCM%d_jetRadius%.1f_pTHat_all_decays%s_MPI%s_n%dk.pdf",outDir.c_str(),ECM,jetRadius,(decayOn ? "On" : "Off"),(MPIon ? "On" : "Off"),nEvents));
    /*
    //TCanvas *s1 = new TCanvas ("s1", "Slope Comparison Histogram");
    //slopeHist->Draw("P");
    TH1D *peakHist = new TH1D("peakHist", "Peak Histogram;Initiating Parton;Slope", 4, -0.5, 3.5);
    for(int i=0; i<4; i++)
    {
      peakHist->Fill(i,peaks[i]);
      peakHist->GetXaxis()->SetBinLabel(i+1,partSpeciesNames[i].c_str());
    }*/
    //TCanvas *p1 = new TCanvas ("p1", "Peak Comparison Histogram");
    //peakHist->Draw("P");
}

