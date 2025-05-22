#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"

// double gaussian
double twogausfun(double *x, double *par) {
    return par[0] * TMath::Gaus(x[0], par[1], par[2]) + par[3] * TMath::Gaus(x[0], par[4], par[5]);
}

TF1 *twogausfit(TH1F *hist, double *fitrange, double *startvalues, double *parlimitslo,
                double *parlimitshi, double *fitparames, double *fiterrors, double *ChiSqr, int *NDF) {
    int i;
    char FunName[100];

    sprintf(FunName, "Fitfcn_%s", hist->GetName());

    TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
    if (ffitold) delete ffitold;

    // two gaussian
    TF1 *ffit = new TF1(FunName, twogausfun, fitrange[0], fitrange[1], 6);
    ffit->SetParameters(startvalues);
    ffit->SetParNames("Constant1", "Mean1", "Sigma1", "Constant2", "Mean2", "Sigma2");

    for(i=0; i<6; i++) {
        ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
    }

    hist->Fit(FunName, "RB0");

    ffit->GetParameters(fitparames);
    for(i=0; i<6; i++) {
        fiterrors[i] = ffit->GetParError(i);
    }
    ChiSqr[0] = ffit->GetChisquare();   
    NDF[0] = ffit->GetNDF();

    return (ffit);
}