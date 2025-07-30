double langau_expo_fun(double *x, double *par) {

    // --- Langau Part (your original langaufun logic) ---
    double invsq2pi = 0.3989422804014;
    double mpshift  = -0.22278298;
    double np = 100.0;
    double sc = 5.0;
    double xx, mpc, fland;
    double sum = 0.0;
    double xlow, xupp, step;

    mpc = par[1] - mpshift * par[0];
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    step = (xupp - xlow) / np;

    int len = np / 2;
    for(int i = 0; i <= len; i++) {
        xx = xlow + (i - 0.5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);

        xx = xupp - (i - 0.5) * step;
        fland = TMath::Landau(xx, mpc, par[0]) / par[0];
        sum += fland * TMath::Gaus(x[0], xx, par[3]);
    }
    double langau_val = (par[2] * step * sum * invsq2pi / par[3]);

    // --- Exponential Part ---
    double expo_val = par[4] * TMath::Exp(par[5] * x[0]);

    // --- Return Combined Function ---
    return langau_val + expo_val;
}

TF1 *langau_expo_fit(TH1F *his, double *fitrange, double *startvalues, double *parlimitslo, double *parlimitshi, double *fitparams, double *fiterrors, double *ChiSqr, int *NDF)
{
    char FunName[100];
    sprintf(FunName, "Fitfcn_%s", his->GetName());

    TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
    if (ffitold) delete ffitold;

    TF1 *ffit = new TF1(FunName, langau_expo_fun, fitrange[0], fitrange[1], 6);
    ffit->SetParameters(startvalues);
    ffit->SetParNames("Width", "MP", "Area", "GSigma", "E_Amp", "E_Slope");

    for (int i = 0; i < 6; i++) {
        ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
    }

    his->Fit(FunName, "RB0");

    ffit->GetParameters(fitparams);
    for (int i = 0; i < 6; i++) {
        fiterrors[i] = ffit->GetParError(i);
    }
    ChiSqr[0] = ffit->GetChisquare();
    NDF[0] = ffit->GetNDF();

    return (ffit);
}