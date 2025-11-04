#ifndef PLOT_CHECK_H
#define PLOT_CHECK_H

#include <string>
#include <vector>
#include <utility>
#include <functional>
#include <cmath>
#include <algorithm>
#include <map>

// ROOT headers
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TColor.h>
#include <TArrow.h>
#include <TMarker.h>
#include <TList.h>
#include <TObjString.h>
#include <TProfile2D.h>
#include <TSystem.h>
#include <TCollection.h>

// My function headers
#include "Messenger/Messenger.h"
#include "plot_histogram/plot_histogram.h"
#include "plot_ExperimentData/plot_ExperimentData.h"
#include "plot_BeamTest/plot_BeamTest.h"

// Other function headers
#include "tools/json.hpp"
#include "tools/cxxopts.hpp"

struct DatasetConfig {
    std::string name; // plot_title
    std::string suffix; // _gap

    std::string expFileBase; // file base name
    double expScale; // experiment scaling factor

    std::string simFile; // simulation file base name
    double simScale;
};


class plot_Check {
public:
    
    plot_Check();
    ~plot_Check();

    TH1D* scalingHistogram(TH1D* hist, double scale_factor, std::string title);
    TH1D* scalingClusterSize(TH1D* hist);

    void run_Check();
    void run_plotCheck();

private:
    std::map<std::string, TFile*> m_fileCache;
    TFile* openFile(const std::string& fileName);
    void cleanupResources();

    void setHistStyle(TH1D* hist, int color, int marker, float alpha = 0.2f);
    TH1D* getScaledHist(const std::string& file, const std::string& histName, double scaleFactor, const std::string& title, bool isMerged = false);

    void plotAllCharge(TCanvas* c, TLegend* l, const DatasetConfig& ds);

    void plotClusterSize(TCanvas* c, TLegend* l, const DatasetConfig& ds);

    void plotComparison(TCanvas* c, TLegend* l, const DatasetConfig& ds,
                        const std::string& expHistName, const std::string& simHistName,
                        const std::string& expLegend, const std::string& simLegend,
                        const std::string& outName, const std::string& plotTitle,
                        int rebin = 1, double xMin = 0, double xMax = 4000,
                        bool normalizedToMax = true);

    void plotSeedChargeByCS(TCanvas* c, TLegend* l, const DatasetConfig& ds);

    void plotDepositedCharge(TCanvas* c, TLegend* l);
    
};

#endif // PLOT_CHECK_H
