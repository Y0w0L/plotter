import ROOT

def perform_root_fit(file_path):
    # Open the ROOT file in read mode
    file = ROOT.TFile.Open(file_path, "READ")
    if not file or file.IsZombie():
        print(f"Error: Could not open file {file_path}")
        return None

    # Enable fit parameters display in the legend/stat box
    # 1111 means: probabilities, errors, mean, and sigma are shown
    ROOT.gStyle.SetOptFit(1111)

    target_hists = ["trackAngleX", "trackAngleY"]
    results = {}

    for hist_name in target_hists:
        hist = file.Get(f"Tracking4D/{hist_name}")
        
        if not hist:
            print(f"Error: Histogram '{hist_name}' not found.")
            continue

        canvas = ROOT.TCanvas(f"c_{hist_name}", f"Fit - {hist_name}", 800, 600)
        
        # --- Range Adjustment ---
        # Narrow down the display range to focus on the peak
        # Adjust these values based on your specific distribution
        hist.GetXaxis().SetRangeUser(-0.002, 0.002)

        # --- Fitting ---
        # Option "S" to return TFitResultPtr
        fit_result = hist.Fit("gaus", "S")

        if fit_result.IsValid():
            sigma_val = fit_result.Parameter(2)
            results[hist_name] = sigma_val
            
            # Print to console for verification
            print(f"[{hist_name}] Sigma: {sigma_val:.4e}")

        # --- Visualization ---
        hist.Draw()
        
        # Update the canvas to ensure the stat box/legend is generated
        canvas.Modified()
        canvas.Update()

        # Save as image
        canvas.SaveAs(f"fit_result_{hist_name}.pdf")
        
        # Prevent the canvas from being garbage collected
        ROOT.SetOwnership(canvas, False)

    return results

# --- Execution ---
file_name = "/home/towa/alice3/hist/sps202404/sps202404_15_gap_4V_SeedThd500e_NeighborThd500e_2.root"
sigmas = perform_root_fit(file_name)

if sigmas:
    print("\nSummary of Sigmas:")
    for key, val in sigmas.items():
        print(f"{key}: {val:.4e}")