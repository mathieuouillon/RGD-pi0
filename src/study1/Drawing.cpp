#include <study1/Drawing.hpp>

namespace study1 {

Drawing::Drawing(Histograms& histograms, const toml::parse_result& config)
    : m_histograms(histograms), m_config(config) {
}

auto Drawing::fit_time_photon(const std::shared_ptr<TH1D>& hist1D_time_cut_fit) const -> void {
    auto canvas = Plotting::make_canvas();

    hist1D_time_cut_fit->Draw("hist");
    hist1D_time_cut_fit->SetLineColor(Plotting::Color::kBlue);

    TF1 fit("fit_time", "gausn", -1, 1);
    hist1D_time_cut_fit->Fit(&fit, "R");
    fit.SetLineColor(Plotting::Color::kRed);
    fit.SetLineWidth(2);
    fit.Draw("same");

    hist1D_time_cut_fit->GetXaxis()->SetTitle("#Deltat [ns]");
    hist1D_time_cut_fit->GetXaxis()->SetTitleOffset(0.8f);
    hist1D_time_cut_fit->GetXaxis()->SetTitleSize(0.05f);

    Plotting::save_canvas(canvas, m_path_photon, "time_fit");
}

auto Drawing::fit_opening_angle(const std::shared_ptr<TH1D>& hist1D_opening_angle_cut_fit) const -> void {
    auto canvas = Plotting::make_canvas();

    hist1D_opening_angle_cut_fit->Draw("hist");
    hist1D_opening_angle_cut_fit->SetLineColor(Plotting::Color::kBlue);

    TF1 fit("fit_opening_angle", "gausn", 20, 40);
    hist1D_opening_angle_cut_fit->Fit(&fit, "R");
    fit.SetLineColor(Plotting::Color::kRed);
    fit.SetLineWidth(2);
    fit.Draw("same");

    hist1D_opening_angle_cut_fit->GetXaxis()->SetTitle("#theta_{e,#gamma} [deg.]");
    hist1D_opening_angle_cut_fit->GetXaxis()->SetTitleOffset(0.8f);
    hist1D_opening_angle_cut_fit->GetXaxis()->SetTitleSize(0.05f);

    Plotting::save_canvas(canvas, m_path_photon, "opening_angle_fit");
}

auto Drawing::draw_electron() -> void {

    // Cuts ---------------------------------------------------------------------------------------
    const double vz_min_e = m_config["electron"]["vz_min"].value_or(NaN);
    const double vz_max_e = m_config["electron"]["vz_max"].value_or(NaN);
    const double chi2_min_e = m_config["electron"]["chi2_min"].value_or(NaN);
    const double chi2_max_e = m_config["electron"]["chi2_max"].value_or(NaN);
    const double p_min = m_config["electron"]["p_min"].value_or(NaN);
    //---------------------------------------------------------------------------------------------

    // Merge 1D histograms ------------------------------------------------------------------------
    const std::shared_ptr<TH1D> hist1D_p = m_histograms.electron.hist1D_p->Merge();
    const std::shared_ptr<TH1D> hist1D_chi2 = m_histograms.electron.hist1D_chi2->Merge();
    const std::shared_ptr<TH1D> hist1D_phi = m_histograms.electron.hist1D_phi->Merge();
    const std::shared_ptr<TH1D> hist1D_theta = m_histograms.electron.hist1D_theta->Merge();
    const std::shared_ptr<TH1D> hist1D_vz = m_histograms.electron.hist1D_vz->Merge();

    const std::shared_ptr<TH1D> hist1D_p_cut = m_histograms.electron.hist1D_p_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_chi2_cut = m_histograms.electron.hist1D_chi2_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_phi_cut = m_histograms.electron.hist1D_phi_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_theta_cut = m_histograms.electron.hist1D_theta_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_vz_cut = m_histograms.electron.hist1D_vz_cut->Merge();
    // --------------------------------------------------------------------------------------------

    Plotting::draw_hist1D(hist1D_p, hist1D_p_cut, m_path_electron, {.cuts = {p_min}, .label = "p_{e} [GeV/c]"});
    Plotting::draw_hist1D(hist1D_chi2, hist1D_chi2_cut, m_path_electron, {.cuts = {chi2_min_e, chi2_max_e}, .label = "chi2pid_{e}"});
    Plotting::draw_hist1D(hist1D_phi, hist1D_phi_cut, m_path_electron, {.label = "#phi_{e} [deg.]"});
    Plotting::draw_hist1D(hist1D_theta, hist1D_theta_cut, m_path_electron, {.label = "#theta_{e} [deg.]"});
    Plotting::draw_hist1D(hist1D_vz, hist1D_vz_cut, m_path_electron, {.cuts = {vz_min_e, vz_max_e}, .log_y = true, .label = "Vz_{e} [cm]"});
}

auto Drawing::draw_photon() -> void {
    // Cuts values --------------------------------------------------------------------------------
    const double E_min = m_config["photon"]["E_min"].value_or(NaN);
    const double time_min = m_config["photon"]["time_min"].value_or(NaN);
    const double time_max = m_config["photon"]["time_max"].value_or(NaN);
    const double opening_angle_min = m_config["photon"]["opening_angle_min"].value_or(NaN);
    //---------------------------------------------------------------------------------------------

    // Merge 1D histograms ------------------------------------------------------------------------
    const std::shared_ptr<TH1D> hist1D_E = m_histograms.photon.hist1D_E->Merge();
    const std::shared_ptr<TH1D> hist1D_beta = m_histograms.photon.hist1D_beta->Merge();
    const std::shared_ptr<TH1D> hist1D_time = m_histograms.photon.hist1D_time->Merge();
    const std::shared_ptr<TH1D> hist1D_vz = m_histograms.photon.hist1D_vz->Merge();
    const std::shared_ptr<TH1D> hist1D_opening_angle = m_histograms.photon.hist1D_opening_angle->Merge();
    const std::shared_ptr<TH1D> hist1D_delta_vz = m_histograms.photon.hist1D_delta_vz->Merge();

    const std::shared_ptr<TH2D> hist2D_phi_vs_theta = m_histograms.photon.hist2D_phi_vs_theta->Merge();
    const std::shared_ptr<TH2D> hist2D_E_vs_vz = m_histograms.photon.hist2D_E_vs_vz->Merge();

    const std::shared_ptr<TH1D> hist1D_E_cut = m_histograms.photon.hist1D_E_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_beta_cut = m_histograms.photon.hist1D_beta_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_time_cut = m_histograms.photon.hist1D_time_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_vz_cut = m_histograms.photon.hist1D_vz_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_opening_angle_cut = m_histograms.photon.hist1D_opening_angle_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_delta_vz_cut = m_histograms.photon.hist1D_delta_vz_cut->Merge();

    const std::shared_ptr<TH2D> hist2D_phi_vs_theta_cut = m_histograms.photon.hist2D_phi_vs_theta_cut->Merge();
    const std::shared_ptr<TH2D> hist2D_E_vs_vz_cut = m_histograms.photon.hist2D_E_vs_vz_cut->Merge();

    const std::shared_ptr<TH1D> hist1D_opening_angle_cut_fit = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(hist1D_opening_angle_cut->Clone()));
    const std::shared_ptr<TH1D> hist1D_time_cut_fit = std::shared_ptr<TH1D>(dynamic_cast<TH1D*>(hist1D_time_cut->Clone()));
    // --------------------------------------------------------------------------------------------

    Plotting::draw_hist1D(hist1D_E, hist1D_E_cut, m_path_photon, {.cuts = {E_min}, .label = "E_{#gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_beta, hist1D_beta_cut, m_path_photon, {.label = "#beta_{#gamma}"});
    Plotting::draw_hist1D(hist1D_time, hist1D_time_cut, m_path_photon, {.cuts = {time_min, time_max}, .label = "#Deltat [ns]"});
    Plotting::draw_hist1D(hist1D_vz, hist1D_vz_cut, m_path_photon, {.label = "Vz_{#gamma} [cm]"});
    Plotting::draw_hist1D(hist1D_opening_angle, hist1D_opening_angle_cut, m_path_photon, {.cuts = {opening_angle_min}, .label = "#theta_{e,#gamma} [deg.]"});

    Plotting::draw_hist2D(hist2D_phi_vs_theta, m_path_photon, {.label_x = "#phi_{#gamma} [deg.]", .label_y = "#theta_{#gamma} [deg.]"});
    Plotting::draw_hist2D(hist2D_phi_vs_theta_cut, m_path_photon, {.label_x = "#phi_{#gamma} [deg.]", .label_y = "#theta_{#gamma} [deg.]"});
    Plotting::draw_hist2D(hist2D_E_vs_vz, m_path_photon, {.label_x = "E_{#gamma} [GeV]", .label_y = "Vz_{#gamma} [cm]"});
    Plotting::draw_hist2D(hist2D_E_vs_vz_cut, m_path_photon, {.label_x = "E_{#gamma} [GeV]", .label_y = "Vz_{#gamma} [cm]"});

    fit_time_photon(hist1D_time_cut_fit);
    fit_opening_angle(hist1D_opening_angle_cut_fit);
}

auto Drawing::draw_event() -> void {
    // Cuts ---------------------------------------------------------------------------------------
    const double Q2_min = m_config["event"]["Q2_min"].value_or(NaN);
    const double W_min = m_config["event"]["W_min"].value_or(NaN);
    const double zh_min = m_config["event"]["zh_min"].value_or(NaN);
    const double zh_max = m_config["event"]["zh_max"].value_or(NaN);
    const double y_max = m_config["event"]["y_max"].value_or(NaN);
    //---------------------------------------------------------------------------------------------

    // Merge 1D histograms ------------------------------------------------------------------------
    const std::shared_ptr<TH1D> hist1D_W = m_histograms.event.hist1D_W->Merge();
    const std::shared_ptr<TH1D> hist1D_Q2 = m_histograms.event.hist1D_Q2->Merge();
    const std::shared_ptr<TH1D> hist1D_nu = m_histograms.event.hist1D_nu->Merge();
    const std::shared_ptr<TH1D> hist1D_zh = m_histograms.event.hist1D_zh->Merge();
    const std::shared_ptr<TH1D> hist1D_y = m_histograms.event.hist1D_y->Merge();
    const std::shared_ptr<TH1D> hist1D_invariant_mass = m_histograms.event.hist1D_invariant_mass->Merge();

    const std::shared_ptr<TH1D> hist1D_W_cut = m_histograms.event.hist1D_W_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_Q2_cut = m_histograms.event.hist1D_Q2_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_nu_cut = m_histograms.event.hist1D_nu_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_zh_cut = m_histograms.event.hist1D_zh_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_y_cut = m_histograms.event.hist1D_y_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_invariant_mass_cut = m_histograms.event.hist1D_invariant_mass_cut->Merge();

    const std::shared_ptr<TH1D> hist1D_invariant_mass_mix_cut = m_histograms.event_mix.hist1D_invariant_mass_cut->Merge();
    // --------------------------------------------------------------------------------------------

    Plotting::draw_hist1D(hist1D_W, hist1D_W_cut, m_path_event, {.cuts = {W_min}, .label = "W [GeV]"});
    Plotting::draw_hist1D(hist1D_Q2, hist1D_Q2_cut, m_path_event, {.cuts = {Q2_min}, .label = "Q2 [GeV^2]"});
    Plotting::draw_hist1D(hist1D_nu, hist1D_nu_cut, m_path_event, {.label = "#nu [GeV]"});
    Plotting::draw_hist1D(hist1D_zh, hist1D_zh_cut, m_path_event, {.cuts = {zh_min, zh_max}, .label = "zh"});
    Plotting::draw_hist1D(hist1D_y, hist1D_y_cut, m_path_event, {.cuts = {y_max}, .label = "y"});
    Plotting::draw_hist1D(hist1D_invariant_mass, hist1D_invariant_mass_cut, m_path_event, {.label = "Invariant mass [GeV]"});
    
    
    
    double sidebandData = hist1D_invariant_mass_cut->Integral(hist1D_invariant_mass_cut->FindBin(0.19), hist1D_invariant_mass_cut->FindBin(0.3));
    double sidebandMixed = hist1D_invariant_mass_mix_cut->Integral(hist1D_invariant_mass_mix_cut->FindBin(0.19), hist1D_invariant_mass_mix_cut->FindBin(0.3));
    double scaleFactor = sidebandData / sidebandMixed;
    fmt::println("scale: {}", scaleFactor);
    
    Plotting::draw_hist1D(hist1D_invariant_mass_cut, hist1D_invariant_mass_mix_cut, m_path_event, {.file_name = "invariant_mass_mix", .scale2 = scaleFactor, .label = "Invariant mass [GeV]"});
}

}  // namespace study1