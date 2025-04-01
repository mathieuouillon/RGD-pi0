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
    //---------------------------------------------------------------------------------------------

    // Merge 1D histograms ------------------------------------------------------------------------
    const std::shared_ptr<TH1D> hist1D_W = m_histograms.event.hist1D_W->Merge();
    const std::shared_ptr<TH1D> hist1D_Q2 = m_histograms.event.hist1D_Q2->Merge();
    const std::shared_ptr<TH1D> hist1D_nu = m_histograms.event.hist1D_nu->Merge();
    const std::shared_ptr<TH1D> hist1D_zh = m_histograms.event.hist1D_zh->Merge();
    const std::shared_ptr<TH1D> hist1D_invariant_mass = m_histograms.event.hist1D_invariant_mass->Merge();

    const std::shared_ptr<TH1D> hist1D_W_cut = m_histograms.event.hist1D_W_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_Q2_cut = m_histograms.event.hist1D_Q2_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_nu_cut = m_histograms.event.hist1D_nu_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_zh_cut = m_histograms.event.hist1D_zh_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_invariant_mass_cut = m_histograms.event.hist1D_invariant_mass_cut->Merge();
    // --------------------------------------------------------------------------------------------

    Plotting::draw_hist1D(hist1D_W, hist1D_W_cut, m_path_event, {.cuts = {W_min}, .label = "W [GeV]"});
    Plotting::draw_hist1D(hist1D_Q2, hist1D_Q2_cut, m_path_event, {.cuts = {Q2_min}, .label = "Q2 [GeV^2]"});
    Plotting::draw_hist1D(hist1D_nu, hist1D_nu_cut, m_path_event, {.label = "#nu [GeV]"});
    Plotting::draw_hist1D(hist1D_zh, hist1D_zh_cut, m_path_event, {.cuts = {zh_min, zh_max}, .label = "zh"});
    Plotting::draw_hist1D(hist1D_invariant_mass, hist1D_invariant_mass_cut, m_path_event, {.label = "Invariant mass [GeV]"});

    // Merge 2D histograms ------------------------------------------------------------------------
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_Q2_cut = m_histograms.event.hist2D_invariant_mass_cut_vs_Q2_cut->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_zh_cut = m_histograms.event.hist2D_invariant_mass_cut_vs_zh_cut->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_nu_cut = m_histograms.event.hist2D_invariant_mass_cut_vs_nu_cut->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_W_cut = m_histograms.event.hist2D_invariant_mass_cut_vs_W_cut->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_y = m_histograms.event.hist2D_invariant_mass_cut_vs_y->Merge();

    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_opening_angle_least = m_histograms.event.hist2D_invariant_mass_cut_vs_opening_angle_least->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_opening_angle_high = m_histograms.event.hist2D_invariant_mass_cut_vs_opening_angle_high->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_E_least = m_histograms.event.hist2D_invariant_mass_cut_vs_E_least->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_E_high = m_histograms.event.hist2D_invariant_mass_cut_vs_E_high->Merge();

    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_p_e = m_histograms.event.hist2D_invariant_mass_cut_vs_p_e->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_theta_e = m_histograms.event.hist2D_invariant_mass_cut_vs_theta_e->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_phi_e = m_histograms.event.hist2D_invariant_mass_cut_vs_phi_e->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_chi2_e = m_histograms.event.hist2D_invariant_mass_cut_vs_chi2_e->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_vz_e = m_histograms.event.hist2D_invariant_mass_cut_vs_vz_e->Merge();

    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_beta_g_least = m_histograms.event.hist2D_invariant_mass_cut_vs_beta_g_least->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_beta_g_high = m_histograms.event.hist2D_invariant_mass_cut_vs_beta_g_high->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_time_g_least = m_histograms.event.hist2D_invariant_mass_cut_vs_time_g_least->Merge();
    const std::shared_ptr<TH2D> hist2D_invariant_mass_cut_vs_time_g_high = m_histograms.event.hist2D_invariant_mass_cut_vs_time_g_high->Merge();
    // --------------------------------------------------------------------------------------------

    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_Q2_cut, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "Q2 [GeV^2]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_zh_cut, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "zh"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_nu_cut, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "#nu [GeV]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_W_cut, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "W [GeV]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_y, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "y"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_opening_angle_least, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "#theta_{e,#gamma} [deg.]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_opening_angle_high, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "#theta_{e,#gamma} [deg.]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_E_least, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "E_{#gamma} [GeV]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_E_high, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "E_{#gamma} [GeV]"});

    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_p_e, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "p_{e} [GeV/c]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_theta_e, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "#theta_{e} [deg.]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_phi_e, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "#phi_{e} [deg.]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_chi2_e, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "chi2pid_{e}"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_vz_e, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "Vz_{e} [cm]"});

    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_beta_g_least, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "#beta_{#gamma}"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_beta_g_high, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "#beta_{#gamma}"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_time_g_least, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "#Deltat [ns]"});
    Plotting::draw_hist2D(hist2D_invariant_mass_cut_vs_time_g_high, m_path_event, {.label_x = "Invariant mass [GeV]", .label_y = "#Deltat [ns]"});


    // Draw the photon information
    draw_photon_info(m_histograms.event.photon_info, m_path_photon_info);
    draw_photon_info(m_histograms.event.photon_info_off, m_path_photon_info_off);


}


auto Drawing::draw_photon_info(const PhotonInfo& info, const std::string& path) -> void {

    // Photon information
    const std::shared_ptr<TH1D> hist1D_E_least = info.hist1D_E_least->Merge();
    const std::shared_ptr<TH1D> hist1D_E_high = info.hist1D_E_high->Merge();
    const std::shared_ptr<TH1D> hist1D_beta_least = info.hist1D_beta_least->Merge();
    const std::shared_ptr<TH1D> hist1D_beta_high = info.hist1D_beta_high->Merge();
    const std::shared_ptr<TH1D> hist1D_time_least = info.hist1D_time_least->Merge();
    const std::shared_ptr<TH1D> hist1D_time_high = info.hist1D_time_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Epcal_high = info.hist1D_Epcal_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Epcal_least = info.hist1D_Epcal_least->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_in_high = info.hist1D_Ecalo_in_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_in_least = info.hist1D_Ecalo_in_least->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_out_high = info.hist1D_Ecalo_out_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_out_least = info.hist1D_Ecalo_out_least->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m2u_high = info.hist1D_Ecalo_m2u_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m2u_least = info.hist1D_Ecalo_m2u_least->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m2v_high = info.hist1D_Ecalo_m2v_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m2v_least = info.hist1D_Ecalo_m2v_least->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m2w_high = info.hist1D_Ecalo_m2w_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m2w_least = info.hist1D_Ecalo_m2w_least->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m3u_high = info.hist1D_Ecalo_m3u_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m3u_least = info.hist1D_Ecalo_m3u_least->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m3v_high = info.hist1D_Ecalo_m3v_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m3v_least = info.hist1D_Ecalo_m3v_least->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m3w_high = info.hist1D_Ecalo_m3w_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_m3w_least = info.hist1D_Ecalo_m3w_least->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_tot_high = info.hist1D_Ecalo_tot_high->Merge();
    const std::shared_ptr<TH1D> hist1D_Ecalo_tot_least = info.hist1D_Ecalo_tot_least->Merge();

    Plotting::draw_hist1D(hist1D_E_least, path, {.label = "E_{#gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_E_high, path, {.label = "E_{#gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_beta_least, path, {.label = "#beta_{#gamma}"});
    Plotting::draw_hist1D(hist1D_beta_high, path, {.label = "#beta_{#gamma}"});
    Plotting::draw_hist1D(hist1D_time_least, path, {.label = "#Deltat [ns]"});
    Plotting::draw_hist1D(hist1D_time_high, path, {.label = "#Deltat [ns]"});
    Plotting::draw_hist1D(hist1D_Epcal_least, path, {.label = "E_{PCAL, #gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_Epcal_high, path, {.label = "E_{PCAL, #gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_Ecalo_in_least, path, {.label = "E_{IN, #gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_Ecalo_in_high, path, {.label = "E_{IN, #gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_Ecalo_out_least, path, {.label = "E_{OUT, #gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_Ecalo_out_high, path, {.label = "E_{OUT, #gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_Ecalo_m2u_least, path, {.label = "m2u_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m2u_high, path, {.label = "m2u_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m2v_least, path, {.label = "m2v_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m2v_high, path, {.label = "m2v_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m2w_least, path, {.label = "m2w_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m2w_high, path, {.label = "m2w_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m3u_least, path, {.label = "m3u_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m3u_high, path, {.label = "m3u_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m3v_least, path, {.label = "m3v_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m3v_high, path, {.label = "m3v_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m3w_least, path, {.label = "m3w_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_m3w_high, path, {.label = "m3w_{#gamma}"});
    Plotting::draw_hist1D(hist1D_Ecalo_tot_least, path, {.label = "E_{CAL,tot, #gamma} [GeV]"});
    Plotting::draw_hist1D(hist1D_Ecalo_tot_high, path, {.label = "E_{CAL,tot, #gamma} [GeV]"});

}
}  // namespace study1