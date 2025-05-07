#include <study1/Drawing.hpp>

namespace study1 {

Drawing::Drawing(Histograms& histograms, const toml::parse_result& config)
    : m_histograms(histograms), m_config(config) {
}

auto Drawing::draw_electron_kinematics() -> void {

    // Cuts ---------------------------------------------------------------------------------------
    const double vz_min_e = m_config["electron"]["vz_min"].value_or(NaN);
    const double vz_max_e = m_config["electron"]["vz_max"].value_or(NaN);
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

    Plotting::draw_hist1D(hist1D_p, hist1D_p_cut, m_path_electron, {.label = "p_{e} [GeV/c]"});
    Plotting::draw_hist1D(hist1D_chi2, hist1D_chi2_cut, m_path_electron, {.label = "chi2pid_{e}"});
    Plotting::draw_hist1D(hist1D_phi, hist1D_phi_cut, m_path_electron, {.label = "#phi_{e} [deg.]"});
    Plotting::draw_hist1D(hist1D_theta, hist1D_theta_cut, m_path_electron, {.label = "#theta_{e} [deg.]"});
    Plotting::draw_hist1D(hist1D_vz, hist1D_vz_cut, m_path_electron, {.cuts = {vz_min_e, vz_max_e}, .label = "Vz_{e} [cm]"});

    // Merge 2D histograms ------------------------------------------------------------------------
    const std::shared_ptr<TH2D> hist2D_p_vs_SF = m_histograms.electron.hist2D_p_vs_SF->Merge();
    const std::shared_ptr<TH2D> hist2D_EPCal_vs_ECal = m_histograms.electron.hist2D_EPCal_vs_ECal->Merge();
    const std::shared_ptr<TH2D> hist2D_EInp_vs_EPCalp = m_histograms.electron.hist2D_EInp_vs_EPCalp->Merge();
    const std::shared_ptr<TH1D> hist1D_v_pcal = m_histograms.electron.hist1D_v_pcal->Merge();
    const std::shared_ptr<TH1D> hist1D_w_pcal = m_histograms.electron.hist1D_w_pcal->Merge();

    const std::shared_ptr<TH2D> hist2D_p_vs_SF_cut = m_histograms.electron.hist2D_p_vs_SF_cut->Merge();
    const std::shared_ptr<TH2D> hist2D_EPCal_vs_ECal_cut = m_histograms.electron.hist2D_EPCal_vs_ECal_cut->Merge();
    const std::shared_ptr<TH2D> hist2D_EInp_vs_EPCalp_cut = m_histograms.electron.hist2D_EInp_vs_EPCalp_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_v_pcal_cut = m_histograms.electron.hist1D_v_pcal_cut->Merge();
    const std::shared_ptr<TH1D> hist1D_w_pcal_cut = m_histograms.electron.hist1D_w_pcal_cut->Merge();
    // --------------------------------------------------------------------------------------------

    Plotting::draw_hist2D(hist2D_p_vs_SF, m_path_electron, {.label_x = "p_{e} [GeV/c]", .label_y = "SF"});
    Plotting::draw_hist2D(hist2D_EPCal_vs_ECal, m_path_electron, {.label_x = "E_{PCal} [GeV]", .label_y = "E_{inner} + E_{outer} [GeV]"});
    Plotting::draw_hist2D(hist2D_EInp_vs_EPCalp, m_path_electron, {.label_x = "E_{inner} / p_{e}", .label_y = "E_{PCal} / p_{e}"});
    Plotting::draw_hist2D(hist2D_p_vs_SF_cut, m_path_electron, {.label_x = "p_{e} [GeV/c]", .label_y = "SF"});
    Plotting::draw_hist2D(hist2D_EPCal_vs_ECal_cut, m_path_electron, {.label_x = "E_{PCal} [GeV]", .label_y = "E_{inner} + E_{outer} [GeV]"});
    Plotting::draw_hist2D(hist2D_EInp_vs_EPCalp_cut, m_path_electron, {.label_x = "E_{inner} / p_{e}", .label_y = "E_{PCal} / p_{e}"});
    Plotting::draw_hist1D(hist1D_v_pcal, hist1D_v_pcal_cut, m_path_electron, {.label = "v_{PCal} [cm]"});
    Plotting::draw_hist1D(hist1D_w_pcal, hist1D_w_pcal_cut, m_path_electron, {.label = "w_{PCal} [cm]"});

    // DC Fiducial cuts ---------------------------------------------------------------------------
    const std::shared_ptr<TH2D> hist2D_x_vs_y_region1_dc = m_histograms.electron.hist2D_x_vs_y_region1_dc->Merge();
    const std::shared_ptr<TH2D> hist2D_x_vs_y_region1_dc_cut = m_histograms.electron.hist2D_x_vs_y_region1_dc_cut->Merge();
    const std::shared_ptr<TH2D> hist2D_x_vs_y_region1_dc_remove = m_histograms.electron.hist2D_x_vs_y_region1_dc_remove->Merge();
    // --------------------------------------------------------------------------------------------

    Plotting::draw_hist2D(hist2D_x_vs_y_region1_dc, m_path_electron, {.label_x = "x [cm]", .label_y = "y [cm]"});
    Plotting::draw_hist2D(hist2D_x_vs_y_region1_dc_cut, m_path_electron, {.label_x = "x [cm]", .label_y = "y [cm]"});
    Plotting::draw_hist2D(hist2D_x_vs_y_region1_dc_remove, m_path_electron, {.label_x = "x [cm]", .label_y = "y [cm]"});
}


}  // namespace study1