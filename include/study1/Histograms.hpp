#pragma once

#include <TH1.h>
#include <TH2.h>
#include <ROOT/TThreadedObject.hxx>

using TTO_TH1D = ROOT::TThreadedObject<TH1D>;
using TTO_TH2D = ROOT::TThreadedObject<TH2D>;
using up_TTO_TH1D = const std::unique_ptr<TTO_TH1D>;
using up_TTO_TH2D = const std::unique_ptr<TTO_TH2D>;

namespace study1 {

struct Electron {
    up_TTO_TH1D hist1D_p = std::make_unique<TTO_TH1D>("p_e", "", 200, 0, 11);
    up_TTO_TH1D hist1D_phi = std::make_unique<TTO_TH1D>("phi_e", "", 200, -190, 190);
    up_TTO_TH1D hist1D_theta = std::make_unique<TTO_TH1D>("theta_e", "", 200, 0, 35);
    up_TTO_TH1D hist1D_chi2 = std::make_unique<TTO_TH1D>("chi2_e", "", 200, -6, 6);
    up_TTO_TH1D hist1D_vz = std::make_unique<TTO_TH1D>("vz_e", "", 200, -30, 20);

    up_TTO_TH1D hist1D_p_cut = std::make_unique<TTO_TH1D>("p_e_cut", "", 200, 0, 11);
    up_TTO_TH1D hist1D_phi_cut = std::make_unique<TTO_TH1D>("phi_e_cut", "", 200, -190, 190);
    up_TTO_TH1D hist1D_theta_cut = std::make_unique<TTO_TH1D>("theta_e_cut", "", 200, 0, 35);
    up_TTO_TH1D hist1D_chi2_cut = std::make_unique<TTO_TH1D>("chi2_e_cut", "", 200, -6, 6);
    up_TTO_TH1D hist1D_vz_cut = std::make_unique<TTO_TH1D>("vz_e_cut", "", 200, -30, 20);

    up_TTO_TH2D hist2D_p_vs_SF = std::make_unique<TTO_TH2D>("p_vs_SF", "", 200, 0, 11, 200, 0, 0.4);
    up_TTO_TH2D hist2D_p_vs_SF_cut = std::make_unique<TTO_TH2D>("p_vs_SF_cut", "", 200, 0, 11, 200, 0, 0.4);
    up_TTO_TH2D hist2D_EPCal_vs_ECal = std::make_unique<TTO_TH2D>("EPCal_vs_ECal", "", 200, 0, 0.2, 200, 0, 0.35);
    up_TTO_TH2D hist2D_EPCal_vs_ECal_cut = std::make_unique<TTO_TH2D>("EPCal_vs_ECal_cut", "", 200, 0, 0.2, 200, 0, 0.35);
    up_TTO_TH2D hist2D_EInp_vs_EPCalp = std::make_unique<TTO_TH2D>("EInp_vs_EPCalp", "", 200, 0, 0.3, 200, 0, 0.35);
    up_TTO_TH2D hist2D_EInp_vs_EPCalp_cut = std::make_unique<TTO_TH2D>("EInp_vs_EPCalp_cut", "", 200, 0, 0.3, 200, 0, 0.35);

    up_TTO_TH1D hist1D_v_pcal = std::make_unique<TTO_TH1D>("v_pcal", "", 200, 0, 50);
    up_TTO_TH1D hist1D_v_pcal_cut = std::make_unique<TTO_TH1D>("v_pcal_cut", "", 200, 0, 50);
    up_TTO_TH1D hist1D_w_pcal = std::make_unique<TTO_TH1D>("w_pcal", "", 200, 0, 50);
    up_TTO_TH1D hist1D_w_pcal_cut = std::make_unique<TTO_TH1D>("w_pcal_cut", "", 200, 0, 50);

    up_TTO_TH2D hist2D_x_vs_y_region1_dc = std::make_unique<TTO_TH2D>("x_vs_y_region1_dc", "", 200, -150, 150, 200, -150, 150);
    up_TTO_TH2D hist2D_x_vs_y_region1_dc_cut = std::make_unique<TTO_TH2D>("x_vs_y_region1_dc_cut", "", 200, -150, 150, 200, -150, 150);
    up_TTO_TH2D hist2D_x_vs_y_region1_dc_remove = std::make_unique<TTO_TH2D>("x_vs_y_region1_dc_remove", "", 200, -150, 150, 200, -150, 150);

};

struct Histograms {
    Electron electron;
};

}  // namespace study1