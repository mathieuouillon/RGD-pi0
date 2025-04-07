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
};

struct Photon {

    up_TTO_TH1D hist1D_E = std::make_unique<TTO_TH1D>("E_g", "", 200, 0.0, 5.0);
    up_TTO_TH1D hist1D_beta = std::make_unique<TTO_TH1D>("beta_g", "", 200, 0.8, 1.3);
    up_TTO_TH1D hist1D_time = std::make_unique<TTO_TH1D>("time", "", 200, -5, 5);
    up_TTO_TH1D hist1D_vz = std::make_unique<TTO_TH1D>("vz", "", 200, -30, 20);
    up_TTO_TH1D hist1D_opening_angle = std::make_unique<TTO_TH1D>("opening_angle", "", 200, 0, 55);
    up_TTO_TH1D hist1D_delta_vz = std::make_unique<TTO_TH1D>("delta_vz", "", 200, -5.0, 5.0);

    up_TTO_TH2D hist2D_phi_vs_theta = std::make_unique<TTO_TH2D>("phi_vs_theta", "", 200, -190, 190, 200, 0, 35);
    up_TTO_TH2D hist2D_E_vs_vz = std::make_unique<TTO_TH2D>("E_vs_vz", "", 200, 0, 11, 200, -30, 20);

    up_TTO_TH1D hist1D_E_cut = std::make_unique<TTO_TH1D>("E_g_cut", "", 200, 0.0, 5.0);
    up_TTO_TH1D hist1D_beta_cut = std::make_unique<TTO_TH1D>("beta_g_cut", "", 200, 0.8, 1.3);
    up_TTO_TH1D hist1D_delta_sector_cut = std::make_unique<TTO_TH1D>("delta_sector_cut", "", 50, -10, 10);
    up_TTO_TH1D hist1D_time_cut = std::make_unique<TTO_TH1D>("time_cut", "", 200, -5, 5);
    up_TTO_TH1D hist1D_vz_cut = std::make_unique<TTO_TH1D>("vz_cut", "", 200, -30, 20);
    up_TTO_TH1D hist1D_opening_angle_cut = std::make_unique<TTO_TH1D>("opening_angle_cut", "", 200, 0, 55);
    up_TTO_TH1D hist1D_delta_vz_cut = std::make_unique<TTO_TH1D>("delta_vz_cut", "", 200, -5.0, 5.0);

    up_TTO_TH2D hist2D_phi_vs_theta_cut = std::make_unique<TTO_TH2D>("phi_vs_theta_cut", "", 200, -190, 190, 200, 0, 35);
    up_TTO_TH2D hist2D_E_vs_vz_cut = std::make_unique<TTO_TH2D>("E_vs_vz_cut", "", 200, 0, 11, 200, -30, 20);

};


struct Event {

    up_TTO_TH1D hist1D_invariant_mass = std::make_unique<TTO_TH1D>("invariant_mass", "", 200, 0.01, 0.5);
    up_TTO_TH1D hist1D_invariant_mass_cut = std::make_unique<TTO_TH1D>("invariant_mass_cut", "", 200, 0.01, 0.5);

    up_TTO_TH1D hist1D_W = std::make_unique<TTO_TH1D>("W", "", 200, 0, 5);
    up_TTO_TH1D hist1D_Q2 = std::make_unique<TTO_TH1D>("Q2", "", 200, 0, 10);
    up_TTO_TH1D hist1D_nu = std::make_unique<TTO_TH1D>("nu", "", 200, 0, 10);
    up_TTO_TH1D hist1D_zh = std::make_unique<TTO_TH1D>("zh", "", 200, 0, 2);
    up_TTO_TH1D hist1D_y = std::make_unique<TTO_TH1D>("y", "", 200, 0, 2);

    up_TTO_TH1D hist1D_W_cut = std::make_unique<TTO_TH1D>("W_cut", "", 200, 0, 5);
    up_TTO_TH1D hist1D_Q2_cut = std::make_unique<TTO_TH1D>("Q2_cut", "", 200, 0, 10);
    up_TTO_TH1D hist1D_nu_cut = std::make_unique<TTO_TH1D>("nu_cut", "", 200, 0, 10);
    up_TTO_TH1D hist1D_y_cut = std::make_unique<TTO_TH1D>("y_cut", "", 200, 0, 2);
    up_TTO_TH1D hist1D_zh_cut = std::make_unique<TTO_TH1D>("zh_cut", "", 200, 0, 2);

    up_TTO_TH2D hist2D_angle_vs_invariant_mass = std::make_unique<TTO_TH2D>("angle_vs_invariant_mass", "", 200, 0, 70, 200, 0, 0.5);


};

struct Histograms {
    Electron electron;
    Photon photon;
    Event event;
    Event event_mix;
};

}  // namespace study1