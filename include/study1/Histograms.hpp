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

struct PhotonInfo {
    up_TTO_TH1D hist1D_E_least = std::make_unique<TTO_TH1D>("E_g_least", "", 200, 0.0, 5.0);
    up_TTO_TH1D hist1D_E_high = std::make_unique<TTO_TH1D>("E_g_high", "", 200, 0.0, 5.0);
    up_TTO_TH1D hist1D_beta_least = std::make_unique<TTO_TH1D>("beta_g_least", "", 200, 0.8, 1.3);
    up_TTO_TH1D hist1D_beta_high = std::make_unique<TTO_TH1D>("beta_g_high", "", 200, 0.8, 1.3);
    up_TTO_TH1D hist1D_time_least = std::make_unique<TTO_TH1D>("time_least", "", 200, -5, 5);
    up_TTO_TH1D hist1D_time_high = std::make_unique<TTO_TH1D>("time_high", "", 200, -5, 5);

    up_TTO_TH1D hist1D_Epcal_high = std::make_unique<TTO_TH1D>("Epcal_high", "", 200, 0, 1);
    up_TTO_TH1D hist1D_Epcal_least = std::make_unique<TTO_TH1D>("Epcal_least", "", 200, 0, 1);
    up_TTO_TH1D hist1D_Ecalo_in_high = std::make_unique<TTO_TH1D>("Ecalo_in_high", "", 200, 0, 1);
    up_TTO_TH1D hist1D_Ecalo_in_least = std::make_unique<TTO_TH1D>("Ecalo_in_least", "", 200, 0, 1);
    up_TTO_TH1D hist1D_Ecalo_out_high = std::make_unique<TTO_TH1D>("Ecalo_out_high", "", 200, 0, 1);
    up_TTO_TH1D hist1D_Ecalo_out_least = std::make_unique<TTO_TH1D>("Ecalo_out_least", "", 200, 0, 1);
    up_TTO_TH1D hist1D_Ecalo_m2u_high = std::make_unique<TTO_TH1D>("Ecalo_m2u_high", "", 200, 0, 50);
    up_TTO_TH1D hist1D_Ecalo_m2u_least = std::make_unique<TTO_TH1D>("Ecalo_m2u_least", "", 200, 0, 50);
    up_TTO_TH1D hist1D_Ecalo_m2v_high = std::make_unique<TTO_TH1D>("Ecalo_m2v_high", "", 200, 0, 50);
    up_TTO_TH1D hist1D_Ecalo_m2v_least = std::make_unique<TTO_TH1D>("Ecalo_m2v_least", "", 200, 0, 50);
    up_TTO_TH1D hist1D_Ecalo_m2w_high = std::make_unique<TTO_TH1D>("Ecalo_m2w_high", "", 200, 0, 50);
    up_TTO_TH1D hist1D_Ecalo_m2w_least = std::make_unique<TTO_TH1D>("Ecalo_m2w_least", "", 200, 0, 50);
    up_TTO_TH1D hist1D_Ecalo_m3u_high = std::make_unique<TTO_TH1D>("Ecalo_m3u_high", "", 200, -10, 10);
    up_TTO_TH1D hist1D_Ecalo_m3u_least = std::make_unique<TTO_TH1D>("Ecalo_m3u_least", "", 200, -10, 10);
    up_TTO_TH1D hist1D_Ecalo_m3v_high = std::make_unique<TTO_TH1D>("Ecalo_m3v_high", "", 200, -10, 10);
    up_TTO_TH1D hist1D_Ecalo_m3v_least = std::make_unique<TTO_TH1D>("Ecalo_m3v_least", "", 200, -10, 10);
    up_TTO_TH1D hist1D_Ecalo_m3w_high = std::make_unique<TTO_TH1D>("Ecalo_m3w_high", "", 200, -10, 10);
    up_TTO_TH1D hist1D_Ecalo_m3w_least = std::make_unique<TTO_TH1D>("Ecalo_m3w_least", "", 200, -10, 10);

    up_TTO_TH1D hist1D_Ecalo_tot_high = std::make_unique<TTO_TH1D>("Ecalo_tot_high", "", 200, 0, 4);
    up_TTO_TH1D hist1D_Ecalo_tot_least = std::make_unique<TTO_TH1D>("Ecalo_tot_least", "", 200, 0, 4);
};


struct Event {

    up_TTO_TH1D hist1D_invariant_mass = std::make_unique<TTO_TH1D>("invariant_mass", "", 200, 0.01, 0.5);
    up_TTO_TH1D hist1D_invariant_mass_cut = std::make_unique<TTO_TH1D>("invariant_mass_cut", "", 200, 0.01, 0.5);

    up_TTO_TH1D hist1D_W = std::make_unique<TTO_TH1D>("W", "", 200, 0, 5);
    up_TTO_TH1D hist1D_Q2 = std::make_unique<TTO_TH1D>("Q2", "", 200, 0, 10);
    up_TTO_TH1D hist1D_nu = std::make_unique<TTO_TH1D>("nu", "", 200, 0, 10);
    up_TTO_TH1D hist1D_zh = std::make_unique<TTO_TH1D>("zh", "", 200, 0, 2);

    up_TTO_TH1D hist1D_W_cut = std::make_unique<TTO_TH1D>("W_cut", "", 200, 0, 5);
    up_TTO_TH1D hist1D_Q2_cut = std::make_unique<TTO_TH1D>("Q2_cut", "", 200, 0, 10);
    up_TTO_TH1D hist1D_nu_cut = std::make_unique<TTO_TH1D>("nu_cut", "", 200, 0, 10);
    up_TTO_TH1D hist1D_zh_cut = std::make_unique<TTO_TH1D>("zh_cut", "", 200, 0, 2);

    // 2D histograms
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_Q2_cut = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_Q2_cut", "", 200, 0, 0.5, 200, 0, 10);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_zh_cut = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_zh_cut", "", 200, 0, 0.5, 200, 0, 2);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_nu_cut = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_nu_cut", "", 200, 0, 0.5, 200, 0, 10);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_W_cut = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_W", "", 200, 0, 0.5, 200, 0.0, 5);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_y = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_y", "", 200, 0, 0.5, 200, 0.0, 1);
    
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_opening_angle_least = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_opening_angle_least", "", 200, 0, 0.5, 200, 0.0, 55);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_opening_angle_high = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_opening_angle_high", "", 200, 0, 0.5, 200, 0.0, 55);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_E_least = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_E_least", "", 200, 0, 0.5, 200, 0.0, 10);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_E_high = std::make_unique<TTO_TH2D>("invariant_mass_cut_vsEe_high", "", 200, 0, 0.5, 200, 0.0, 10);

    up_TTO_TH2D hist2D_invariant_mass_cut_vs_p_e = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_p_e", "", 200, 0, 0.5, 200, 0.0, 10);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_theta_e = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_theta_e", "", 200, 0, 0.5, 200, 0.0, 45);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_phi_e = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_phi_e", "", 200, 0, 0.5, 200, -180, 180);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_chi2_e = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_chi2_e", "", 200, 0, 0.5, 200, -6, 6);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_vz_e = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_vz_e", "", 200, 0, 0.5, 200, -30, 20);

    up_TTO_TH2D hist2D_invariant_mass_cut_vs_beta_g_least = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_beta_g_least", "", 200, 0, 0.5, 200, 0.8, 1.3);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_beta_g_high = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_beta_g_high", "", 200, 0, 0.5, 200, 0.8, 1.3);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_time_g_least = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_time_g_least", "", 200, 0, 0.5, 200, -5, 5);
    up_TTO_TH2D hist2D_invariant_mass_cut_vs_time_g_high = std::make_unique<TTO_TH2D>("invariant_mass_cut_vs_time_g_high", "", 200, 0, 0.5, 200, -5, 5);


    // Photon information
    PhotonInfo photon_info;
    PhotonInfo photon_info_off;




};

struct Histograms {
    Electron electron;
    Photon photon;
    Event event;
};

}  // namespace study1