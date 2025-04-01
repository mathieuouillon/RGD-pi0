#pragma once

// C++ headers
#include <limits>

// ROOT headers
#include <TF1.h>

// toml headers
#include <toml++/toml.h>

// Project headers
#include <Plotting/Draw.hpp>
#include <study1/Histograms.hpp>

namespace study1 {
class Drawing {

   private:
    // ****** private members
    Histograms& m_histograms;
    const toml::parse_result& m_config;
    const std::string m_path_electron = "../plots/electron/";
    const std::string m_path_photon = "../plots/photon/";
    const std::string m_path_event = "../plots/event/";
    const std::string m_path_photon_info = "../plots/photon_info/";
    const std::string m_path_photon_info_off = "../plots/photon_info_off/";
    

    const double NaN = std::numeric_limits<double>::quiet_NaN();

    // ****** private methods
    auto fit_time_photon(const std::shared_ptr<TH1D>& hist1D_time_cut_fit) const -> void;
    auto fit_opening_angle(const std::shared_ptr<TH1D>& hist1D_opening_angle_cut_fit) const -> void;
    auto draw_photon_info(const PhotonInfo& info, const std::string& path) -> void;

   public:
    // ****** constructors and destructor
    explicit Drawing(Histograms& histograms, const toml::parse_result& config);

    // ****** public methods
    auto draw_electron() -> void;
    auto draw_photon() -> void;
    auto draw_event() -> void;
};

}  // namespace study1