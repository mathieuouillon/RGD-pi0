#pragma once

// C++ headers
#include <set>
#include <string>
#include <deque>

// <fmt> headers
#include <fmt/core.h>

// toml++ headers
#include <toml++/toml.h>

// ROOT headers
#include <Math/VectorUtil.h>

// Project headers
#include <hipo4/reader.h>
#include <hipo4/writer.h>
#include <Core/Constantes.hpp>
#include <Core/Helpers.hpp>
#include <Core/Particle.hpp>
#include <Core/ReadBank.hpp>
#include <study1/Histograms.hpp>
#include <study1/models/RGA_inbending_pass1.hpp>
#include <study1/models/RGA_inbending_pass2.hpp>
#include <study1/models/RGA_outbending_pass1.hpp>
#include <study1/models/RGA_outbending_pass2.hpp>
#include <study1/models/RGC_Summer2022_pass1.hpp>

namespace study1 {

class Reader {
   private:
    struct Topology {
        std::vector<Core::Particle> electrons;
        std::vector<Core::Particle> photons;
    };

    struct calo_row_data {
        double pcal_x = 0;
        double pcal_y = 0;
        double pcal_z = 0;
        double ecin_x = 0;
        double ecin_y = 0;
        double ecin_z = 0;
        double ecout_x = 0;
        double ecout_y = 0;
        double ecout_z = 0;
        double pcal_e = 0;
        double pcal_m2u = 0;
        double pcal_m2v = 0;
        double ecin_e = 0;
        double ecin_m2u = 0;
        double ecin_m2v = 0;
        double ecout_e = 0;
        double ecout_m2u = 0;
        double ecout_m2v = 0;
    };

    // clang-format off
    // Map for the GBT Models to use depending on pass and run number
    std::map<std::tuple<int, int, int>, std::function<double(std::vector<float> const &)>> const modelMap = {
      {{5032, 5332, 1}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass1(data); }}, // Fall2018 RGA Inbending
      {{5032, 5332, 2}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass2(data); }}, // Fall2018 RGA Inbending
      {{5333, 5666, 1}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_outbending_pass1(data); }}, // Fall2018 RGA Outbending
      {{5333, 5666, 2}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_outbending_pass2(data); }}, // Fall2018 RGA Outbending
      {{6616, 6783, 1}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass1(data); }}, // Spring2019 RGA Inbending
      {{6616, 6783, 2}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass2(data); }}, // Spring2019 RGA Inbending
      {{6156, 6603, 1}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass1(data); }}, // Spring2019 RGB Inbending
      {{6156, 6603, 2}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass2(data); }}, // Spring2019 RGB Inbending
      {{11093, 11283, 1}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_outbending_pass1(data); }}, // Fall2019 RGB Outbending
      {{11093, 11283, 2}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_outbending_pass2(data); }}, // Fall2019 RGB Outbending
      {{11284, 11300, 1}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass1(data); }}, // Fall2019 RGB BAND Inbending
      {{11284, 11300, 2}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass2(data); }}, // Fall2019 RGB BAND Inbending
      {{11323, 11571, 1}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass1(data); }}, // Spring2020 RGB Inbending
      {{11323, 11571, 2}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGA_inbending_pass2(data); }}, // Spring2020 RGB Inbending
      {{16042, 16772, 1}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGC_Summer2022_pass1(data); }}, // Summer2022 RGC Inbending
      {{16042, 16772, 2}, [](std::vector<float> const &data) { return ApplyCatboostModel_RGC_Summer2022_pass1(data); }} // Summer2022 RGC Inbending (no pass2 currently)
    };
    // clang-format on

    
    const double o_threshold = 0.78; // Threshold value for model predictions
    const int o_pass = 1; // Integer for the event reconstruction pass
    const int number_of_events = 50;

    Histograms& m_histograms;
    const toml::parse_result& m_config;

    const double NaN = std::numeric_limits<double>::quiet_NaN();

    // ****** private methods
    auto get_topology(const hipo::bank& REC_Particle) const -> Topology;
    auto select_electron(const Core::Particle& electron, const hipo::bank& REC_Calorimeter, const hipo::bank& REC_Cherenkov) const -> bool;
    auto select_photons(const std::vector<Core::Particle>& photons, const Core::Particle& electron, const hipo::bank& REC_Calorimeter, double start_time) const -> std::vector<Core::Particle>;
    auto select_photons_with_AI(const hipo::bank& RUN_config, const hipo::bank& REC_Particle, const hipo::bank& REC_Calorimeter) const -> std::vector<Core::Particle>;
    std::map<int, calo_row_data> GetCaloMap(hipo::bank const& bank) const;
    bool Filter(hipo::bank const& particleBank, hipo::bank const& caloBank, std::map<int, calo_row_data> calo_map, int const row, int const runnum) const;
    bool PidPurityPhotonFilter(float const E, float const Epcal, float const theta) const;
    bool ForwardDetectorFilter(float const theta) const;
    ROOT::Math::XYZVector GetParticleCaloVector(calo_row_data calo_row) const;
    bool ClassifyPhoton(std::vector<float> const& input_data, int const runnum) const;
    std::function<double(std::vector<float> const&)> getModelFunction(int runnum) const;

   public:
    // ****** constructors and destructor
    explicit Reader(Histograms& histograms, const toml::parse_result& config);
    ~Reader();

    // ****** public methods
    auto operator()(const std::string& file) const -> void;
};

}  // namespace study1