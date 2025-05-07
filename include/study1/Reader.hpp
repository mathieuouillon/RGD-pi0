#pragma once

// C++ headers
#include <string>

// <fmt> headers
#include <fmt/core.h>

// toml++ headers
#include <toml++/toml.h>

// Project headers
#include <hipo4/hipoeventiterator.h>
#include <hipo4/reader.h>
#include <Core/Constantes.hpp>
#include <Core/Helpers.hpp>
#include <Core/Particle.hpp>
#include <Core/ReadBank.hpp>
#include <study1/Histograms.hpp>
#include <vector>

namespace study1 {

class Reader {
   private:
    struct Parameter {
        double a;
        double b;
        double c;
    };

    struct Parameters {
        Parameter mu;
        Parameter sigma;
    };

    // ****** private variables
    Histograms& m_histograms;
    const toml::parse_result& m_config;
    std::unordered_map<int, std::vector<Core::Particle>> particle_collections;
    const int m_num_sectors = 6;
    Parameters m_parameters;

    // ****** private constants
    static constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

    // ****** private methods
    auto select_electron(const Core::Particle& electron, const hipo::bank& REC_Calorimeter, const hipo::bank& REC_Traj, const hipo::bank& REC_Track) const -> bool;
    auto get_topology(const hipo::bank& REC_Particle) -> void;

   public:
    // ****** constructors and destructor
    explicit Reader(Histograms& histograms, const toml::parse_result& config, const std::vector<int>& pids, const toml::parse_result& sf_cfg);
    ~Reader();

    // ****** public methods
    auto operator()(const std::string& file) -> void;
};

}  // namespace study1