#include <study1/Reader.hpp>
#include <vector>

namespace study1 {

Reader::Reader(Histograms& histograms, const toml::parse_result& config, const std::vector<int>& pids, const toml::parse_result& sf_cfg)
    : m_histograms(histograms), m_config(config) {
    for (const int pid : pids) {
        particle_collections[pid] = {};
    }

    
        auto sec = sf_cfg["all"];
        auto& mu = m_parameters.mu;
        auto& sig = m_parameters.sigma;

        mu.a = sec["mu"]["a"].value_or(NaN);
        mu.b = sec["mu"]["b"].value_or(NaN);
        mu.c = sec["mu"]["c"].value_or(NaN);

        sig.a = sec["sigma"]["a"].value_or(NaN);
        sig.b = sec["sigma"]["b"].value_or(NaN);
        sig.c = sec["sigma"]["c"].value_or(NaN);
    
}

Reader::~Reader() = default;

auto Reader::operator()(const std::string& file) -> void {

    hipo::hipoeventfile events(file);

    int count = 0;
    for (auto event : events) {
        count++;
        if (count > 10000) break;
        
        // Clear the particle collections for each event
        for (auto& [key, vec] : particle_collections) {
            vec.clear();
        }

        hipo::bank& REC_Particle = event.get_bank("REC::Particle");
        hipo::bank& REC_Calorimeter = event.get_bank("REC::Calorimeter");
        hipo::bank& REC_Cherenkov = event.get_bank("REC::Cherenkov");
        hipo::bank& REC_Event = event.get_bank("REC::Event");
        hipo::bank& REC_Traj = event.get_bank("REC::Traj");
        hipo::bank& REC_Track = event.get_bank("REC::Track");

        if (REC_Particle.getRows() == 0) continue;

        get_topology(REC_Particle);
        const std::vector<Core::Particle>& electrons = particle_collections[11];
        const std::vector<Core::Particle>& photons = particle_collections[22];

        auto possible_electron = Core::find_trigger_electron(electrons);
        if (!possible_electron.has_value()) continue;
        Core::Particle electron = possible_electron.value();

        bool pass_electron_cuts = select_electron(electron, REC_Calorimeter, REC_Traj, REC_Track);
        if (!pass_electron_cuts) continue;

        
    }
}

auto Reader::get_topology(const hipo::bank& REC_Particle) -> void {
    const int rows = REC_Particle.getRows();

    // Pre-allocate for each particle types
    for (auto& [key, vec] : particle_collections) {
        vec.reserve(rows);
    }

    for (int i = 0; i < REC_Particle.getRows(); i++) {
        const int pid = REC_Particle.get<int>("pid", i);

        if (!particle_collections.contains(pid)) continue;

        // Get particle properties only after we know we want this particle
        const int status = REC_Particle.get<int>("status", i);
        const int charge = REC_Particle.get<int>("charge", i);
        const double px = REC_Particle.get<double>("px", i);
        const double py = REC_Particle.get<double>("py", i);
        const double pz = REC_Particle.get<double>("pz", i);
        const double vx = REC_Particle.get<double>("vx", i);
        const double vy = REC_Particle.get<double>("vy", i);
        const double vz = REC_Particle.get<double>("vz", i);
        const double vt = REC_Particle.get<double>("vt", i);
        const double beta = REC_Particle.get<double>("beta", i);
        const double chi2pid = REC_Particle.get<double>("chi2pid", i);

        const double mass = Core::get_mass(pid);
        const double E = Core::compute_energy(px, py, pz, pid);

        // Add particle to the appropriate collection
        particle_collections[pid].emplace_back(pid, status, i, charge, mass, px, py, pz, E, vx, vy, vz, vt, beta, chi2pid);
    }
}

auto Reader::select_electron(const Core::Particle& electron, const hipo::bank& REC_Calorimeter, const hipo::bank& REC_Traj, const hipo::bank& REC_Track) const -> bool {
    bool pass = false;

    const Core::CalorimeterBank calorimeter = Core::read_Calorimeter_bank(REC_Calorimeter, electron.index());
    const Core::DCTrajBank traj = Core::read_Traj_bank(REC_Traj, electron.index());
    const Core::DCTrackBank track_region1 = Core::read_Track_bank(REC_Track, traj.region1.index);
    const Core::DCTrackBank track_region2 = Core::read_Track_bank(REC_Track, traj.region2.index);
    const Core::DCTrackBank track_region3 = Core::read_Track_bank(REC_Track, traj.region3.index);

    double ETot = calorimeter.pcal.energy + calorimeter.inner.energy + calorimeter.outer.energy;
    double SF = ETot / electron.p();

    // Cuts ---------------------------------------------------------------------------------------
    const double vz_min_e = m_config["electron"]["vz_min"].value_or(NaN);
    const double vz_max_e = m_config["electron"]["vz_max"].value_or(NaN);
    //---------------------------------------------------------------------------------------------

    // Get variables for electrons kinematics to apply the cuts -----------------------------------
    // const int sector = calorimeter.pcal.sector;
    const double p_e = electron.p();
    const double chi2_e = electron.chi2pid();
    const double vz_e = electron.vz();
    const double phi_e = electron.phi();
    const double theta_e = electron.theta();
    //---------------------------------------------------------------------------------------------

    // Get the parameters for the SF fit ----------------------------------------------------------
    const double mu = m_parameters.mu.a + m_parameters.mu.b / p_e + m_parameters.mu.c / (p_e * p_e);
    const double sigma = m_parameters.sigma.a + m_parameters.sigma.b / p_e + m_parameters.sigma.c / (p_e * p_e);

    // DC Fiducial volume cuts --------------------------------------------------------------------
    double x = traj.region1.x;
    double y = traj.region1.y;
    double z = traj.region1.z;
    double edge = traj.region1.edge;
    int sec = track_region1.sector;

    // Rotate Z by PI/3*(sector-1)
    double angle = (sec - 1) * M_PI / 3;
    double s = std::sin(angle);
    double c = std::cos(angle);
    double x0 = x;
    double y0 = y;
    x = c * x0 - s * y0;
    y = s * x0 + c * y0;

    double theta_DC = std::atan2(std::sqrt(x * x + y * y), z) * 180.0 / M_PI;
    // --------------------------------------------------------------------------------------------
    
    // Cuts ---------------------------------------------------------------------------------------
    const bool cut_vz = vz_min_e <= vz_e && vz_e <= vz_max_e;
    const bool cut_SF = mu - 3 * sigma < SF && SF < mu + 3 * sigma;
    const bool cut_EInp_vs_EPcalp = (p_e > 4.5 && (calorimeter.inner.energy / p_e + calorimeter.pcal.energy / p_e > 0.2)) || (p_e <= 4.5);
    const bool cut_EPcal = calorimeter.pcal.energy > 0.07;
    const bool cut_v = calorimeter.pcal.lv > 9;
    const bool cut_w = calorimeter.pcal.lw > 9;
    const bool cut_DC_region1 = edge > (2.23035 + 0.04925 * theta_DC + 0.00231 * theta_DC * theta_DC);
    const bool cuts = cut_SF && cut_EInp_vs_EPcalp && cut_EPcal && cut_vz && cut_v && cut_w && cut_DC_region1;
    // --------------------------------------------------------------------------------------------

    // Fill the histograms before any cuts are applied --------------------------------------------
    m_histograms.electron.hist1D_p->Get()->Fill(p_e);
    m_histograms.electron.hist1D_phi->Get()->Fill(phi_e);
    m_histograms.electron.hist1D_theta->Get()->Fill(theta_e);
    m_histograms.electron.hist1D_chi2->Get()->Fill(chi2_e);
    m_histograms.electron.hist1D_vz->Get()->Fill(vz_e);
    m_histograms.electron.hist2D_p_vs_SF->Get()->Fill(p_e, SF);
    if (p_e > 4.5) m_histograms.electron.hist2D_EInp_vs_EPCalp->Get()->Fill(calorimeter.inner.energy / p_e, calorimeter.pcal.energy / p_e);
    m_histograms.electron.hist2D_EPCal_vs_ECal->Get()->Fill(calorimeter.pcal.energy, calorimeter.inner.energy + calorimeter.outer.energy);
    m_histograms.electron.hist1D_v_pcal->Get()->Fill(calorimeter.pcal.lv);
    m_histograms.electron.hist1D_w_pcal->Get()->Fill(calorimeter.pcal.lw);
    m_histograms.electron.hist2D_x_vs_y_region1_dc->Get()->Fill(traj.region1.x, traj.region1.y);
    //---------------------------------------------------------------------------------------------

    // Fill histograms after cuts: ----------------------------------------------------------------
    if (cut_SF && cut_EInp_vs_EPcalp && cut_EPcal) m_histograms.electron.hist1D_vz_cut->Get()->Fill(vz_e);

    if (cut_SF && cut_EInp_vs_EPcalp && cut_EPcal && cut_vz && cut_v && cut_w && !cut_DC_region1)
        m_histograms.electron.hist2D_x_vs_y_region1_dc_remove->Get()->Fill(traj.region1.x, traj.region1.y);

    if (cuts) {
        m_histograms.electron.hist1D_vz_cut->Get()->Fill(vz_e);
        m_histograms.electron.hist1D_p_cut->Get()->Fill(p_e);
        m_histograms.electron.hist1D_phi_cut->Get()->Fill(phi_e);
        m_histograms.electron.hist1D_theta_cut->Get()->Fill(theta_e);
        m_histograms.electron.hist1D_chi2_cut->Get()->Fill(chi2_e);

        m_histograms.electron.hist2D_p_vs_SF_cut->Get()->Fill(p_e, SF);
        m_histograms.electron.hist2D_EInp_vs_EPCalp_cut->Get()->Fill(calorimeter.inner.energy / p_e, calorimeter.pcal.energy / p_e);
        m_histograms.electron.hist2D_EPCal_vs_ECal_cut->Get()->Fill(calorimeter.pcal.energy, calorimeter.inner.energy + calorimeter.outer.energy);
        m_histograms.electron.hist1D_v_pcal_cut->Get()->Fill(calorimeter.pcal.lv);
        m_histograms.electron.hist1D_w_pcal_cut->Get()->Fill(calorimeter.pcal.lw);

        m_histograms.electron.hist2D_x_vs_y_region1_dc_cut->Get()->Fill(traj.region1.x, traj.region1.y);

        pass = true;
    }
    //---------------------------------------------------------------------------------------------

    return pass;
}

}  // namespace study1