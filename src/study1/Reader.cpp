#include <study1/Reader.hpp>

namespace study1 {

Reader::Reader(Histograms& histograms, const toml::parse_result& config)
    : m_histograms(histograms), m_config(config) {
}

Reader::~Reader() = default;

auto Reader::operator()(const std::string& file) const -> void {

    bool mode_write = false;

    auto dict = hipo::dictionary();
    auto reader = hipo::reader(file, dict);
    auto hipo_event = hipo::event();

    auto writer = hipo::writer();
    auto dict_writer = hipo::dictionary();
    reader.readDictionary(dict_writer);

    for (const std::string& s : dict_writer.getSchemaList())
        writer.getDictionary().addSchema(dict_writer.getSchema(s.c_str()));

    hipo::schema schemaPart("ANAL::Particle", 100, 1);
    schemaPart.parse("pid/S,px/F,py/F,pz/F");
    writer.getDictionary().addSchema(schemaPart);

    int pos = file.find("skim_run_");
    std::vector<std::string> dst_schema = {"RUN::config", "REC::Event", "REC::Particle", "REC::Calorimeter",
                                           "REC::ForwardTagger", "REC::Scintillator", "REC::Track", "REC::CovMat",
                                           "REC::Traj", "REC::Cherenkov", "RUN::config", "RUN::rf"};
    if (mode_write) writer.open(std::string("../data/photon_" + file.substr(pos + 5)).c_str());

    auto REC_Event = hipo::bank(dict.getSchema("REC::Event"));
    auto RUN_config = hipo::bank(dict.getSchema("RUN::config"));
    auto REC_Particle = hipo::bank(dict.getSchema("REC::Particle"));
    auto REC_Calorimeter = hipo::bank(dict.getSchema("REC::Calorimeter"));
    auto REC_Cherenkov = hipo::bank(dict.getSchema("REC::Cherenkov"));

    while (reader.next(hipo_event, REC_Event, RUN_config, REC_Particle, REC_Calorimeter, REC_Cherenkov)) {
        if (REC_Particle.getRows() == 0) continue;

        double start_time = REC_Event.getFloat("startTime", 0);

        Topology topology = get_topology(REC_Particle);

        auto possible_electron = Core::find_trigger_electron(topology.electrons);
        if (!possible_electron.has_value()) continue;
        Core::Particle electron = possible_electron.value();

        bool pass_electron_cuts = select_electron(electron, REC_Calorimeter, REC_Cherenkov);
        if (!pass_electron_cuts) continue;

        // Select photons using AI
        std::vector<Core::Particle> good_photons = select_photons_with_AI(RUN_config, REC_Particle, REC_Calorimeter);
        if (good_photons.size() < 2) continue;

        // std::vector<Core::Particle> good_photons = select_photons(topology.photons, electron, REC_Calorimeter, start_time);

        for (Core::Particle& photon : good_photons) {
            Core::CalorimeterBank calo_info = Core::read_Calorimeter_bank(REC_Calorimeter, photon.index());
            double time = calo_info.pcal.time - start_time - calo_info.pcal.path / 30;

            double opening_angle = electron.theta(photon);

            m_histograms.photon.hist1D_E->Get()->Fill(photon.E());
            m_histograms.photon.hist1D_beta->Get()->Fill(photon.beta());
            m_histograms.photon.hist1D_time->Get()->Fill(time);
            m_histograms.photon.hist1D_vz->Get()->Fill(photon.vz());
            m_histograms.photon.hist1D_opening_angle->Get()->Fill(opening_angle);
            m_histograms.photon.hist1D_delta_vz->Get()->Fill(electron.vz() - photon.vz());

            m_histograms.photon.hist2D_phi_vs_theta->Get()->Fill(photon.phi(), photon.theta());
            m_histograms.photon.hist2D_E_vs_vz->Get()->Fill(photon.E(), photon.vz());

            m_histograms.photon.hist1D_time_cut->Get()->Fill(time);
            m_histograms.photon.hist1D_E_cut->Get()->Fill(photon.E());
            m_histograms.photon.hist1D_opening_angle_cut->Get()->Fill(opening_angle);

            m_histograms.photon.hist1D_vz_cut->Get()->Fill(photon.vz());
            m_histograms.photon.hist1D_beta_cut->Get()->Fill(photon.beta());

            m_histograms.photon.hist2D_phi_vs_theta_cut->Get()->Fill(photon.phi(), photon.theta());
            m_histograms.photon.hist2D_E_vs_vz_cut->Get()->Fill(photon.E(), photon.vz());

            m_histograms.photon.hist1D_delta_vz_cut->Get()->Fill(electron.vz() - photon.vz());
        }

        // Get the photon from another event
        std::vector<Core::Particle> photon_mix_queue;
        int count = 0;
        while (photon_mix_queue.size() < number_of_events && reader.eventNumber() < reader.getEntries() - number_of_events) {
            auto REC_Particle_main = hipo::bank(dict.getSchema("REC::Particle"));
            auto REC_Calorimeter_main = hipo::bank(dict.getSchema("REC::Calorimeter"));

            REC_Particle_main.setRows(REC_Particle.getRows());
            REC_Calorimeter_main.setRows(REC_Calorimeter.getRows());

            // Read the next event
            for (int row = 0; row < REC_Particle.getRows(); row++) {
                int pid = REC_Particle.getInt("pid", row);
                short status = REC_Particle.getShort("status", row);
                short charge = REC_Particle.getByte("charge", row);

                float px = REC_Particle.getFloat("px", row);
                float py = REC_Particle.getFloat("py", row);
                float pz = REC_Particle.getFloat("pz", row);
                float vx = REC_Particle.getFloat("vx", row);
                float vy = REC_Particle.getFloat("vy", row);
                float vz = REC_Particle.getFloat("vz", row);
                float vt = REC_Particle.getFloat("vt", row);
                float beta = REC_Particle.getFloat("beta", row);
                float chi2pid = REC_Particle.getFloat("chi2pid", row);

                REC_Particle_main.putInt("pid", row, pid);
                REC_Particle_main.putShort("status", row, status);
                REC_Particle_main.putByte("charge", row, charge);
                REC_Particle_main.putFloat("px", row, px);
                REC_Particle_main.putFloat("py", row, py);
                REC_Particle_main.putFloat("pz", row, pz);
                REC_Particle_main.putFloat("vx", row, vx);
                REC_Particle_main.putFloat("vy", row, vy);
                REC_Particle_main.putFloat("vz", row, vz);
                REC_Particle_main.putFloat("vt", row, vt);
                REC_Particle_main.putFloat("beta", row, beta);
                REC_Particle_main.putFloat("chi2pid", row, chi2pid);
            }

            for (int row = 0; row < REC_Calorimeter.getRows(); row++) {
                short index = REC_Calorimeter.getShort("index", row);
                short pindex = REC_Calorimeter.getShort("pindex", row);
                short detector = REC_Calorimeter.getByte("detector", row);
                short sector = REC_Calorimeter.getByte("sector", row);
                short layer = REC_Calorimeter.getByte("layer", row);
                float energy = REC_Calorimeter.getFloat("energy", row);
                float time = REC_Calorimeter.getFloat("time", row);
                float path = REC_Calorimeter.getFloat("path", row);
                float chi2 = REC_Calorimeter.getFloat("chi2", row);
                float x = REC_Calorimeter.getFloat("x", row);
                float y = REC_Calorimeter.getFloat("y", row);
                float z = REC_Calorimeter.getFloat("z", row);
                float hx = REC_Calorimeter.getFloat("hx", row);
                float hy = REC_Calorimeter.getFloat("hy", row);
                float hz = REC_Calorimeter.getFloat("hz", row);
                float lu = REC_Calorimeter.getFloat("lu", row);
                float lv = REC_Calorimeter.getFloat("lv", row);
                float lw = REC_Calorimeter.getFloat("lw", row);
                float du = REC_Calorimeter.getFloat("du", row);
                float dv = REC_Calorimeter.getFloat("dv", row);
                float dw = REC_Calorimeter.getFloat("dw", row);
                float m2u = REC_Calorimeter.getFloat("m2u", row);
                float m2v = REC_Calorimeter.getFloat("m2v", row);
                float m2w = REC_Calorimeter.getFloat("m2w", row);
                float m3u = REC_Calorimeter.getFloat("m3u", row);
                float m3v = REC_Calorimeter.getFloat("m3v", row);
                float m3w = REC_Calorimeter.getFloat("m3w", row);
                short status = REC_Calorimeter.getShort("status", row);

                REC_Calorimeter_main.putShort("index", row, index);
                REC_Calorimeter_main.putShort("pindex", row, pindex);
                REC_Calorimeter_main.putByte("detector", row, detector);
                REC_Calorimeter_main.putByte("sector", row, sector);
                REC_Calorimeter_main.putByte("layer", row, layer);
                REC_Calorimeter_main.putFloat("energy", row, energy);
                REC_Calorimeter_main.putFloat("time", row, time);
                REC_Calorimeter_main.putFloat("path", row, path);
                REC_Calorimeter_main.putFloat("chi2", row, chi2);
                REC_Calorimeter_main.putFloat("x", row, x);
                REC_Calorimeter_main.putFloat("y", row, y);
                REC_Calorimeter_main.putFloat("z", row, z);
                REC_Calorimeter_main.putFloat("hx", row, hx);
                REC_Calorimeter_main.putFloat("hy", row, hy);
                REC_Calorimeter_main.putFloat("hz", row, hz);
                REC_Calorimeter_main.putFloat("lu", row, lu);
                REC_Calorimeter_main.putFloat("lv", row, lv);
                REC_Calorimeter_main.putFloat("lw", row, lw);
                REC_Calorimeter_main.putFloat("du", row, du);
                REC_Calorimeter_main.putFloat("dv", row, dv);
                REC_Calorimeter_main.putFloat("dw", row, dw);
                REC_Calorimeter_main.putFloat("m2u", row, m2u);
                REC_Calorimeter_main.putFloat("m2v", row, m2v);
                REC_Calorimeter_main.putFloat("m2w", row, m2w);
                REC_Calorimeter_main.putFloat("m3u", row, m3u);
                REC_Calorimeter_main.putFloat("m3v", row, m3v);
                REC_Calorimeter_main.putFloat("m3w", row, m3w);
                REC_Calorimeter_main.putShort("status", row, status);
            }

            // fmt::println("----------------------------------------------------------------------");
            // fmt::println("----------------------------------------------------------------------");
            // REC_Calorimeter.show();
            reader.next(hipo_event, REC_Event, RUN_config, REC_Particle, REC_Calorimeter, REC_Cherenkov);
            Topology topology_mix = get_topology(REC_Particle);

            // Get all photons
            std::vector<int> all_photons_rows;
            for (int i = 0; i < REC_Particle.getRows(); i++) {
                int pid = REC_Particle.getInt("pid", i);
                if (pid == 22) {
                    all_photons_rows.push_back(i);
                }
            }

            int nb_rows_calo = 0;
            for (int i = 0; i < REC_Calorimeter.getRows(); i++) {
                for (const auto row : all_photons_rows) {
                    if (row == REC_Calorimeter.getInt("pindex", i)) {
                        nb_rows_calo++;
                    }
                }
            }

            // fmt::println("----------------------------------------------------------------------");

            auto REC_Particle_main2 = hipo::bank(dict.getSchema("REC::Particle"));
            auto REC_Calorimeter_main2 = hipo::bank(dict.getSchema("REC::Calorimeter"));

            REC_Particle_main2.setRows(REC_Particle_main.getRows() + all_photons_rows.size());
            REC_Calorimeter_main2.setRows(REC_Calorimeter_main.getRows() + nb_rows_calo);

            for (int row = 0; row < REC_Particle_main.getRows(); row++) {
                int pid = REC_Particle_main.getInt("pid", row);
                short status = REC_Particle_main.getShort("status", row);
                short charge = REC_Particle_main.getByte("charge", row);

                float px = REC_Particle_main.getFloat("px", row);
                float py = REC_Particle_main.getFloat("py", row);
                float pz = REC_Particle_main.getFloat("pz", row);
                float vx = REC_Particle_main.getFloat("vx", row);
                float vy = REC_Particle_main.getFloat("vy", row);
                float vz = REC_Particle_main.getFloat("vz", row);
                float vt = REC_Particle_main.getFloat("vt", row);
                float beta = REC_Particle_main.getFloat("beta", row);
                float chi2pid = REC_Particle_main.getFloat("chi2pid", row);

                REC_Particle_main2.putInt("pid", row, pid);
                REC_Particle_main2.putShort("status", row, status);
                REC_Particle_main2.putByte("charge", row, charge);
                REC_Particle_main2.putFloat("px", row, px);
                REC_Particle_main2.putFloat("py", row, py);
                REC_Particle_main2.putFloat("pz", row, pz);
                REC_Particle_main2.putFloat("vx", row, vx);
                REC_Particle_main2.putFloat("vy", row, vy);
                REC_Particle_main2.putFloat("vz", row, vz);
                REC_Particle_main2.putFloat("vt", row, vt);
                REC_Particle_main2.putFloat("beta", row, beta);
                REC_Particle_main2.putFloat("chi2pid", row, chi2pid);
            }

            for (int row = 0; row < REC_Calorimeter_main.getRows(); row++) {
                short index = REC_Calorimeter_main.getShort("index", row);
                short pindex = REC_Calorimeter_main.getShort("pindex", row);
                short detector = REC_Calorimeter_main.getByte("detector", row);
                short sector = REC_Calorimeter_main.getByte("sector", row);
                short layer = REC_Calorimeter_main.getByte("layer", row);
                float energy = REC_Calorimeter_main.getFloat("energy", row);
                float time = REC_Calorimeter_main.getFloat("time", row);
                float path = REC_Calorimeter_main.getFloat("path", row);
                float chi2 = REC_Calorimeter_main.getFloat("chi2", row);
                float x = REC_Calorimeter_main.getFloat("x", row);
                float y = REC_Calorimeter_main.getFloat("y", row);
                float z = REC_Calorimeter_main.getFloat("z", row);
                float hx = REC_Calorimeter_main.getFloat("hx", row);
                float hy = REC_Calorimeter_main.getFloat("hy", row);
                float hz = REC_Calorimeter_main.getFloat("hz", row);
                float lu = REC_Calorimeter_main.getFloat("lu", row);
                float lv = REC_Calorimeter_main.getFloat("lv", row);
                float lw = REC_Calorimeter_main.getFloat("lw", row);
                float du = REC_Calorimeter_main.getFloat("du", row);
                float dv = REC_Calorimeter_main.getFloat("dv", row);
                float dw = REC_Calorimeter_main.getFloat("dw", row);
                float m2u = REC_Calorimeter_main.getFloat("m2u", row);
                float m2v = REC_Calorimeter_main.getFloat("m2v", row);
                float m2w = REC_Calorimeter_main.getFloat("m2w", row);
                float m3u = REC_Calorimeter_main.getFloat("m3u", row);
                float m3v = REC_Calorimeter_main.getFloat("m3v", row);
                float m3w = REC_Calorimeter_main.getFloat("m3w", row);
                short status = REC_Calorimeter_main.getShort("status", row);

                REC_Calorimeter_main2.putShort("index", row, index);
                REC_Calorimeter_main2.putShort("pindex", row, pindex);
                REC_Calorimeter_main2.putByte("detector", row, detector);
                REC_Calorimeter_main2.putByte("sector", row, sector);
                REC_Calorimeter_main2.putByte("layer", row, layer);
                REC_Calorimeter_main2.putFloat("energy", row, energy);
                REC_Calorimeter_main2.putFloat("time", row, time);
                REC_Calorimeter_main2.putFloat("path", row, path);
                REC_Calorimeter_main2.putFloat("chi2", row, chi2);
                REC_Calorimeter_main2.putFloat("x", row, x);
                REC_Calorimeter_main2.putFloat("y", row, y);
                REC_Calorimeter_main2.putFloat("z", row, z);
                REC_Calorimeter_main2.putFloat("hx", row, hx);
                REC_Calorimeter_main2.putFloat("hy", row, hy);
                REC_Calorimeter_main2.putFloat("hz", row, hz);
                REC_Calorimeter_main2.putFloat("lu", row, lu);
                REC_Calorimeter_main2.putFloat("lv", row, lv);
                REC_Calorimeter_main2.putFloat("lw", row, lw);
                REC_Calorimeter_main2.putFloat("du", row, du);
                REC_Calorimeter_main2.putFloat("dv", row, dv);
                REC_Calorimeter_main2.putFloat("dw", row, dw);
                REC_Calorimeter_main2.putFloat("m2u", row, m2u);
                REC_Calorimeter_main2.putFloat("m2v", row, m2v);
                REC_Calorimeter_main2.putFloat("m2w", row, m2w);
                REC_Calorimeter_main2.putFloat("m3u", row, m3u);
                REC_Calorimeter_main2.putFloat("m3v", row, m3v);
                REC_Calorimeter_main2.putFloat("m3w", row, m3w);
                REC_Calorimeter_main2.putShort("status", row, status);
            }

            // REC_Particle.show();
            // REC_Calorimeter.show();
            // REC_Calorimeter_main2.show();
            int ii = REC_Particle_main.getRows();
            int iii = REC_Calorimeter_main.getRows();
            for (const auto row : all_photons_rows) {
                // fmt::println("ii: {}", ii);
                int pid = REC_Particle.getInt("pid", row);
                short status = REC_Particle.getShort("status", row);
                short charge = REC_Particle.getByte("charge", row);

                float px = REC_Particle.getFloat("px", row);
                float py = REC_Particle.getFloat("py", row);
                float pz = REC_Particle.getFloat("pz", row);
                float vx = REC_Particle.getFloat("vx", row);
                float vy = REC_Particle.getFloat("vy", row);
                float vz = REC_Particle.getFloat("vz", row);
                float vt = REC_Particle.getFloat("vt", row);
                float beta = REC_Particle.getFloat("beta", row);
                float chi2pid = REC_Particle.getFloat("chi2pid", row);

                REC_Particle_main2.putInt("pid", ii, pid);
                REC_Particle_main2.putShort("status", ii, status);
                REC_Particle_main2.putByte("charge", ii, charge);
                REC_Particle_main2.putFloat("px", ii, px);
                REC_Particle_main2.putFloat("py", ii, py);
                REC_Particle_main2.putFloat("pz", ii, pz);
                REC_Particle_main2.putFloat("vx", ii, vx);
                REC_Particle_main2.putFloat("vy", ii, vy);
                REC_Particle_main2.putFloat("vz", ii, vz);
                REC_Particle_main2.putFloat("vt", ii, vt);
                REC_Particle_main2.putFloat("beta", ii, beta);
                REC_Particle_main2.putFloat("chi2pid", ii, chi2pid);
                ii++;
            }
            

            for (int i = 0; i < REC_Calorimeter.getRows(); i++) {
                for (const auto row : all_photons_rows) {
                    if (row == REC_Calorimeter.getInt("pindex", i)) {
                        short index = REC_Calorimeter.getShort("index", i);
                        short pindex = REC_Calorimeter.getShort("pindex", i);
                        short detector = REC_Calorimeter.getByte("detector", i);
                        short sector = REC_Calorimeter.getByte("sector", i);
                        short layer = REC_Calorimeter.getByte("layer", i);
                        float energy = REC_Calorimeter.getFloat("energy", i);
                        float time = REC_Calorimeter.getFloat("time", i);
                        float path = REC_Calorimeter.getFloat("path", i);
                        float chi2 = REC_Calorimeter.getFloat("chi2", i);
                        float x = REC_Calorimeter.getFloat("x", i);
                        float y = REC_Calorimeter.getFloat("y", i);
                        float z = REC_Calorimeter.getFloat("z", i);
                        float hx = REC_Calorimeter.getFloat("hx", i);
                        float hy = REC_Calorimeter.getFloat("hy", i);
                        float hz = REC_Calorimeter.getFloat("hz", i);
                        float lu = REC_Calorimeter.getFloat("lu", i);
                        float lv = REC_Calorimeter.getFloat("lv", i);
                        float lw = REC_Calorimeter.getFloat("lw", i);
                        float du = REC_Calorimeter.getFloat("du", i);
                        float dv = REC_Calorimeter.getFloat("dv", i);
                        float dw = REC_Calorimeter.getFloat("dw", i);
                        float m2u = REC_Calorimeter.getFloat("m2u", i);
                        float m2v = REC_Calorimeter.getFloat("m2v", i);
                        float m2w = REC_Calorimeter.getFloat("m2w", i);
                        float m3u = REC_Calorimeter.getFloat("m3u", i);
                        float m3v = REC_Calorimeter.getFloat("m3v", i);
                        float m3w = REC_Calorimeter.getFloat("m3w", i);
                        short status = REC_Calorimeter.getShort("status", i);
        
                        REC_Calorimeter_main2.putShort("index", iii, index);
                        REC_Calorimeter_main2.putShort("pindex", iii, pindex);
                        REC_Calorimeter_main2.putByte("detector", iii, detector);
                        REC_Calorimeter_main2.putByte("sector", iii, sector);
                        REC_Calorimeter_main2.putByte("layer", iii, layer);
                        REC_Calorimeter_main2.putFloat("energy", iii, energy);
                        REC_Calorimeter_main2.putFloat("time", iii, time);
                        REC_Calorimeter_main2.putFloat("path", iii, path);
                        REC_Calorimeter_main2.putFloat("chi2", iii, chi2);
                        REC_Calorimeter_main2.putFloat("x", iii, x);
                        REC_Calorimeter_main2.putFloat("y", iii, y);
                        REC_Calorimeter_main2.putFloat("z", iii, z);
                        REC_Calorimeter_main2.putFloat("hx", iii, hx);
                        REC_Calorimeter_main2.putFloat("hy", iii, hy);
                        REC_Calorimeter_main2.putFloat("hz", iii, hz);
                        REC_Calorimeter_main2.putFloat("lu", iii, lu);
                        REC_Calorimeter_main2.putFloat("lv", iii, lv);
                        REC_Calorimeter_main2.putFloat("lw", iii, lw);
                        REC_Calorimeter_main2.putFloat("du", iii, du);
                        REC_Calorimeter_main2.putFloat("dv", iii, dv);
                        REC_Calorimeter_main2.putFloat("dw", iii, dw);
                        REC_Calorimeter_main2.putFloat("m2u", iii, m2u);
                        REC_Calorimeter_main2.putFloat("m2v", iii, m2v);
                        REC_Calorimeter_main2.putFloat("m2w", iii, m2w);
                        REC_Calorimeter_main2.putFloat("m3u", iii, m3u);
                        REC_Calorimeter_main2.putFloat("m3v", iii, m3v);
                        REC_Calorimeter_main2.putFloat("m3w", iii, m3w);
                        REC_Calorimeter_main2.putShort("status", iii, status);



                    }
                }
            }
            
            // REC_Calorimeter_main2.show();
            
            
            std::vector<Core::Particle> good_photons_mix = select_photons_with_AI(RUN_config, REC_Particle_main2, REC_Calorimeter_main2);
            
            // Check if the photons are not for the main event
            std::vector<Core::Particle> good_photons_mix_2;
            for (const auto& photon : good_photons_mix) {
                for (const auto& photon2 : good_photons) {
                    if (photon.index() == photon2.index()) {
                        good_photons_mix_2.push_back(photon);
                    }
                }
            }
            
            
            if (good_photons_mix_2.size() < 2) continue;

            for (Core::Particle& photon : good_photons_mix_2) {
                photon_mix_queue.push_back(photon);
            }
            count++;
        }

        for (size_t i = 0; i < count; i++) {
            reader.previous(hipo_event, REC_Event, RUN_config, REC_Particle, REC_Calorimeter, REC_Cherenkov);
        }

        std::vector<std::pair<Core::Particle, Core::Particle>> pair_photon = Core::generate_unique_pairs(good_photons);

        std::set<Core::Particle> set_photons_pass_event_cuts;
        for (const auto& [photon1, photon2] : pair_photon) {
            // Cuts ---------------------------------------------------------------------------------------
            const double Q2_min = m_config["event"]["Q2_min"].value_or(NaN);
            const double W_min = m_config["event"]["W_min"].value_or(NaN);
            const double zh_min = m_config["event"]["zh_min"].value_or(NaN);
            const double zh_max = m_config["event"]["zh_max"].value_or(NaN);
            const double y_max = m_config["event"]["y_max"].value_or(NaN);
            //---------------------------------------------------------------------------------------------

            ROOT::Math::PxPyPzEVector k1 = {0, 0, 10.53, 10.53};
            ROOT::Math::PxPyPzEVector q1 = k1 - electron.PxPyPzEVector();
            ROOT::Math::PxPyPzEVector p1 = {0, 0, 0, Core::Constantes::ProtonMass};
            ROOT::Math::PxPyPzEVector h = photon1.PxPyPzEVector() + photon2.PxPyPzEVector();

            double M = Core::Constantes::ProtonMass;
            double nu = k1.E() - electron.PxPyPzEVector().E();
            double y = nu / k1.E();
            double Q2 = 2. * q1.E() * electron.E() * (1 - std::cos(electron.PxPyPzEVector().Theta()));
            double W = std::sqrt(M * M + 2 * M * nu - Q2);
            double zh = h.E() / nu;
            double xB = Q2 / (2. * M * nu);
            double pT = h.Pt();
            double m = ROOT::Math::VectorUtil::InvariantMass(photon1.PxPyPzEVector(), photon2.PxPyPzEVector());

            Core::CalorimeterBank calo_info_1 = Core::read_Calorimeter_bank(REC_Calorimeter, photon1.index());
            int sector_1 = calo_info_1.pcal.sector;
            Core::CalorimeterBank calo_info_2 = Core::read_Calorimeter_bank(REC_Calorimeter, photon2.index());
            int sector_2 = calo_info_2.pcal.sector;

            m_histograms.event.hist1D_W->Get()->Fill(W);
            m_histograms.event.hist1D_Q2->Get()->Fill(Q2);
            m_histograms.event.hist1D_nu->Get()->Fill(nu);
            m_histograms.event.hist1D_zh->Get()->Fill(zh);
            m_histograms.event.hist1D_y->Get()->Fill(y);
            m_histograms.event.hist1D_invariant_mass->Get()->Fill(m);

            const bool cut_Q2 = Q2 > Q2_min;
            const bool cut_W = W > W_min;
            const bool cut_zh = zh_min < zh && zh < zh_max;
            const bool cut_y = y < y_max;
            const bool cut_sector = sector_1 != sector_2;
            const bool cuts = cut_Q2 && cut_W && cut_zh && cut_y;

            if (cut_Q2 && cut_W && cut_y) m_histograms.event.hist1D_zh_cut->Get()->Fill(zh);
            if (cut_Q2 && cut_zh && cut_y) m_histograms.event.hist1D_W_cut->Get()->Fill(W);
            if (cut_W && cut_zh && cut_y) m_histograms.event.hist1D_Q2_cut->Get()->Fill(Q2);
            if (cut_Q2 && cut_W && cut_zh) m_histograms.event.hist1D_y_cut->Get()->Fill(nu);

            if (cuts) {
                m_histograms.event.hist1D_nu_cut->Get()->Fill(nu);
                m_histograms.event.hist1D_invariant_mass_cut->Get()->Fill(m);
                m_histograms.event.hist2D_angle_vs_invariant_mass->Get()->Fill(ROOT::Math::VectorUtil::Angle(photon1.PxPyPzEVector(), photon2.PxPyPzEVector()) * 180. / TMath::Pi(), m);

                set_photons_pass_event_cuts.insert(photon1);
                set_photons_pass_event_cuts.insert(photon2);
            }
        }

        // Generate pair of photons between good_photon and photon_mix_queue
        std::vector<std::pair<Core::Particle, Core::Particle>> pair_photon_mix;
        for (Core::Particle& p : good_photons) {
            for (Core::Particle& q : photon_mix_queue) {
                std::pair<Core::Particle, Core::Particle> pair = {p, q};
                pair_photon_mix.push_back(pair);
            }
        }

        for (const auto& [photon1, photon2] : pair_photon_mix) {
            // Cuts ---------------------------------------------------------------------------------------
            const double Q2_min = m_config["event"]["Q2_min"].value_or(NaN);
            const double W_min = m_config["event"]["W_min"].value_or(NaN);
            const double zh_min = m_config["event"]["zh_min"].value_or(NaN);
            const double zh_max = m_config["event"]["zh_max"].value_or(NaN);
            const double y_max = m_config["event"]["y_max"].value_or(NaN);
            //---------------------------------------------------------------------------------------------

            ROOT::Math::PxPyPzEVector k1 = {0, 0, 10.53, 10.53};
            ROOT::Math::PxPyPzEVector q1 = k1 - electron.PxPyPzEVector();
            ROOT::Math::PxPyPzEVector p1 = {0, 0, 0, Core::Constantes::ProtonMass};
            ROOT::Math::PxPyPzEVector h = photon1.PxPyPzEVector() + photon2.PxPyPzEVector();

            double M = Core::Constantes::ProtonMass;
            double nu = k1.E() - electron.PxPyPzEVector().E();
            double y = nu / k1.E();
            double Q2 = 2. * q1.E() * electron.E() * (1 - std::cos(electron.PxPyPzEVector().Theta()));
            double W = std::sqrt(M * M + 2 * M * nu - Q2);
            double zh = h.E() / nu;
            double xB = Q2 / (2. * M * nu);
            double pT = h.Pt();
            double m = ROOT::Math::VectorUtil::InvariantMass(photon1.PxPyPzEVector(), photon2.PxPyPzEVector());

            Core::CalorimeterBank calo_info_1 = Core::read_Calorimeter_bank(REC_Calorimeter, photon1.index());
            int sector_1 = calo_info_1.pcal.sector;
            Core::CalorimeterBank calo_info_2 = Core::read_Calorimeter_bank(REC_Calorimeter, photon2.index());
            int sector_2 = calo_info_2.pcal.sector;

            m_histograms.event_mix.hist1D_W->Get()->Fill(W);
            m_histograms.event_mix.hist1D_Q2->Get()->Fill(Q2);
            m_histograms.event_mix.hist1D_nu->Get()->Fill(nu);
            m_histograms.event_mix.hist1D_zh->Get()->Fill(zh);
            m_histograms.event_mix.hist1D_y->Get()->Fill(y);
            m_histograms.event_mix.hist1D_invariant_mass->Get()->Fill(m);

            const bool cut_Q2 = Q2 > Q2_min;
            const bool cut_W = W > W_min;
            const bool cut_zh = zh_min < zh && zh < zh_max;
            const bool cut_y = y < y_max;
            const bool cut_sector = sector_1 != sector_2;
            const bool cuts = cut_Q2 && cut_W && cut_zh && cut_y;

            // if (ROOT::Math::VectorUtil::Angle(photon1.PxPyPzEVector(), photon2.PxPyPzEVector()) * 180. / TMath::Pi() < 4) continue;

            if (cut_Q2 && cut_W && cut_y) m_histograms.event_mix.hist1D_zh_cut->Get()->Fill(zh);
            if (cut_Q2 && cut_zh && cut_y) m_histograms.event_mix.hist1D_W_cut->Get()->Fill(W);
            if (cut_W && cut_zh && cut_y) m_histograms.event_mix.hist1D_Q2_cut->Get()->Fill(Q2);
            if (cut_Q2 && cut_W && cut_zh) m_histograms.event_mix.hist1D_y_cut->Get()->Fill(nu);

            if (cuts) {
                m_histograms.event_mix.hist1D_nu_cut->Get()->Fill(nu);
                m_histograms.event_mix.hist1D_invariant_mass_cut->Get()->Fill(m);

                m_histograms.event_mix.hist2D_angle_vs_invariant_mass->Get()->Fill(ROOT::Math::VectorUtil::Angle(photon1.PxPyPzEVector(), photon2.PxPyPzEVector()) * 180. / TMath::Pi(), m);

                set_photons_pass_event_cuts.insert(photon1);
                set_photons_pass_event_cuts.insert(photon2);
            }
        }

        // Loop over the pairs of photons

        if (mode_write) {
            hipo::bank partBank(schemaPart, 1 + set_photons_pass_event_cuts.size());

            if (!(set_photons_pass_event_cuts.empty())) {
                // Write events with good photons and electron
                hipo::event outEvent = hipo_event;
                outEvent.reset();

                for (const std::string& s : dst_schema) {
                    auto bank = hipo::bank(dict.getSchema(s.c_str()));
                    hipo_event.getStructure(bank);
                    if (bank.getRows() > 0) {
                        outEvent.addStructure(bank);
                    }
                }

                partBank.putInt("pid", 0, 11);
                partBank.putFloat("px", 0, electron.px());
                partBank.putFloat("py", 0, electron.py());
                partBank.putFloat("pz", 0, electron.pz());

                for (int j = 1; const auto& photon : set_photons_pass_event_cuts) {
                    partBank.putInt("pid", j, 22);
                    partBank.putFloat("px", j, photon.px());
                    partBank.putFloat("py", j, photon.py());
                    partBank.putFloat("pz", j, photon.pz());
                    j++;
                }

                outEvent.addStructure(partBank);

                writer.addEvent(outEvent);
            }
        }
    }

    // Close the writer
    if (mode_write) {
        writer.close();
        writer.showSummary();
    }
}

auto Reader::get_topology(const hipo::bank& REC_Particle) const -> Topology {
    std::vector<Core::Particle> electrons;
    std::vector<Core::Particle> photons;

    for (int i = 0; i < REC_Particle.getRows(); i++) {
        int pid = REC_Particle.get<int>("pid", i);
        int status = REC_Particle.get<int>("status", i);
        int charge = REC_Particle.get<int>("charge", i);

        double px = REC_Particle.get<double>("px", i);
        double py = REC_Particle.get<double>("py", i);
        double pz = REC_Particle.get<double>("pz", i);
        double vx = REC_Particle.get<double>("vx", i);
        double vy = REC_Particle.get<double>("vy", i);
        double vz = REC_Particle.get<double>("vz", i);
        double vt = REC_Particle.get<double>("vt", i);
        double beta = REC_Particle.get<double>("beta", i);
        double chi2pid = REC_Particle.get<double>("chi2pid", i);
        double E = Core::compute_energy(px, py, pz, pid);

        if (pid == 11) electrons.emplace_back(pid, status, i, charge, Core::Constantes::ElectronMass, px, py, pz, E, vx, vy, vz, vt, beta, chi2pid);
        if (pid == 22 && E > 0.0001) photons.emplace_back(pid, status, i, charge, 0, px, py, pz, E, vx, vy, vz, vt, beta, chi2pid);
    }

    return {.electrons = electrons, .photons = photons};
}

auto Reader::select_electron(const Core::Particle& electron, const hipo::bank& REC_Calorimeter, const hipo::bank& REC_Cherenkov) const -> bool {
    bool pass = false;

    Core::CherenkovBank cherenkov = Core::read_Cherenkov_bank(REC_Cherenkov, electron.index());
    Core::CalorimeterBank calorimeter = Core::read_Calorimeter_bank(REC_Calorimeter, electron.index());

    double ETot = calorimeter.pcal.energy + calorimeter.inner.energy + calorimeter.outer.energy;
    double SF = ETot / electron.p();

    // Cuts ---------------------------------------------------------------------------------------
    const double vz_min_e = m_config["electron"]["vz_min"].value_or(NaN);
    const double vz_max_e = m_config["electron"]["vz_max"].value_or(NaN);
    const double chi2_min_e = m_config["electron"]["chi2_min"].value_or(NaN);
    const double chi2_max_e = m_config["electron"]["chi2_max"].value_or(NaN);
    const double p_min = m_config["electron"]["p_min"].value_or(NaN);
    //---------------------------------------------------------------------------------------------

    // Get variables for electrons kinematics to apply the cuts -----------------------------------
    const double p_e = electron.p();
    const double chi2_e = electron.chi2pid();
    const double vz_e = electron.vz();
    const double phi_e = electron.phi();
    const double theta_e = electron.theta();
    //---------------------------------------------------------------------------------------------

    // Cuts ---------------------------------------------------------------------------------------
    const bool cut_p = p_e > p_min;
    const bool cut_chi2 = chi2_min_e < chi2_e && chi2_e < chi2_max_e;
    const bool cut_vz = vz_min_e <= vz_e && vz_e <= vz_max_e;
    const bool cuts = cut_chi2 && cut_vz && cut_p;
    // --------------------------------------------------------------------------------------------

    // Fill the histograms before any cuts are applied --------------------------------------------
    m_histograms.electron.hist1D_p->Get()->Fill(p_e);
    m_histograms.electron.hist1D_phi->Get()->Fill(phi_e);
    m_histograms.electron.hist1D_theta->Get()->Fill(theta_e);
    m_histograms.electron.hist1D_chi2->Get()->Fill(chi2_e);
    m_histograms.electron.hist1D_vz->Get()->Fill(vz_e);

    //---------------------------------------------------------------------------------------------

    // Fill histograms after cuts: ----------------------------------------------------------------
    if (cut_chi2 && cut_p) m_histograms.electron.hist1D_vz_cut->Get()->Fill(vz_e);
    if (cut_vz && cut_p) m_histograms.electron.hist1D_chi2_cut->Get()->Fill(chi2_e);
    if (cut_chi2 && cut_vz) m_histograms.electron.hist1D_p_cut->Get()->Fill(p_e);

    if (cuts) {
        m_histograms.electron.hist1D_phi_cut->Get()->Fill(phi_e);
        m_histograms.electron.hist1D_theta_cut->Get()->Fill(theta_e);
        pass = true;
    }
    //---------------------------------------------------------------------------------------------

    return pass;
}

auto Reader::select_photons(const std::vector<Core::Particle>& photons, const Core::Particle& electron, const hipo::bank& REC_Calorimeter, double start_time) const -> std::vector<Core::Particle> {
    // Cuts values --------------------------------------------------------------------------------
    const double E_min = m_config["photon"]["E_min"].value_or(NaN);
    const double time_min = m_config["photon"]["time_min"].value_or(NaN);
    const double time_max = m_config["photon"]["time_max"].value_or(NaN);
    const double opening_angle_min = m_config["photon"]["opening_angle_min"].value_or(NaN);
    //---------------------------------------------------------------------------------------------

    std::vector<Core::Particle> selected_photons;

    for (const auto& photon : photons) {
        Core::CalorimeterBank calo_info = Core::read_Calorimeter_bank(REC_Calorimeter, photon.index());
        double time = calo_info.pcal.time - start_time - calo_info.pcal.path / 30;

        double opening_angle = electron.theta(photon);

        // Cuts ---------------------------------------------------------------------------------------
        const bool cut_E = photon.E() > E_min;
        const bool cut_time = time_min < time && time < time_max;
        const bool cut_opening_angle = opening_angle > opening_angle_min;
        const bool cuts = cut_E && cut_time && cut_opening_angle;
        // --------------------------------------------------------------------------------------------

        m_histograms.photon.hist1D_E->Get()->Fill(photon.E());
        m_histograms.photon.hist1D_beta->Get()->Fill(photon.beta());
        m_histograms.photon.hist1D_time->Get()->Fill(time);
        m_histograms.photon.hist1D_vz->Get()->Fill(photon.vz());
        m_histograms.photon.hist1D_opening_angle->Get()->Fill(opening_angle);
        m_histograms.photon.hist1D_delta_vz->Get()->Fill(electron.vz() - photon.vz());

        m_histograms.photon.hist2D_phi_vs_theta->Get()->Fill(photon.phi(), photon.theta());
        m_histograms.photon.hist2D_E_vs_vz->Get()->Fill(photon.E(), photon.vz());

        if (cut_E && cut_opening_angle) m_histograms.photon.hist1D_time_cut->Get()->Fill(time);
        if (cut_time && cut_opening_angle) m_histograms.photon.hist1D_E_cut->Get()->Fill(photon.E());
        if (cut_E && cut_time) m_histograms.photon.hist1D_opening_angle_cut->Get()->Fill(opening_angle);

        if (cuts) {

            m_histograms.photon.hist1D_vz_cut->Get()->Fill(photon.vz());
            m_histograms.photon.hist1D_beta_cut->Get()->Fill(photon.beta());

            m_histograms.photon.hist2D_phi_vs_theta_cut->Get()->Fill(photon.phi(), photon.theta());
            m_histograms.photon.hist2D_E_vs_vz_cut->Get()->Fill(photon.E(), photon.vz());

            m_histograms.photon.hist1D_delta_vz_cut->Get()->Fill(electron.vz() - photon.vz());

            selected_photons.push_back(photon);
        }
    }

    return selected_photons;
}

// Stuff for GBTFilter
auto Reader::select_photons_with_AI(const hipo::bank& RUN_config, const hipo::bank& REC_Particle, const hipo::bank& REC_Calorimeter) const -> std::vector<Core::Particle> {
    // ---------------------------
    //     Stuff for GBTFilter
    // ---------------------------
    int run_number = RUN_config.getInt("run", 0);
    std::vector<int> good_photon_rows;
    // Get CaloMap for the event
    auto calo_map = GetCaloMap(REC_Calorimeter);
    for (int i = 0; i < REC_Particle.getRows(); i++) {
        int pid = REC_Particle.getInt("pid", i);
        if (pid == 22) {
            if (Filter(REC_Particle, REC_Calorimeter, calo_map, i, run_number)) {
                good_photon_rows.push_back(i);
            }
        }
    }

    // Get good photons
    std::vector<Core::Particle> good_photons;
    for (const auto& row : good_photon_rows) {
        int pid = REC_Particle.getInt("pid", row);
        int status = REC_Particle.getInt("status", row);
        int charge = REC_Particle.getInt("charge", row);

        double px = REC_Particle.getDouble("px", row);
        double py = REC_Particle.getDouble("py", row);
        double pz = REC_Particle.getDouble("pz", row);
        double vx = REC_Particle.getDouble("vx", row);
        double vy = REC_Particle.getDouble("vy", row);
        double vz = REC_Particle.getDouble("vz", row);
        double vt = REC_Particle.getDouble("vt", row);
        double beta = REC_Particle.getDouble("beta", row);
        double chi2pid = REC_Particle.getDouble("chi2pid", row);
        double E = Core::compute_energy(px, py, pz, pid);

        good_photons.emplace_back(pid, status, row, charge, 0, px, py, pz, E, vx, vy, vz, vt, beta, chi2pid);
    }

    return good_photons;
    // End ---------------------------------------
}

std::map<int, Reader::calo_row_data> Reader::GetCaloMap(hipo::bank const& bank) const {
    std::map<int, Reader::calo_row_data> calo_map;
    // Loop over REC::Calorimeter rows
    // Here we use bank.getRows() to purposefully ignore upstream filters
    for (int row = 0; row < bank.getRows(); row++) {
        auto pindex = bank.getShort("pindex", row);
        auto x = bank.getFloat("x", row);
        auto y = bank.getFloat("y", row);
        auto z = bank.getFloat("z", row);
        auto m2u = bank.getFloat("m2u", row);
        auto m2v = bank.getFloat("m2v", row);
        auto layer = bank.getInt("layer", row);
        auto e = bank.getFloat("energy", row);

        // Ensure an entry exists in the map for the given pindex
        if (calo_map.find(pindex) == calo_map.end()) {
            calo_map[pindex] = Reader::calo_row_data();
        }

        switch (layer) {
            case 1:  // pcal
                calo_map[pindex].pcal_x = x;
                calo_map[pindex].pcal_y = y;
                calo_map[pindex].pcal_z = z;
                calo_map[pindex].pcal_e = e;
                calo_map[pindex].pcal_m2u = m2u;
                calo_map[pindex].pcal_m2v = m2v;
                break;
            case 4:  // ecin
                calo_map[pindex].ecin_x = x;
                calo_map[pindex].ecin_y = y;
                calo_map[pindex].ecin_z = z;
                break;
            case 7:  // ecout
                calo_map[pindex].ecout_x = x;
                calo_map[pindex].ecout_y = y;
                calo_map[pindex].ecout_z = z;
                break;
        }
    }
    return calo_map;
}

bool Reader::Filter(hipo::bank const& particleBank, hipo::bank const& caloBank, std::map<int, Reader::calo_row_data> calo_map, int const row, int const runnum) const {

    // Set variables native to the photon we are classifying
    double gPx = particleBank.getFloat("px", row);
    double gPy = particleBank.getFloat("py", row);
    double gPz = particleBank.getFloat("pz", row);

    // Set ML features intrinsic to the photon of interest
    double gE = sqrt(gPx * gPx + gPy * gPy + gPz * gPz);
    double gTheta = acos(gPz / gE);
    double gEpcal = calo_map[row].pcal_e;
    double gm2u = calo_map[row].pcal_m2u;
    double gm2v = calo_map[row].pcal_m2v;

    // Apply PID purity cuts on the photon
    // If they do not pass, then these photons are incompatible with the trained GBT model
    if (PidPurityPhotonFilter(gE, gEpcal, gTheta) == false) return false;

    //Define the variables "m_g" , "m_ch" , "m_nh"
    //Should not be changed because the model was trained with this specific set of inputs
    int const m_g = 3;   // Number of neighboring gammas
    int const m_ch = 2;  // Number of neighboring charged hadrons (protons, pions, kaons)
    int const m_nh = 2;  // Number of neighboring neutral hadrons (neutrons)

    double R_e = 0;
    double dE_e = 0;

    double R_gamma[m_g];      // Angular distance between calo shower centers
    double dE_gamma[m_g];     // Energy difference
    double Epcal_gamma[m_g];  // Energy deposited in the pcal
    double m2u_gamma[m_g];    // Shower shape variables
    double m2v_gamma[m_g];    // Shower shape variables

    double R_ch[m_ch];      // Angular distance between calo shower centers
    double dE_ch[m_ch];     // Energy difference
    double Epcal_ch[m_ch];  // Energy deposited in the pcal
    double m2u_ch[m_ch];    // Shower shape variables
    double m2v_ch[m_ch];    // Shower shape variables

    double R_nh[m_nh];      // Angular distance between calo shower centers
    double dE_nh[m_nh];     // Energy difference
    double Epcal_nh[m_nh];  // Energy deposited in the pcal
    double m2u_nh[m_nh];    // Shower shape variables
    double m2v_nh[m_nh];    // Shower shape variables

    double num_photons_0_1, num_photons_0_2, num_photons_0_35;

    //Initialize the arrays
    for (int i = 0; i < m_g; ++i) {
        R_gamma[i] = 0;
        dE_gamma[i] = 0;
        Epcal_gamma[i] = 0;
        m2u_gamma[i] = 0;
        m2v_gamma[i] = 0;
    }
    for (int i = 0; i < m_ch; ++i) {
        R_ch[i] = 0;
        dE_ch[i] = 0;
        Epcal_ch[i] = 0;
        m2u_ch[i] = 0;
        m2v_ch[i] = 0;
    }
    for (int i = 0; i < m_nh; ++i) {
        R_nh[i] = 0;
        dE_nh[i] = 0;
        Epcal_nh[i] = 0;
        m2u_nh[i] = 0;
        m2v_nh[i] = 0;
    }

    //Set the number of photons within R<0.1, R<0.2, R<0.35
    num_photons_0_1 = 0;
    num_photons_0_2 = 0;
    num_photons_0_35 = 0;

    // Get 3-vector that points to the photon of interest's location in the calorimeter
    auto calo_POI = calo_map.at(row);
    ROOT::Math::XYZVector vPOI = GetParticleCaloVector(calo_POI);

    // Build nearest neighbor event structure
    // Loop over particles in the event
    // Here we loop over particleBank.getRows(), which ignores upstream filters
    // This is critical as the GBTs were trained on identifying nearest neighbors for the whole REC::Particle bank
    // Only considering nearest neighbor particles that pass upstream filters would call the accuracy of the model into question
    for (int inner_row = 0; inner_row < particleBank.getRows(); inner_row++) {
        // Skip over the particle if it is photon we are trying to classify
        if (inner_row == row) continue;

        // Check if the key exists in the calo_map
        // This skips over REC::Particle entries without a REC::Calorimeter entry
        if (calo_map.find(inner_row) == calo_map.end()) continue;
        auto calo_PART = calo_map.at(inner_row);

        auto pid = particleBank.getInt("pid", inner_row);
        auto mass = Core::get_mass(pid);

        // Skip over particle if its mass was undefined
        if (mass == std::numeric_limits<double>::quiet_NaN()) continue;
        auto px = particleBank.getFloat("px", inner_row);
        auto py = particleBank.getFloat("py", inner_row);
        auto pz = particleBank.getFloat("pz", inner_row);
        auto p = sqrt(px * px + py * py + pz * pz);
        auto E = sqrt(p * p + mass * mass);
        auto th = acos(pz / p);
        // Skip over particle if it is not in the forward detector (necessary for model compatibility)
        if (ForwardDetectorFilter(th) == false) continue;

        // Get 3-vector that points to the neighboring particle's location in the calorimeter
        ROOT::Math::XYZVector vPART = GetParticleCaloVector(calo_PART);

        // Get angular distance between photon of interest and particle
        double R = ROOT::Math::VectorUtil::Angle(vPOI, vPART);

        // Logic for filling nearest neighbor variables
        if (pid == 22) {  //photon

            // Apply Photon Purity Cuts to ensure this neighbor can be used in classification
            if (PidPurityPhotonFilter(E, calo_PART.pcal_e, th) == false) continue;

            if (R < 0.1) num_photons_0_1++;
            if (R < 0.2) num_photons_0_2++;
            if (R < 0.35) num_photons_0_35++;

            for (int i = 0; i < m_g; ++i) {
                if (R < R_gamma[i] || R_gamma[i] == 0) {
                    int j = m_g - 1;
                    while (j > i) {
                        R_gamma[j] = R_gamma[j - 1];
                        dE_gamma[j] = dE_gamma[j - 1];
                        Epcal_gamma[j] = Epcal_gamma[j - 1];
                        m2u_gamma[j] = m2u_gamma[j - 1];
                        m2v_gamma[j] = m2v_gamma[j - 1];
                        j--;
                    }
                    R_gamma[i] = R;
                    dE_gamma[i] = gE - E;
                    Epcal_gamma[i] = calo_PART.pcal_e;
                    m2u_gamma[i] = calo_PART.pcal_m2u;
                    m2v_gamma[i] = calo_PART.pcal_m2v;
                    break;
                }
            }
        } else if (pid == 11) {  //electron
            if (R < R_e || R_e == 0) {
                R_e = R;
                dE_e = gE - E;
            }
        } else if (pid == 2212 || pid == -2212 || pid == 211 || pid == -211 || pid == 321 || pid == -321) {  //charged hadron
            for (int i = 0; i < m_ch; ++i) {
                if (R < R_ch[i] || R_ch[i] == 0) {
                    int j = m_ch - 1;
                    while (j > i) {
                        R_ch[j] = R_ch[j - 1];
                        dE_ch[j] = dE_ch[j - 1];
                        Epcal_ch[j] = Epcal_ch[j - 1];
                        m2u_ch[j] = m2u_ch[j - 1];
                        m2v_ch[j] = m2v_ch[j - 1];
                        j--;
                    }
                    R_ch[i] = R;
                    dE_ch[i] = gE - E;
                    Epcal_ch[i] = calo_PART.pcal_e;
                    m2u_ch[i] = calo_PART.pcal_m2u;
                    m2v_ch[i] = calo_PART.pcal_m2v;
                    break;
                }
            }
        } else if (pid == 2112 || pid == -2112) {  //neutral hadron
            for (int i = 0; i < m_nh; ++i) {
                if (R < R_nh[i] || R_nh[i] == 0) {
                    int j = m_nh - 1;
                    while (j > i) {
                        R_nh[j] = R_nh[j - 1];
                        dE_nh[j] = dE_nh[j - 1];
                        Epcal_nh[j] = Epcal_nh[j - 1];
                        m2u_nh[j] = m2u_nh[j - 1];
                        m2v_nh[j] = m2v_nh[j - 1];
                        j--;
                    }
                    R_nh[i] = R;
                    dE_nh[i] = gE - E;
                    Epcal_nh[i] = calo_PART.pcal_e;
                    m2u_nh[i] = calo_PART.pcal_m2u;
                    m2v_nh[i] = calo_PART.pcal_m2v;
                    break;
                }
            }
        } else {  //unrecognized OR uncompatible particle type for the trained model
            continue;
        }
    }

    // Create and populate input_data vector for the ML model
    std::vector<float> input_data = {
        static_cast<float>(gE), static_cast<float>(gEpcal), static_cast<float>(gTheta),
        static_cast<float>(gm2u), static_cast<float>(gm2v), static_cast<float>(R_e),
        static_cast<float>(dE_e)};

    for (int i = 0; i < m_g; ++i)
        input_data.push_back(static_cast<float>(R_gamma[i]));
    for (int i = 0; i < m_g; ++i)
        input_data.push_back(static_cast<float>(dE_gamma[i]));
    for (int i = 0; i < m_g; ++i)
        input_data.push_back(static_cast<float>(Epcal_gamma[i]));
    for (int i = 0; i < m_g; ++i)
        input_data.push_back(static_cast<float>(m2u_gamma[i]));
    for (int i = 0; i < m_g; ++i)
        input_data.push_back(static_cast<float>(m2v_gamma[i]));

    for (int i = 0; i < m_ch; ++i)
        input_data.push_back(static_cast<float>(R_ch[i]));
    for (int i = 0; i < m_ch; ++i)
        input_data.push_back(static_cast<float>(dE_ch[i]));
    for (int i = 0; i < m_ch; ++i)
        input_data.push_back(static_cast<float>(Epcal_ch[i]));
    for (int i = 0; i < m_ch; ++i)
        input_data.push_back(static_cast<float>(m2u_ch[i]));
    for (int i = 0; i < m_ch; ++i)
        input_data.push_back(static_cast<float>(m2v_ch[i]));

    for (int i = 0; i < m_nh; ++i)
        input_data.push_back(static_cast<float>(R_nh[i]));
    for (int i = 0; i < m_nh; ++i)
        input_data.push_back(static_cast<float>(dE_nh[i]));
    for (int i = 0; i < m_nh; ++i)
        input_data.push_back(static_cast<float>(Epcal_nh[i]));
    for (int i = 0; i < m_nh; ++i)
        input_data.push_back(static_cast<float>(m2u_nh[i]));
    for (int i = 0; i < m_nh; ++i)
        input_data.push_back(static_cast<float>(m2v_nh[i]));

    input_data.push_back(static_cast<float>(num_photons_0_1));
    input_data.push_back(static_cast<float>(num_photons_0_2));
    input_data.push_back(static_cast<float>(num_photons_0_35));

    return ClassifyPhoton(input_data, runnum);
}

bool Reader::PidPurityPhotonFilter(float const E, float const Epcal, float const theta) const {
    // Apply standard pid cuts on the photon, compatible to how the models were trained
    // 1. Minimum photon energy cut of 200 MeV
    // 2. Photon must have deposited energy in the PCal
    // 3. Photon must be in the Forward Detector
    if (E < 0.2 || Epcal <= 0 || !(ForwardDetectorFilter(theta))) return false;
    return true;
}

bool Reader::ForwardDetectorFilter(float const theta) const {
    // Apply forward detector cut
    if (theta * 180 / M_PI < 5 || theta * 180 / M_PI > 35) return false;
    return true;
}

ROOT::Math::XYZVector Reader::GetParticleCaloVector(Reader::calo_row_data calo_row) const {
    // Determine the 3-vector location of where the photon of interest's calo deposition is
    // First we check the pcal coords, then ecin, then ecout
    ROOT::Math::XYZVector v;
    if (calo_row.pcal_x == 0) {
        if (calo_row.ecin_x == 0) {
            v.SetXYZ(calo_row.ecout_x, calo_row.ecout_y, calo_row.ecout_z);
        } else {
            v.SetXYZ(calo_row.ecin_x, calo_row.ecin_y, calo_row.ecin_z);
        }
    } else {
        v.SetXYZ(calo_row.pcal_x, calo_row.pcal_y, calo_row.pcal_z);
    }
    return v;
}

bool Reader::ClassifyPhoton(std::vector<float> const& input_data, int const runnum) const {
    auto modelFunction = getModelFunction(runnum);
    double sigmoid_x = modelFunction(input_data);
    double prediction = 1 / (1 + exp(-sigmoid_x));
    return (prediction > o_threshold);
}

std::function<double(std::vector<float> const&)> Reader::getModelFunction(int runnum) const {

    for (const auto& entry : modelMap) {
        if (runnum >= std::get<0>(entry.first) && runnum <= std::get<1>(entry.first) && o_pass == std::get<2>(entry.first)) {
            return entry.second;
        }
    }

    return [](std::vector<float> const& data) {
        return ApplyCatboostModel_RGA_inbending_pass1(data);
    };
}

}  // namespace study1