#include "DNA2BeadInteraction.h"
#include "BaseInteraction.h"
#include "DNA2Interaction.h"
#include "../Utilities/TopologyParser.h"
#include "../Utilities/oxDNAException.h"
#include "../Utilities/Utils.h"
#include "../Particles/DNANucleotide.h"

#include <fstream>
#include <sstream>

std::vector<input_file> DNA2BeadInteraction::get_beadblocks() {
    std::vector<input_file> bead_blocks;

    if (!include_beads) {
        return bead_blocks;
    }

    std::ifstream beads(bead_file.c_str());
    if (!beads.good()) {
        throw oxDNAException("Can't read bead file '%s'", bead_file.c_str());
    }

    OX_LOG(Logger::LOG_INFO, "Parsing bead file %s", bead_file.c_str());

    int open = 0;
    int justopen = 0;
    int a = beads.get();

    bool is_commented = false;
    std::stringstream bead_str;

    while (beads.good()) {
        justopen = 0;

        switch (a) {
            case '#':
                is_commented = true;
                break;
            case '\n':
                is_commented = false;
                break;
            case '{':
                if (!is_commented) {
                    open++;
                    justopen = 1;
                }
                break;
            case '}':
                if (!is_commented) {
                    if (justopen) throw oxDNAException("Syntax error in '%s': nothing between parentheses", bead_file.c_str());
                    open--;
                }
                break;
            default:
                break;
        }

        if (!is_commented) bead_str << (char) a;
        if (open > 1 || open < 0) throw oxDNAException("Syntax error in '%s': parentheses do not match", bead_file.c_str());
        a = beads.get();
    }

    bead_str.clear();
    bead_str.seekg(0, std::ios::beg);
    a = beads.get();

    while (bead_str.good()) {
        while (a != '{' && bead_str.good()) {
            a = bead_str.get();
        }
        if (!bead_str.good()) {
            break;
        }

        a = bead_str.get();
        std::string input_string("");

        while (a != '}' && bead_str.good()) {
            input_string += a;
            a = bead_str.get();
        }

        input_file inp;
        inp.init_from_string(input_string);

        bead_blocks.emplace_back(inp);
    }

    return bead_blocks;
}

number DNA2BeadInteraction::potential(MicroBead *p, BaseParticle *q) {
    MicroBead::AnchorResult result = p->get_anchor(q->index);

    number energy = 0.0;

    if (result.attachment_exists) {
        LR_vector dr = q->pos - result.anchor;

        energy += 0.5 * p->anchor_stiffness * dr.norm();
    }

    number r = (p->pos - q->pos).module();
    number overlap = p->radius - r;

    if (overlap > 0.0) {
        energy += 0.5 * p->ev_stiffness * overlap * overlap;
    }

    return energy;
}

number DNA2BeadInteraction::dna_bead_interactions(MicroBead *p, BaseParticle *q, bool update_forces) {
    MicroBead::AnchorResult result = p->get_anchor(q->index);

    number energy = 0.0;

    // Tethers
    if (result.attachment_exists) {
        LR_vector dr = q->pos - result.anchor;

        energy += 0.5 * p->anchor_stiffness * dr.norm();

        if (update_forces) {
            LR_vector F = dr * (-p->anchor_stiffness);
            q->force += F;
            p->force -= F;

            LR_vector lever = result.anchor - p->pos;
            LR_vector torque_global = lever.cross(-F);

            p->torque += p->orientationT * torque_global;
        }
    }

    // Excluded volume
    number r = (p->pos - q->pos).module();
}

void DNA2BeadInteraction::get_settings(input_file &inp) {
    DNA2Interaction::get_settings(inp);

    getInputBool(&inp, "beads", &include_beads, 0);

    if (include_beads) {
        getInputString(&inp, "bead_file", bead_file, 1);
    }
}

void DNA2BeadInteraction::init() {
    DNA2Interaction::init();
}

void DNA2BeadInteraction::read_topology(int *N_strands, std::vector<BaseParticle *> &particles) {
    BaseInteraction::read_topology(N_strands, particles);

    TopologyParser parser(_topology_filename);

    int N_from_topology = 0;

    if (parser.is_new_topology()) {
        parser.parse_new_topology();

        int strand_n = 0;
        int current_ind = 0;

        for (input_file &seq_input : parser.lines()) {
            bool is_circular = false;
            getInputBool(&seq_input, "circular", &is_circular, 0);

            std::string type;
            std::string seq;

            getInputString(&seq_input, "specs", seq, 1);

            if (getInputString(&seq_input, "type", type, 0) == KEY_FOUND) {
                if (type != "DNA" && type != "BEAD") {
                    throw oxDNAException("topology file, strand %d (line %d): the DNA2Bead interaction only supports type=DNA or type=BEAD strands.", strand_n, strand_n + 1);
                }
            }

            std::vector<int> btypes;
            try {
                btypes = Utils::btypes_from_bead_seq(seq);
            }
            catch (oxDNAException &e) {
                throw oxDNAException("topology file, strand %d (line %d): %s", strand_n, strand_n + 1, e.what());
            }

            int N_in_strand = btypes.size();
            if (current_ind + N_in_strand > parser.N()) {
                throw oxDNAException("Too many particles found in the topology file (should be %d), aborting", parser.N());
            }

            for (int i = 0; i < N_in_strand; i++, current_ind++) {
                BaseParticle *p = particles[current_ind];

                p->strand_id = strand_n;
                p->btype = btypes[i];
                p->n3 = P_VIRTUAL;
                p->n5 = P_VIRTUAL;

                if (type == "DNA") {
                    p->type = (p->btype < 0) ? 3 - ((3 - p->btype) % 4) : p->btype % 4;

                    if (i > 0) {
                        p->n5 = particles[current_ind - 1];
                        p->affected.push_back(ParticlePair(p->n5, p));
                    }
                    if (i < N_in_strand - 1) {
                        p->n3 = particles[current_ind + 1];
                        p->affected.push_back(ParticlePair(p->n3, p));
                    }
                    else if (is_circular) {
                        BaseParticle *first = particles[current_ind - i];
                        p->n3 = first;
                        p->affected.push_back(ParticlePair(first, p));
                        first->n5 = p;
                        first->affected.push_back(ParticlePair(p, first));
                    }
                }
                else {
                    p->type = B_B;
                }
            }

            strand_n++;
        }
        N_from_topology = current_ind;
    }
    else {
        N_from_topology = parser.parse_old_bead_topology(particles);
    }

    if (N_from_topology != (int) particles.size()) {
        throw oxDNAException("The number of particles specified in the header of the topology file (%d) and the number of particles found in the topology file (%d) don't match. Aborting.", N_from_topology, particles.size());
    }

    *N_strands = parser.N_strands();
}

void DNA2BeadInteraction::allocate_particles(std::vector<BaseParticle *> &particles) {
    std::vector<input_file> bead_blocks = get_beadblocks();

    for (size_t i = 0; i < particles.size() - bead_blocks.size(); i++) {
        particles[i] = new DNANucleotide(_grooving);
    }

    size_t ind = 0;
    for (size_t i = particles.size() - bead_blocks.size(); i < particles.size(); i++, ind++) {
        input_file bead_block = bead_blocks[ind];
        particles[i] = new MicroBead(bead_block);
    }
}
