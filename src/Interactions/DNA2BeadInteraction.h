/*
* @brief Adds the interactions between nucleotides and MicroBead particles to the existing oxDNA2 implementation
*
* Created on: 17/nov/2026
*   Author: Joshua DeCoteau
*/

#ifndef DNA2BEAD_INTERACTION_H
#define DNA2BEAD_INTERACTION_H

#include "DNA2Interaction.h"
#include "../Particles/BaseParticle.h"
#include "../Particles/MicroBead.h"
#include "DNAInteraction.h"

class DNA2BeadInteraction : public DNA2Interaction {
protected:
    bool include_beads;
    std::string bead_file;

    std::vector<input_file> get_beadblocks();

    number dna_bead_interactions(MicroBead *p, BaseParticle *q, bool update_forces);
public:
    DNA2BeadInteraction() : DNAInteraction() {};
    ~DNA2BeadInteraction() {};

    void get_settings(input_file &inp) override;
    void init() override;
    void read_topology(int *N_strands, std::vector<BaseParticle *> &particles) override;

    void allocate_particles(std::vector<BaseParticle *> &particles) override;

    number pair_interaction_nonbonded(BaseParticle *p, BaseParticle *q, bool compute_r = true, bool update_forces = false) override;
};

#endif
