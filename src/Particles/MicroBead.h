/* 
* MicroBead.h
*
* Created on: 17/nov/2026
*   Author: Joshua DeCoteau
*/

#include "BaseParticle.h"

/*
 * @brief Represents a plastic bead commonly seen in optical tweezers experiments
*/

class MicroBead: public BaseParticle {
protected:
    unsigned int id;
    unsigned int n_midpoints;
    number tolerance = 1e-12;

    std::vector<LR_vector> rot_axes;

    void add_site(const LR_vector &position, unsigned int target_nucleotide);

    void calc_rot_axes();
public:
    struct InteractionSite {
        LR_vector local_position;
        unsigned int target_nucleotide;
    };

    number radius;
    number anchor_stiffness;
    number rot_stiffness;
    bool constrained;

    std::vector<InteractionSite> sites;

    MicroBead(input_file &inp);
    virtual ~MicroBead() = default;

    bool is_rigid_body() const override {
        return true;
    }

    void set_positions() override;

    std::vector<LR_vector> get_rot_axes() const;
};
