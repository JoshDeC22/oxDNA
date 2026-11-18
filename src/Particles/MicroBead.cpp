/*
 * MicroBead.cpp
 *
 * Created on: 17/nov/2026
 *  Author: Joshua DeCoteau
 */

#include "MicroBead.h"
#include "../Utilities/Utils.h"
#include "../Utilities/oxDNAException.h"

void MicroBead::add_site(const LR_vector &position, unsigned int target_nucleotide) {
    LR_vector local_position = orientationT * (position - pos);

    // Check whether the tether site lies on the surface of the micro bead.
    if (local_position.module() - radius > tolerance) {
        throw oxDNAException("Error in micro bead %zu: the given tether position (%f, %f, %f) does not lie on the surface of the micro bead.", id, position.x, position.y, position.z);
    }

    InteractionSite site;
    site.local_position = local_position;
    site.target_nucleotide = target_nucleotide;

    sites.emplace_back(site);
}

void MicroBead::calc_rot_axes() {
    if (constrained) {
        for (size_t i = 0; i < n_midpoints; i += 2) {
            LR_vector midpoint = (sites[i].local_position + sites[i + 1].local_position);
            midpoint = midpoint / midpoint.module();

            rot_axes.emplace_back(midpoint);
        }

        for (size_t i = 2 * n_midpoints; i < sites.size(); i++) {
            LR_vector axis = sites[i].local_position / sites[i].local_position.module();
            rot_axes.emplace_back(axis);
        }
    }
}

MicroBead::MicroBead(input_file &inp) : BaseParticle() {
    getInputNumber(&inp, "radius", &radius, 1);
    getInputNumber(&inp, "anchor_stiffness", &anchor_stiffness, 1);
    getInputNumber(&inp, "rotation_stiffness", &rot_stiffness, 1);

    getInputUInt(&inp, "id", &id, 1);

    if (getInputBool(&inp, "constrained", &constrained, 0) == KEY_NOT_FOUND) {
        constrained = false;
    }

    if (getInputUInt(&inp, "num_midpoints", &n_midpoints, 0) == KEY_NOT_FOUND) {
        n_midpoints = 0;
    }

    std::string particles_str;
    getInputString(&inp, "particle", particles_str, 1);

    std::vector<unsigned int> particles = Utils::split_to_uint(particles_str, ",");

    std::string anchor_str;
    getInputString(&inp, "anchor_positions", anchor_str, 1);

    anchor_str.erase(
        std::remove_if(anchor_str.begin(), anchor_str.end(), [](unsigned char c) {
            return std::isspace(c) || c == '(' || c == ')';
        }),
        anchor_str.end()
    );

    std::vector<number> anchors = Utils::split_to_numbers(anchor_str, ",");

    if (anchors.size() % 3 != 0) {
        throw oxDNAException("Error in micro bead %zu: the anchor positions should all have 3 values.", id);
    }

    if (particles.size() != anchors.size() / 3) {
        throw oxDNAException("Error in micro bead %zu: found a different number of anchors and anchored particles. Found %zu anchored particles and %zu anchors.", id, particles.size(), anchors.size() / 3);
    }

    if (n_midpoints > particles.size() / 2) {
        throw oxDNAException("Error in micro bead %zu: too many midpoints. Expected a maximum of %zu, got %zu", particles.size() / 2, n_midpoints);
    }

    size_t ind = 0;
    for (unsigned int particle : particles) {
        LR_vector anchor_position = LR_vector(anchors[ind], anchors[ind + 1], anchors[ind + 2]);
        add_site(anchor_position, particle);

        ind += 3;
    }

    int_centers.resize(sites.size());

    calc_rot_axes();
}

void MicroBead::set_positions() {
    for (size_t i = 0; i < sites.size(); i++) {
        int_centers[i] = orientation * sites[i].local_position;
    }
}

std::vector<LR_vector> MicroBead::get_rot_axes() const {
    std::vector<LR_vector> global_rot;

    if (!constrained) {
        return global_rot;
    }

    global_rot.resize(rot_axes.size());

    for (size_t i = 0; i < rot_axes.size(); i++) {
        global_rot[i] = orientation * rot_axes[i];
    }

    return global_rot;
}
