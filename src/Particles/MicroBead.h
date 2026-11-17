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
public:
    MicroBead() = delete;
    explicit MicroBead(number radius);
    ~MicroBead();

    number radius;

    bool is_rigid_body() const override {
        return true;
    }
};
