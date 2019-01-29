//
// Created by krab1k on 11.12.18.
//

#include  <stdexcept>

#include "element.h"


int Element::valence_electron_count() const {
    if (group_ < 3) {
        return group_;
    } else if (group_ > 12){
        return group_ - 10;
    }
    throw std::runtime_error("Valence electron count undefined");
}
