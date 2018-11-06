//
// Created by krab1k on 6.11.18.
//

#include <vector>
#include <iostream>
#include <QString>

#include "Utility.h"


std::ostream &operator<<(std::ostream &os, const std::vector<QString> &vec) {
    for (const auto &val: vec) {
        os << val.toStdString() << " ";
    }
    return os;
}
