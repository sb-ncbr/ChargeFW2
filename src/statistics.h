//
// Created by krab1k on 4.1.19.
//

#pragma once

#include "charges.h"


double RMSD(const Charges &charges1, const Charges &charges2);

double Pearson2(const Charges &charges1, const Charges &charges2);

double D_max(const Charges &charges1, const Charges &charges2);

double D_avg(const Charges &charges1, const Charges &charges2);
