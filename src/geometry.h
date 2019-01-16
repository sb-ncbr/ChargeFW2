//
// Created by krab1k on 02/01/19.
//

#pragma once

#include "structures/atom.h"
#include "structures/bond.h"


double distance(const Atom &atom1, const Atom &atom2);

double distance(const Atom &atom, const Bond &bond, bool weighted = false);

double distance(const Bond &bond1, const Bond &bond2, bool weighted = false);
