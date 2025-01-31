#pragma once

#include <vector>
#include <memory>
#include <utility>
#include <tuple>
#include <string>

#include "chargefw2.h"
#include "molecule.h"
#include "../method.h"
#include "../parameters.h"


enum class AtomClassifier {
    PLAIN,
    HBO,
    BONDED
};


enum class BondClassifier {
    PLAIN,
    BO
};


template<typename T>
struct deduce_from
{
};

template<>
struct deduce_from<Atom>
{
    typedef atom_t AB_t;
};

template<>
struct deduce_from<Bond>
{
    typedef bond_t AB_t;
};

struct MoleculeSetStats {
    struct AtomTypeCount {
        std::string symbol;
        std::string cls;
        std::string type;
        int count;
    };

    size_t total_molecules;
    size_t total_atoms;
    std::vector<AtomTypeCount> atom_type_counts;
};

class MoleculeSet {
    std::vector<atom_t> atom_types_{};
    std::vector<bond_t> bond_types_{};
    std::unique_ptr<std::vector<Molecule> > molecules_{nullptr};

    template<typename AB, typename AB_t = typename deduce_from<AB>::AB_t>
    void set_type(AB &object, const AB_t &type);

    template<typename AB, typename AB_t = typename deduce_from<AB>::AB_t>
    size_t classify_objects_from_parameters(const Parameters &parameters,
                                            bool remove_unclassified,
                                            bool permissive_types);

public:
    explicit MoleculeSet() = default;

    explicit MoleculeSet(std::unique_ptr<std::vector<Molecule> > molecules);

    void info() const;

    MoleculeSetStats get_stats() const;

    [[nodiscard]] const std::vector<Molecule> &molecules() const { return *molecules_; }

    void classify_atoms(AtomClassifier cls);

    void classify_bonds(BondClassifier cls);

    size_t classify_set_from_parameters(const Parameters &parameters, bool remove_unclassified = true,
                                        bool permissive_types = false);

    [[nodiscard]] std::vector<atom_t> atom_types() const { return atom_types_; }

    [[nodiscard]] std::vector<bond_t> bond_types() const { return bond_types_; }

    void fulfill_requirements(const std::vector<RequiredFeatures> &features);

    [[nodiscard]] bool has_proteins() const;
};
