/**
 * This file manages all pairs of collisions.
 * This means, it performs various approaches
 * for a broad collision detection.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_PAIR_H__
#define __DCA_PAIR_H__

#include "Primitives.h"
#include "Utils.h"

#if BUILD_COMPACT_N_SEARCH

#include <CompactNSearch>
#include <array>  // for windows
#include <exception>

#endif /* BUILD_COMPACT_N_SEARCH */

namespace DCA {
/**
 * This is the base class for all generators.
 */
class PairGenerator {
public:
    /**
     * This function generates the pairs given the primitives.
     * The returned vector consists of pairs, where each pair holds two numbers:
     * The indices of the corresponding primitives which were given.
     * @param primitives All primitives to generate the pairs from.
     * @return A vector of pairs of indices, where each index corresponds to a primitive in the primitives vector.
     */
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives) const = 0;
};

/**
 * This generator creates all possible permutations of pairs.
 * This means, the amount of pairs created is n^2, where n = #primitives.
 */
class PermutationPairGenerator : public PairGenerator {
public:
    /**
     * @copydoc PairGenerator::generate
     */
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives) const override;
};

#if BUILD_COMPACT_N_SEARCH
/**
 * This generator uses the CompactNSearch library 
 * for efficient neighbor searching.
 * Altough the neighborhood is computed new for each call to generate,
 * it serves as a first usefull implementation
 */
class NeighborsPairGenerator : public PairGenerator {
public:
    /**
     * Construct this generator with a given radius
     * @param radius The radius to search other primitives in.
     */
    NeighborsPairGenerator(const double &radius) : m_radius(radius) {}

    /**
     * @copydoc PairGenerator::generate
     */
    virtual std::vector<pair_t> generate(
        const std::vector<primitive_t> &primitives) const override;

private:
    double m_radius;
};
#endif /* BUILD_COMPACT_N_SEARCH */

}  // namespace DCA
#endif /* __DCA_PAIR_H__ */