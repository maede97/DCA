#include <DCA/Pair.h>
#include <DCA/Sphere.h>
#include <DCA/Capsule.h>

namespace DCA {

std::vector<pair_t> PermutationPairGenerator::generate(
    const std::vector<primitive_t> &primitives) const {
    // Create all possible permutations.
    std::vector<pair_t> ret;
    ret.reserve(primitives.size() * primitives.size() - primitives.size());

    for (size_t i = 0; i < primitives.size(); i++) {
        for (size_t j = 0; j < primitives.size(); j++) {
            // skip primitve self-collision
            if (i == j) continue;
            ret.push_back({i, j});
        }
    }

    return ret;
}

#if BUILD_COMPACT_N_SEARCH
std::vector<pair_t> NeighborsPairGenerator::generate(
    const std::vector<primitive_t> &primitives) const {
    std::vector<pair_t> ret;

    std::vector<std::array<double, 3>> positions;
    positions.resize(primitives.size());
    // Create the position vector

    for (size_t i = 0; i < primitives.size(); i++) {
        Vector3d pos = std::visit(
            overloaded{[](const Sphere &cp) { return cp.getPosition(); },
                       [](const Capsule &cp) {
                           return Vector3d(0.5 * (cp.getStartPosition() +
                                                  cp.getEndPosition()));
                       },
                       [](const primitive_t &cp) {
                           throw std::logic_error(
                               "Could not get primitive position. Missing "
                               "overload.");
                           return 0.;
                       }},
            primitives[i]);
        positions[i] = {pos.x(), pos.y(), pos.z()};
    }

    CompactNSearch::NeighborhoodSearch nsearch(m_radius);
    unsigned id =
        nsearch.add_point_set(positions.front().data(), positions.size());

    // Actually perform the computation
    nsearch.find_neighbors();

    // Now get all neighbors and add them to the pairs list
    CompactNSearch::PointSet const &ps = nsearch.point_set(id);
    for (int i = 0; i < ps.n_points(); i++) {
        for (size_t j = 0; j < ps.n_neighbors(id, i); j++) {
            // return point id of the j-th neighbor of the i-th particle in the point set
            const unsigned pid = ps.neighbor(id, i, j);

            ret.push_back({i, pid});
        }
    }

    return ret;
}

#endif

}  // namespace DCA