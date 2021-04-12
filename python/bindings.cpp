#include <pybind11/eigen.h>  //For automatic eigen-numpy conversion
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "DCA/dca.h"

namespace py = pybind11;

namespace DCA {

PYBIND11_MODULE(PythonDCA, m) {
    m.doc() = "Differentiable Collision Avoidance - Python Module";

    py::class_<Sphere>(m, "Sphere")
        .def(py::init<>(),
             "Creates a unit sphere with center (0,0,0) and radius 1.")
        .def(py::init<const Vector3d&, const double&>(),
             "Creates a sphere with given position and radius.",
             py::arg("position") = Vector3d(0, 0, 0), py::arg("radius") = 1.)
        .def("compute_D",
             py::overload_cast<const Sphere&>(&Sphere::compute_D, py::const_),
             "Compute Sphere-Sphere-Distance")
        .def("compute_D",
             py::overload_cast<const Capsule&>(&Sphere::compute_D, py::const_),
             "Compute Sphere-Capsule-Distance")
        .def("compute_dDdP", py::overload_cast<VectorXd&, const Capsule&>(
                                 &Sphere::compute_dDdP, py::const_))
        .def("compute_dDdP", py::overload_cast<VectorXd&, const Sphere&>(
                                 &Sphere::compute_dDdP, py::const_))
        .def("compute_d2DdP2", py::overload_cast<MatrixXd&, const Capsule&>(
                                   &Sphere::compute_d2DdP2, py::const_))
        .def("compute_d2DdP2", py::overload_cast<MatrixXd&, const Sphere&>(
                                   &Sphere::compute_d2DdP2, py::const_))
        .def_property_readonly("position", &Sphere::getPosition)
        .def_property_readonly("radius", &Sphere::getRadius);

    py::class_<Capsule>(m, "Capsule")
        .def(py::init<>())
        .def(py::init<const Vector3d&, const Vector3d&, const double&>())
        .def("compute_D",
             py::overload_cast<const Sphere&>(&Capsule::compute_D, py::const_),
             "Compute Capsule-Sphere-Distance")
        .def("compute_D",
             py::overload_cast<const Capsule&>(&Capsule::compute_D, py::const_),
             "Compute Capsule-Capsule-Distance")
        .def("compute_dDdP", py::overload_cast<VectorXd&, const Capsule&>(
                                 &Capsule::compute_dDdP, py::const_))
        .def("compute_dDdP", py::overload_cast<VectorXd&, const Sphere&>(
                                 &Capsule::compute_dDdP, py::const_))
        .def("compute_d2DdP2", py::overload_cast<MatrixXd&, const Capsule&>(
                                   &Capsule::compute_d2DdP2, py::const_))
        .def("compute_d2DdP2", py::overload_cast<MatrixXd&, const Sphere&>(
                                   &Capsule::compute_d2DdP2, py::const_))

        .def_property_readonly("startPosition", &Capsule::getStartPosition)
        .def_property_readonly("endPosition", &Capsule::getEndPosition)
        .def_property_readonly("radius", &Capsule::getRadius);

    m.def("compute_D", &compute_D,
          "Compute the distance between two primitives, given the indices");

    m.def("compute_dDdP", &compute_dDdP,
          "Compute the gradient between two primitives, given the indices");

    m.def("compute_d2DdP2", &compute_d2DdP2, "Compute the hessian.");

    // Pair bindings
    py::class_<PermutationPairGenerator>(m, "PermutationPairGenerator")
        .def(py::init<>())
        .def("generate", &PermutationPairGenerator::generate);

#if BUILD_COMPACT_N_SEARCH == 1
    py::class_<NeighborsPairGenerator>(m, "NeighborsPairGenerator")
        .def(py::init<>())
        .def(py::init<const double&>())
        .def("generate", &NeighborsPairGenerator::generate);
#endif
}

}  // namespace DCA