![Build](https://github.com/SimiPro/DCA/actions/workflows/build.yml/badge.svg?branch=master)

# DCA - Differentiable Collision Avoidance

We implement differentiable collision avoidance using multiple shape primitives.

## Python Bindings
If the cmake option `BUILD_PYTOHN_BINDINGS` is on, the `make` command will build python bindings for all functionality of this library.

After `make` has been run successfully, run the following command:

```pip3 install -e python/PythonDCA/```

which installs the python library PythonDCA into your local python installation. You have now the module `PythonDCA` available inside python.