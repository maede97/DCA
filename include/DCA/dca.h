/**
 * This file serves as a helper such that a user
 * has only to include a single file to use this library.
 * 
 * @author Matthias Busenhart, Simon Zimmermann, Simon Huber, Stelian Coros
 * CRL Group, ETH Zurich, crl.ethz.ch
 * (c) 2021
 */

#ifndef __DCA_H__
#define __DCA_H__

// All base DCA files
#include "newton.h"
#include "pair.h"
#include "primitives.h"
#include "utils.h"

// And finally implementations of the primitives:
#include "Capsule.h"
#include "Sphere.h"

#endif /* __DCA_H__ */