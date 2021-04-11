//
// Created by maksbh on 4/11/21.
//

#ifndef CY_RSOXS_TESTUTILS_H
#define CY_RSOXS_TESTUTILS_H

#include <Datatypes.h>
#include <limits>
static inline Real computeLinfError(const Real * vec1, const Real * vec2, const BigUINT size){
  Real maxDifference = -std::numeric_limits<Real>::infinity();
  for(BigUINT i = 0; i < size; i++){
    Real difference = fabs(vec1[i] - vec2[i]);
    maxDifference = difference > maxDifference: difference?maxDifference;
  }
  return maxDifference;
}
#endif //CY_RSOXS_TESTUTILS_H
