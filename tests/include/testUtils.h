//
// Created by maksbh on 4/11/21.
//

#ifndef CY_RSOXS_TESTUTILS_H
#define CY_RSOXS_TESTUTILS_H

#include <Datatypes.h>
#include <limits>

static constexpr Real TOLERANCE_CHECK = 1E-3;

static inline Real computeLinfError(const Real * vec1, const Real * vec2, const BigUINT size){
  Real maxDifference = -std::numeric_limits<Real>::infinity();
  for(BigUINT i = 0; i < size; i++){
    Real difference = fabs(vec1[i] - vec2[i]);
    maxDifference = difference > maxDifference ? difference:maxDifference;
  }
  return maxDifference;
}

static inline Complex computeLinfError(const Complex * vec1, const Complex * vec2, const BigUINT size){
  Complex maxDifference;
  maxDifference.x = -std::numeric_limits<Real>::infinity();
  maxDifference.y = -std::numeric_limits<Real>::infinity();
  for(BigUINT i = 0; i < size; i++){
    Complex difference{fabs(vec1[i].x - vec2[i].x),fabs(vec1[i].y - vec2[i].y)};
    maxDifference.x = difference.x > maxDifference.x? difference.x:maxDifference.x;
    maxDifference.y = difference.y > maxDifference.y? difference.x:maxDifference.y;
  }
  return maxDifference;
}

template<typename T>
static inline void readFile(T * vec1, const std::string & filename, const BigUINT size){
  FILE *fstream = fopen(filename.c_str(),"rb");
  fread(vec1, sizeof(T),size,fstream);
  fclose(fstream);

}


#endif //CY_RSOXS_TESTUTILS_H
