//
// Created by maksbh on 4/10/21.
//

#ifndef CY_RSOXS_SAMPLE_H
#define CY_RSOXS_SAMPLE_H

#include <Datatypes.h>
#include <cudaMain.h>
double returnVal(){
  return 1.0;
}


TEST(CyRSoXS, testFramework) {
warmup();

EXPECT_EQ(1.0,returnVal());

}
#endif //CY_RSOXS_SAMPLE_H
