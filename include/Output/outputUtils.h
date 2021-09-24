/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 - 2021 Iowa State University
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//////////////////////////////////////////////////////////////////////////////////

#ifndef CY_RSOXS_OUTPUTUTILS_H
#define CY_RSOXS_OUTPUTUTILS_H
#include <iomanip>
#include <iostream>
#include <errno.h> //Allows use of error numbers
#include <fcntl.h> //Specified in man 2 open
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <errno.h> //Allows use of error numbers
#include <fcntl.h> //Specified in man 2 open
#include <stdio.h>
#include <sys/stat.h>
#include <cstring>

/**
 * @brief Creates the directory.
 * @param dirname The name of the directory
 */
static void createDirectory(const std::string & dirName){

  int ierr = mkdir(dirName.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (ierr != 0 && errno != EEXIST) {
    std::cout  << "Could not create folder for storing results (" <<  strerror(errno) << "\n";
    exit(EXIT_FAILURE);
  }
}
#endif //CY_RSOXS_OUTPUTUTILS_H
