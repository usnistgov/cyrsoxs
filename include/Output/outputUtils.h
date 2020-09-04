//
// Created by maksbh on 9/4/20.
//

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
