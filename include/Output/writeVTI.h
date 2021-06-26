/////////////////////////////////////////////////////////////////////////////////
// MIT License
//
//Copyright (c) 2019 -2020 Iowa State University
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

#ifndef CUDA_BASE_WRITEVTI_H
#define CUDA_BASE_WRITEVTI_H

#include <Datatypes.h>
#include <fstream>
#include <Input/Input.h>
#include <iostream>
#include <vector>
#include <vector_types.h>
#include <Output/cencode.h>
#include <assert.h>
#include <omp.h>
#include <chrono>
#include <cstring>


/// These function are writing in vti file format. Refer to
/// <a href="https://www.paraview.org/Wiki/ParaView/ParaView_Readers_and_Parallel_Data_Distribution">Paraview</a>
/// for more details. The file writing is provided through both C++ and C style.
/// In general, C style of file writing gave almost 1.5X speedup.

namespace VTI{
/**
 * @brief This function writes data to vti using base64 encoding and fstream object.
 * @param fout fstream object
 * @param numeric_data the pointer to numeric data
 * @param byte_length
 * @return 0 for success and -1 for failure
 */
static int vtk_write_binary (std::ofstream & fout, char *numeric_data, size_t byte_length)
{


  size_t              chunks, chunksize, remaining, writenow;
  size_t              code_length, base_length;
  uint32_t            int_header;
  char                *base_data;
  base64_encodestate  encode_state;

  /* VTK format used 32bit header info */
  assert (byte_length <= (size_t) UINT32_MAX);

  /* This value may be changed although this is not tested with VTK */
  chunksize = (size_t) 1u << 15; /* 32768 */
  int_header = (uint32_t) byte_length;

  /* Allocate sufficient memory for base64 encoder */
  code_length = 2 * std::max (chunksize, sizeof (int_header));
  code_length = std::max (code_length, (size_t)4) + 1;
  base_data = (char*)calloc(code_length,sizeof(char));// (code_length*sizeof(char));

  base64_init_encodestate (&encode_state);

  base_length =base64_encode_block ((char *) &int_header, sizeof (int_header), base_data,&encode_state);

  assert (base_length < code_length);
  base_data[base_length] = '\0';
  for(int i = 0; i < base_length;i++){
    fout << base_data[i] ;
  }
  chunks = 0;
  remaining = byte_length;
  while (remaining > 0) {
    writenow = std::min (remaining, chunksize);

    base_length = base64_encode_block (numeric_data + chunks * chunksize,writenow, base_data, &encode_state);

    assert (base_length < code_length);
    base_data[base_length] = '\0';

    for(int i = 0; i < base_length;i++){
      fout << base_data[i] ;
    }
    remaining -= writenow;
    ++chunks;
  }

  base_length = base64_encode_blockend (base_data, &encode_state);
  assert (base_length < code_length);
  base_data[base_length] = '\0';
  free(base_data);

  if(not(fout.good())){
    return -1;
  }
  return 0;
}

/**
 * @brief This function writes data to vti using base64 encoding and file pointer.
 * @param fp file pointer
 * @param numeric_data the pointer to numeric data
 * @param byte_length
 * @return 0 for success and -1 for failure
 */
static int vtk_write_binary (FILE * fp, char *numeric_data, size_t byte_length)
{


  size_t              chunks, chunksize, remaining, writenow;
  size_t              code_length, base_length;
  uint32_t            int_header;
  char                *base_data;
  base64_encodestate  encode_state;

  /* VTK format used 32bit header info */
  assert (byte_length <= (size_t) UINT32_MAX);

  /* This value may be changed although this is not tested with VTK */
  chunksize = (size_t) 1u << 15; /* 32768 */
  int_header = (uint32_t) byte_length;

  /* Allocate sufficient memory for base64 encoder */
  code_length = 2 * std::max (chunksize, sizeof (int_header));
  code_length = std::max (code_length, (size_t)4) + 1;
  base_data = (char*)calloc(code_length,sizeof(char));// (code_length*sizeof(char));

  base64_init_encodestate (&encode_state);

  base_length =base64_encode_block ((char *) &int_header, sizeof (int_header), base_data,&encode_state);


  assert (base_length < code_length);
  base_data[base_length] = '\0';

  (void) fwrite (base_data, 1, base_length, fp);

  chunks = 0;
  remaining = byte_length;
  while (remaining > 0) {
    writenow = std::min (remaining, chunksize);

    base_length = base64_encode_block (numeric_data + chunks * chunksize,writenow, base_data, &encode_state);

    assert (base_length < code_length);
    base_data[base_length] = '\0';

    (void) fwrite (base_data, 1, base_length, fp);
    remaining -= writenow;
    ++chunks;
  }

  base_length = base64_encode_blockend (base_data, &encode_state);
  assert (base_length < code_length);
  base_data[base_length] = '\0';
  (void) fwrite (base_data, 1, base_length, fp);
  free(base_data);
  if (ferror (fp)) {
    return -1;
  }
  return 0;
}

/**
 * @brief This function writes the file heade for 3D data
 * @param VoxelSize
 * @param out fstream object
 */
static void writeFileHeader(const UINT *VoxelSize, std::ofstream &out) {

  // The header
  out << "<?xml version=\"1.0\"?>" << std::endl;
  out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  out << "<ImageData WholeExtent=\"";
  out << "0 " << VoxelSize[0] << " ";
  out << "0 " << VoxelSize[1] << " ";
  out << "0 " << VoxelSize[2] << "\" ";
  out << "Origin=\"0 0 0\" ";
  out << "Spacing=\"" << 1 << " " << 1 << " " << 1 << "\">" << std::endl;
  out << "<Piece Extent=\"";
  out << "0 " << VoxelSize[0] << " ";
  out << "0 " << VoxelSize[1] << " ";
  out << "0 " << VoxelSize[2] << "\">" << std::endl;

}

/**
 * @brief This function writes the file heade for 2D data
 * @param VoxelSize
 * @param out fstream object
 */

static void writeFileHeader2D(const UINT *VoxelSize, std::ofstream &out) {

  // The header
  out << "<?xml version=\"1.0\"?>" << std::endl;
  out << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
  out << "<ImageData WholeExtent=\"";
  out << "0 " << VoxelSize[0] << " ";
  out << "0 " << VoxelSize[1] << " ";
  out << "0 " << 1 << "\" ";
  out << "Origin=\"0 0 \" ";
  out << "Spacing=\"" << 1 << " " << 1  << "\">" << std::endl;
  out << "<Piece Extent=\"";
  out << "0 " << VoxelSize[0] << " ";
  out << "0 " << VoxelSize[1] << " ";
  out << "0 " << 2 << "\">" << std::endl;

}
/**
 * @brief This function writes the file header for 2D data
 * @param VoxelSize voxel size
 * @param fp file pointer
 */

static void writeFileHeader2D(const UINT *VoxelSize, FILE * fp) {

  fprintf(fp,"<?xml version=\"1.0\"?>\n");
  fprintf(fp,"<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
  fprintf(fp,"<ImageData WholeExtent=\"0 %d 0 %d 0 1 \" Origin=\" 0 0 \" Spacing=\"1 1 \">\n",VoxelSize[0],VoxelSize[1]);
  fprintf(fp,"<Piece Extent=\" 0 %d 0 %d 0 2 \">\n",VoxelSize[0],VoxelSize[1]);

}


/**
 * @brief This function writes the header for scalar data
 * @param out fstream object
 * @param var_name variable name
 * @param isComplex if the input is a complex or real
 */
static void writeVariableHeaderScalar(std::ofstream &out, const char * var_name, const bool isComplex = false){
  if(isComplex){
  //  out << " <CellData Vectors=\"" << var_name << "\">\n";
#if VTI_BINARY
    out << " <DataArray type=\"Float32\" Name=\"" << var_name << "\" NumberOfComponents=\"2\" format=\"binary\">" << std::endl;
#else

    out << " <DataArray type=\"Float32\" Name=\"" << var_name << "\" NumberOfComponents=\"2\" format=\"ascii\">" << std::endl;
#endif
  }

  else{

#if VTI_BINARY
    out << " <DataArray type=\"Float32\" Name=\"" << var_name << "\" format=\"binary\">" << std::endl;
#else
    out << " <DataArray type=\"Float32\" Name=\"" << var_name << "\" format=\"ascii\">" << std::endl;
#endif

  }

}

/**
 * This function writes the header for scalar data
 * @param fp file pointer
 * @param var_name variable name
 * @param isComplex if the input is a complex or real
 */

static void writeVariableHeaderScalar(FILE * fp, const char * var_name, const bool isComplex = false){
  if(isComplex){
    //  out << " <CellData Vectors=\"" << var_name << "\">\n";
#if VTI_BINARY
    fprintf(fp," <DataArray type=\"Float32\" Name=\"%s \" NumberOfComponents=\"2\" format=\"binary\">\n",var_name);
#else

    fprintf(fp," <DataArray type=\"Float32\" Name=\"%s \" NumberOfComponents=\"2\" format=\"ascii\">\n",var_name);
#endif
  }

  else{
    // out << " <CellData Scalars=\"" << var_name << "\">\n";

#if VTI_BINARY
    fprintf(fp,"<DataArray type=\"Float32\" Name=\"%s\" format=\"binary\">\n",var_name);
#else
    fprintf(fp,"<DataArray type=\"Float32\" Name=\"%s\" format=\"ascii\">\n",var_name);
#endif

  }

}
/**
 * This function writes the header for the variable
 * @param out fstream object
 * @param var_name  variable name
 * @param isInitial if the function is called for the first time.
 */
static void writeVariableHeaderVector(std::ofstream &out, const char * var_name, bool isInitial = true){
  if(isInitial) {

  }
#if VTI_BINARY
  out << " <DataArray type=\"Float32\" Name=\"" << var_name << "\" NumberOfComponents=\"3\" format=\"binary\">" << std::endl;
#else
  out << " <DataArray type=\"Float32\" Name=\"" << var_name << "\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
#endif
}

/**
 * This function writes the footer to the file
 * @param out fstream object
 */
static void writeFileFooter(std::ofstream &out) {
  out << "</Piece>\n"
         "</ImageData>\n"
         "</VTKFile>";
}

/**
 * This function writes the footer to the file
 * @param fp file pointer object
 */
static void writeFileFooter(FILE * fp) {
  fprintf(fp, "</Piece>\n"
         "</ImageData>\n"
         "</VTKFile>\n");
}
/**
 * This function writes the variable footer to the file
 * @param out fstream object
 */
static void writeVariableFooter(std::ofstream &out) {
  out << "</DataArray>\n";
}

/**
 * This function writes the variable footer to the file
 * @param fp file pointer
 */
static void writeVariableFooter(FILE * fp) {
  fprintf(fp,"</DataArray>\n");
}

/**
 * @brief This function writes the voxel scalar data
 * @param data input voxel  data
 * @param voxelSize 3D voxel size
 * @param fname name of files
 * @param varname variable name
 */
static void writeVoxelDataScalar(const Voxel *data,
                          const UINT *voxelSize,
                          const std::string &fname,
                          const char **varname){
  std::string filename = fname+".vti";
  std::ofstream fout(filename.c_str());
  BigUINT totalSize = voxelSize[0]*voxelSize[1]*voxelSize[2];

  writeFileHeader(voxelSize, fout);
#if VTI_BINARY
  Real  * scalardata = new Real[totalSize*1];
#endif
  fout << " <CellData Scalars=\"" << "CellData"  << "\">\n";
  for(int numMat = 0; numMat < NUM_MATERIAL; numMat++) {
    writeVariableHeaderScalar(fout,varname[numMat]);
#if VTI_BINARY
    for(int i = 0; i < totalSize; i++){
     scalardata[i] =   data[i].s1[numMat].w;
    }
    vtk_write_binary(fout, (char *) scalardata, sizeof(Real) * totalSize);
#else
    for(int i = 0; i < totalSize; i++){
      fout  <<   data[numMat*totalSize+i].s1.w << " ";
    }
#endif
    writeVariableFooter(fout);
  }
#if VTI_BINARY
  delete[] scalardata;
#endif
  fout << "</CellData>\n";
  writeFileFooter(fout);

  fout.close();
}

/**
 * @brief This function writes the voxel vector data
 * @param data input voxel  data
 * @param voxelSize 3D voxel size
 * @param fname name of files
 * @param varname variable name
 */
static void writeVoxelDataVector(const Voxel *data,
                          const UINT *voxelSize,
                          const std::string &fname,
                          const char **varname){

  std::string filename = fname+".vti";
  std::ofstream fout(filename.c_str());
  BigUINT totalSize = voxelSize[0]*voxelSize[1]*voxelSize[2];
  writeFileHeader(voxelSize,fout);
#if VTI_BINARY
  Real  * vecdata = new Real[totalSize*3];
#endif
  fout << " <CellData Vectors=\"" << "CellData"  << "\">\n";
  for(int numMat = 0; numMat < NUM_MATERIAL; numMat++) {
    writeVariableHeaderVector(fout,varname[numMat]);
#if VTI_BINARY

    for (int j = 0; j < totalSize; j++) {
      vecdata[3 * j + 0] = data[j].s1[numMat].x;
      vecdata[3 * j + 1] = data[j].s1[numMat].y;
      vecdata[3 * j + 2] = data[j].s1[numMat].z;
    }
    vtk_write_binary(fout, (char *) data, sizeof(Real) * totalSize * 3);
#else
    for(BigUINT j = 0; j < totalSize; j++){

      fout << data[totalSize*numMat + j].s1.x << " " << data[totalSize*numMat + j].s1.y << " " << data[totalSize*numMat + j].s1.z << "\n";
    }
#endif
    writeVariableFooter(fout);
  }
#if VTI_BINARY
  delete [] vecdata;
#endif
  fout << "</CellData>\n";
  writeFileFooter(fout);
  fout.close();
}

/**
 * @brief This function writes the voxel scalar data
 * @param data input voxel  data
 * @param voxelSize 3D voxel size
 * @param fname name of files
 * @param varName variable name
 */
static void writeDataScalar(Real * data,const UINT *voxelSize,
                     const std::string &fname, const char * varName){
  std::string filename = fname+".vti";
  std::ofstream fout(filename.c_str());
  BigUINT totalSize = voxelSize[0]*voxelSize[1]*voxelSize[2];


  writeFileHeader(voxelSize, fout);
  fout << " <CellData Scalars=\"" << "CellData" << "\">\n";
  writeVariableHeaderScalar(fout,varName,false);
#if VTI_BINARY
  vtk_write_binary(fout,(char *)data, sizeof(Real)*totalSize);
#else
  for(int i = 0; i < totalSize; i++){
    fout << data[i]<< " "  << " ";
  }
#endif
  writeVariableFooter(fout);
  fout << "</CellData>\n";
  writeFileFooter(fout);
  fout.close();
}

/**
 * @brief This function writes the scalar 2D data
 * @param data input voxel  data
 * @param voxelSize 3D voxel size
 * @param fname name of files
 * @param varName variable name
 */
static void writeDataScalar2D(Real * data,const UINT *voxelSize,
                       const std::string &fname, const char * varName){

  std::string filename = fname+".vti";
  std::ofstream fout(filename.c_str());
  BigUINT totalSize = voxelSize[0]*voxelSize[1];


  writeFileHeader2D(voxelSize, fout);
  fout << " <CellData Scalars=\"" << "CellData" << "\">\n";
  writeVariableHeaderScalar(fout,varName,false);
#if VTI_BINARY
  int retVal = vtk_write_binary(fout,(char *)data, sizeof(Real)*totalSize);
  if(retVal == -1){
    std::cerr << "[I/O] error in file writing\n";
  }
#else
  for(int i = 0; i < totalSize; i++){
    if(std::isnan(data[i])) {
      fout << 0 << "\n";
    }
    else {
      fout << data[i]<< " "  << " ";
    }

  }
#endif
  writeVariableFooter(fout);
  fout << "</CellData>\n";
  writeFileFooter(fout);
  fout.close();
}

/**
 * @brief This function writes the scalar 3D data
 * @param data input voxel  data
 * @param voxelSize 3D voxel size
 * @param fname name of files
 * @param varName variable name
 */
static void writeDataScalar(Complex * data,const UINT *voxelSize,
                     const std::string &fname, const char * varName){
  std::string filename = fname+".vti";
  std::ofstream fout(filename.c_str());
  BigUINT totalSize = voxelSize[0]*voxelSize[1]*voxelSize[2];


  writeFileHeader(voxelSize, fout);
  fout << " <CellData Vectors=\"" << "CellData" << "\">\n";
  writeVariableHeaderScalar(fout,varName,true);

#if VTI_BINARY
  Real * compData = new Real[totalSize*2];
  for(int i = 0; i < totalSize; i++){
    compData[2*i+0] = data[i].x;
    compData[2*i+1] = data[i].y;
  }
  int retVal = vtk_write_binary(fout,(char *)compData, sizeof(Real)*totalSize*2);
  if(retVal == -1){
    std::cerr << "[I/O] error in file writing\n";
  }
  delete[] compData;
#else
  for(int i = 0; i < totalSize; i++){
    fout << data[i].x << " " << data[i].y << " ";
  }
#endif
  writeVariableFooter(fout);
  fout << "</CellData>\n";
  writeFileFooter(fout);
  fout.close();
}

/**
 * @brief This function writes the scalar 2D data in C style using file pointer
 * @param data input voxel  data
 * @param voxelSize 3D voxel size
 * @param fname name of files
 * @param varName variable name
 */
static void writeDataScalar2DFP(Real * data,const UINT *voxelSize,
                       const std::string &fname, const char * varName){

  std::string filename = fname+".vti";
  FILE * fp = fopen(filename.c_str(),"w");

  BigUINT totalSize = voxelSize[0]*voxelSize[1];


  writeFileHeader2D(voxelSize, fp);
  fprintf(fp, "<CellData Scalars=\"CellData\">\n");
  writeVariableHeaderScalar(fp,varName,false);

#if VTI_BINARY
  int retVal = vtk_write_binary(fp,(char *)data, sizeof(Real)*totalSize);

  if(retVal == -1){
    std::cerr << "[I/O] error in file writing\n";
  }
#else
  for(int i = 0; i < totalSize; i++){
    if(std::isnan(data[i])) {
      fprintf(fp, "%e ", 0.0);
    }
    else {
      fprintf(fp, "%e ", data[i]);
    }
  }
#endif
  writeVariableFooter(fp);
  fprintf(fp,"</CellData>\n");
  writeFileFooter(fp);
  fclose(fp);
}

}
#endif //CUDA_BASE_WRITEVTI_H
