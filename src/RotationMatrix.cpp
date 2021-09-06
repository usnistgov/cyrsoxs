//
// Created by maksbh on 5/4/21.
//
#include <RotationMatrix.h>
RotationMatrix::RotationMatrix(const InputData *inputData)
:inputData_(inputData){

}

void RotationMatrix::initComputation() {
  const UINT sz = inputData_->kVectors.size();
  baseConfig_.resize(inputData_->kVectors.size());
  baseAxis_.resize(inputData_->kVectors.size());
  const auto & kVecs = inputData_->kVectors;

  for(int i = 0; i < sz; i++){
    computeRotationMatrixK(kVecs[i],baseConfig_[i].matrix);
    Matrix mat;
    computeRotationMatrixBaseConfiguration(kVecs[i],baseConfig_[i].matrix,mat,baseConfig_[i].baseRotAngle);
    static constexpr Real3 X{1,0,0};
    static constexpr Real3 Y{0,1,0};
    static constexpr Real3 Z{0,0,1};
    doMatVec<false>(mat,X,baseAxis_[i].X);
    doMatVec<false>(mat,Y,baseAxis_[i].Y);
    doMatVec<false>(mat,Z,baseAxis_[i].Z);
  }


  computeRotationMatrixK(inputData_->detectorCoordinates,detectorRotationMatrix_);

}

void RotationMatrix::printToFile(std::ofstream & fout) const{
  fout << "Detector Rotation Matrix" << "\n";
  fout << "-------------------------------\n";
  detectorRotationMatrix_.printToFile(fout);
  fout << "-------------------------------\n";

  fout << "K Rotation Matrix And Base Configuration \n";
  fout << "-----------------------------------------\n";
  for(UINT i = 0; i < baseConfig_.size(); i++){
    const auto & baseCfg = baseConfig_[i];
    const auto & kVec = inputData_->kVectors[i];
    fout << "-----------------------------------------\n";
    fout << "K = " << kVec.x << " " << kVec.y << " " << kVec.z << "\n";
    baseCfg.matrix.printToFile(fout);
    fout << "RotAngle = " << baseCfg.baseRotAngle << "\n";
    fout << "Base X   = " << baseAxis_[i].X.x << " " << baseAxis_[i].X.y << " " << baseAxis_[i].X.z << "\n";
    fout << "Base Y   = " << baseAxis_[i].Y.x << " " << baseAxis_[i].Y.y << " " << baseAxis_[i].Y.z << "\n";
    fout << "Base Z   = " << baseAxis_[i].Z.x << " " << baseAxis_[i].Z.y << " " << baseAxis_[i].Z.z << "\n";
    fout << "-----------------------------------------\n";
  }
  fout << "-------------------------------\n";
}