import sys
import h5py as h5
import numpy as np

argc: int = len(sys.argv)
if(argc < 4):
    print("Usage : ",sys.argv[0], "oldFileName newFileName numMaterial")

oldFileName = sys.argv[1]
newFileName = sys.argv[2]
numMaterial = int(sys.argv[3])

fOld = h5.File(oldFileName,'r')
fnewZYX = h5.File(newFileName+'ZYX.h5','w')
fnewXYZ = h5.File(newFileName+'XYZ.h5','w')

fnewZYX.create_group('vector_morphology')
fnewXYZ.create_group('vector_morphology')

#
for i in range(1,numMaterial+1):

    datasetName = 'Mat_'+str(i)+'_unaligned'
    valOld = fOld['vector_morphology'][datasetName][()]
    fnewZYX['vector_morphology'].create_dataset(datasetName,data = valOld)
    fnewZYX['vector_morphology'][datasetName].dims[0].label = 'Z'
    fnewZYX['vector_morphology'][datasetName].dims[1].label = 'Y'
    fnewZYX['vector_morphology'][datasetName].dims[2].label = 'X'

    valXYZ = np.swapaxes(valOld,0,2)
    fnewXYZ['vector_morphology'].create_dataset(datasetName,data = valXYZ)

    fnewXYZ['vector_morphology'][datasetName].dims[0].label = 'X'
    fnewXYZ['vector_morphology'][datasetName].dims[1].label = 'Y'
    fnewXYZ['vector_morphology'][datasetName].dims[2].label = 'Z'

for i in range(1,numMaterial+1):
    datasetName = 'Mat_'+str(i)+'_alignment'
    valOld = fOld['vector_morphology'][datasetName][()]
    valOld[:,:,:,[0,2]] = valOld[:,:,:,[2,0]]

    fnewZYX['vector_morphology'].create_dataset(datasetName,data = valOld)
    fnewZYX['vector_morphology'][datasetName].dims[0].label = 'Z'
    fnewZYX['vector_morphology'][datasetName].dims[1].label = 'Y'
    fnewZYX['vector_morphology'][datasetName].dims[2].label = 'X'

    valXYZ = np.swapaxes(valOld,0,2)
    fnewXYZ['vector_morphology'].create_dataset(datasetName,data = valXYZ)
    fnewXYZ['vector_morphology'][datasetName].dims[0].label = 'X'
    fnewXYZ['vector_morphology'][datasetName].dims[1].label = 'Y'
    fnewXYZ['vector_morphology'][datasetName].dims[2].label = 'Z'

fnewXYZ.close()
fnewZYX.close()
fOld.close()

print("New format File Generated :",newFileName+"ZYX/XYZ.h5")
