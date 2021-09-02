import h5py as h5
import numpy as np
import sys
def generateXDMF(filename:str, numMaterial:int):
    hdf5 = h5.File(filename,'r')
    print(hdf5.keys())
    lenHDF5Keys = (len(hdf5['vector_morphology'].keys()))
    assert(numMaterial > 0)
    if(lenHDF5Keys != numMaterial*2):  # Unaligned + Aligned for all material
        print("Please Check the number of material. XDMF not generated")
        return


    vmorph = hdf5['vector_morphology']
    unalignedShape = (vmorph['Mat_1_unaligned'])
    dataType = unalignedShape.dtype
    if(dataType == np.float64):
        precision = 8
    elif(dataType == np.float32):
        precision = 4
    else:
        assert("Wrong DataType")
    unalignedShape = unalignedShape.shape
    assert(len(unalignedShape) == 3)
    alignedShape = (vmorph['Mat_1_alignment'])
    alignedShape = alignedShape.shape
    assert(len(alignedShape) == 4)
    assert(alignedShape[0] == unalignedShape[0])
    assert(alignedShape[1] == unalignedShape[1])
    assert(alignedShape[2] == unalignedShape[2])
    assert(alignedShape[3] == 3)
    for i in range(1,numMaterial):
        unalignedShapeMat = (vmorph['Mat_'+str(i+1)+'_unaligned'])
        unalignedShapeMat = unalignedShapeMat.shape
        alignedShapeMat = (vmorph['Mat_'+str(i+1)+'_alignment'])
        alignedShapeMat = alignedShapeMat.shape
        assert(unalignedShapeMat  == unalignedShape)
        assert(alignedShapeMat  == alignedShape)



    print("Morphology OK.")
    f = open("morphology.xdmf", "w")
    ## Header section
    f.write('<?xml version="1.0" ?>\n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
    f.write('<Xdmf xmlns:xi="http://www.w3.org/2003/XInclude" Version="2.2">\n')


#     ### Grid section
    f.write('\t<Domain>\n')
    f.write('\t<Grid Name="morphology" GridType="Uniform">\n')

    f.write('\t\t<Topology TopologyType="3DCORECTMesh" NumberOfElements=" {} {} {}" />\n'.format(unalignedShape[0],unalignedShape[1],unalignedShape[2]))
    f.write('\t\t<Geometry GeometryType="ORIGIN_DXDYDZ">')


    f.write('\t\t<DataItem Name="origin" Dimensions="3" NumberType="Float" Precision="{}" Format="XML">\n'.format(precision))
    f.write('\t\t0.0 0.0 0.0\n')
    f.write('\t\t</DataItem>\n')
    f.write('\t\t<DataItem Name="spacing" Dimensions="3" NumberType="Float" Precision="{}" Format="XML">\n'.format(precision))
    f.write('\t\t1.0 1.0 1.0\n')
    f.write('\t\t</DataItem>\n')
    f.write('\t\t</Geometry>\n')

    #Attribute section
    for i in range(1,numMaterial+1):
        f.write('\t\t<Attribute Name="Mat_{}_alignment" AttributeType="Vector" Center="Cell">\n'.format(i))
        f.write('\t\t<DataItem Dimensions="{} {} {} {}" NumberType="Float" Precision="{}" Format="HDF">\n'.format(alignedShape[0],alignedShape[1],alignedShape[2],3,precision))
        f.write('\t\t{}:/vector_morphology/Mat_{}_alignment\n'.format(filename,str(i)))
        f.write('\t\t</DataItem>\n')
        f.write('\t\t</Attribute>\n')


    for i in range(1,numMaterial+1):
        f.write('\t\t<Attribute Name="Mat_{}_unaligned" AttributeType="Scalar" Center="Cell">\n'.format(i))
        f.write('\t\t<DataItem Dimensions="{} {} {}" NumberType="Float" Precision="{}" Format="HDF">\n'.format(alignedShape[0],alignedShape[1],alignedShape[2],precision))
        f.write('\t\t{}:/vector_morphology/Mat_{}_unaligned\n'.format(filename,str(i)))
        f.write('\t\t</DataItem>\n')
        f.write('\t\t</Attribute>\n')

    #footer
    f.write('\t</Grid>\n')
    f.write('\t</Domain>\n')
    f.write('</Xdmf>\n')
    hdf5.close()
    print('XDMF Generated as morphology.xdmf')

if __name__ == '__main__':

    argc:int = len(sys.argv)
    if(argc != 3):
        print("Usage : python3 {} HDF5FileName numMaterial".format(sys.argv[0]))
        exit(0)

    HDF5File: str = sys.argv[1]
    numMaterial:int = int(sys.argv[2])
    generateXDMF(HDF5File,numMaterial)
