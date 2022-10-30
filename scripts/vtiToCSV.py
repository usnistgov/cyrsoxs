# trace generated using paraview version 5.5.2
import sys
# import the simple module from the paraview
from paraview.simple import *


# create a new 'XML Image Data Reader'
projectionAverage_refvti = XMLImageDataReader(FileName=sys.argv[1])
projectionAverage_refvti.CellArrayStatus = ['projection']

# save data
SaveData(sys.argv[2], proxy=projectionAverage_refvti,FieldAssociation='Cells')
