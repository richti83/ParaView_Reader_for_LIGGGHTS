input1 = self.GetInputDataObject(0, 0)
input2 = self.GetInputDataObject(0, 1)
output = self.GetOutputDataObject(0)

newPoints = vtk.vtkPoints()
newcells = vtk.vtkCellArray()

P=input1.GetPointData()
B=input2.GetPointData()
BODIES=input1.GetNumberOfPoints()
EDGES=input2.GetNumberOfPoints()
SURFACES=input2.GetNumberOfCells()
c=0
for body in range(0,BODIES):
 r=input1.GetPoint(body)
 M=P.GetAbstractArray('M')
 T11=M.GetComponent(body,0)
 T12=M.GetComponent(body,1)
 T13=M.GetComponent(body,2)
 T21=M.GetComponent(body,3)
 T22=M.GetComponent(body,4)
 T23=M.GetComponent(body,5)
 T31=M.GetComponent(body,6)
 T32=M.GetComponent(body,7)
 T33=M.GetComponent(body,8)
 print T33
 for surface in range(0,SURFACES):
  cell = input2.GetCell(surface)
  COUNT=cell.GetNumberOfPoints()
  for j in range(0,COUNT): #alle punkte der jeweiligen zelle
   id=cell.GetPointId(j)
   #print id
   coord = input2.GetPoint(id)
   x, y, z = coord[:3]
   #print x,y,z
   #print r
   xnew=T11*x+T12*y+T13*z+r[0]
   ynew=T21*x+T22*y+T23*z+r[1]
   znew=T31*x+T32*y+T33*z+r[2]
   newPoints.InsertPoint(c, xnew, ynew, znew)
   c+=1
  newcells.InsertNextCell(4)
  newcells.InsertCellPoint(c-1)
  newcells.InsertCellPoint(c-2)
  newcells.InsertCellPoint(c-3)
  newcells.InsertCellPoint(c-4)
output.SetPoints(newPoints)
output.SetCells(9,newcells)

