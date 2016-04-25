# ShapeComplexityIndex

### How to build

1. Download the source into _SRC_ directory.
2. Create a directory to build the package
3. Run CMAKE and generate Makefiles
  * ccmake SRC/ShapeComplexityIndex
  * make
4. Set the library pathway 
  * export DYLD_LIBRARY_PATH=/../../ShapeComplexityIndex_Geo/lib/
  * setenv DYLD_LIBRARY_PATH /../../ShapeComplexityIndex_Geo/lib/ 

### How to use

 -i input surface[obj|vtk]
 -o output file
 --format ASCII|KWM
    choose the output file format
      ASCII as a textfile
      KWM as a KWMesh
 -t the template MNI object file

