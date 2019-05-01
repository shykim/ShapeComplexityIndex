### How to build

1. Download the source into _SRC_ directory.
2. Create a directory to build the package
3. Run CMAKE and generate Makefiles
  * ccmake SRC/ShapeComplexityIndex
  * make

### How to use

```
 -o output file
 --format ASCII|KWM
    choose the output file format
      ASCII as a textfile
      KWM as a KWMesh
 -t the template MNI object file
```
