
1. Check the Slicer's generateCLP. If not, you have to download and compile.

2. if your surface has 40962 vertex, template.obj has to be changed by template40962.obj

3. compile SCI code
   a. create build folder and cd build folder.
   b. ccmake ../SCI_source_code folder
   c. set the pathway of generateCLP
   d. make

4. usage 
   --> this version can run with a part of vertex to save the processing time.
	example) if you want to devide 4 part and run
	Cal_Complex -i in_surf.obj -o in_surf_SCI_part01.txt -s 3 -k 4 -L -t template.obj -d 4 -p 1
	Cal_Complex -i in_surf.obj -o in_surf_SCI_part02.txt -s 3 -k 4 -L -t template.obj -d 4 -p 2
	Cal_Complex -i in_surf.obj -o in_surf_SCI_part03.txt -s 3 -k 4 -L -t template.obj -d 4 -p 3
	Cal_Complex -i in_surf.obj -o in_surf_SCI_part04.txt -s 3 -k 4 -L -t template.obj -d 4 -p 4

   --> and Murge out_sci between 1 to 4.	

5. merge SCI file (matlab code), if you want.
   usage: MurgeSCI(in_surf.obj)


