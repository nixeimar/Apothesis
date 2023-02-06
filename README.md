<p align="center"> <img width="200" src="https://github.com/nixeimar/Apothesis/blob/master/ApothesisLogo_small.png" alt="Material Bread logo" > </p>
  

Apothesis 
--------------------------------------------------------------------------------------------------------------
Apothesis is an open source software for simulating deposition processes via the kinetic Monte Carlo method. 
This is the first step for creating a generalized kinetic Monte Carlo code 
for surface growth phenomena targeting chemical vapor and atomic layer deposition processes.  
It is still under development but please feel free to contact me if you have something in mind! 

Compile - qmake
--------------------------------------------------------------------------------------------------------------
On linux systems in the Apothesis directory run 
```
qmake
make
```
This should do the trick.
Having QtCreator will make things a lot easier on (mostly) windows and linux OS. 
Since the project is not based on Qt framework you can use any IDE of your choice. 
I will provided a cmake at first instance. 

Compile - cmake
--------------------------------------------------------------------------------------------------------------
On linux systems, navigate to build
``` 
cmake ..
make -j
```

Contact information:
Nikolaos (Nikos) Cheimarios: 
nixeimar@chemeng.ntua.gr
