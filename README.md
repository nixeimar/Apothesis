![My Image](ApothesisLogo.jpg)

General
--------------------------------------------------------------------------------------------------------------
Apothesis is an open source software for simulating deposition processes via the kinetic Monte Carlo method. 
This is the first step for creating a generalized kinetic Monte Carlo code 
for surface growth phenomena targeting chemical vapor and atomic layer deposition processes.  
It is still under development but please feel free to contact me if you have something in mind! 

Compile - qmake
--------------------------------------------------------------------------------------------------------------
On linux systems in the src directory run 
```
qmake
make
```
This should do the trick.
Having QtCreator will make things a lot easier on (mostly) windows and linux OS. 
Since the project is not based on Qt framework (although it started like that - thats why the Qt deps) 
I will provided a cmake at first instance. 

Compile - cmake
--------------------------------------------------------------------------------------------------------------
On linux systems, navigate to src/build
``` 
cmake ..
make -j
```

Contact information:
Nikolaos (Nikos) Cheimarios: 
nixeimar@chemeng.ntua.gr
