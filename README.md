<p align="center"> <img width="200" src="https://github.com/nixeimar/Apothesis/blob/master/ApothesisLogo_small.png" alt="Material Bread logo" > </p>
  

Apothesis 
--------------------------------------------------------------------------------------------------------------
Apothesis is an open source software for simulating deposition processes via the kinetic Monte Carlo method. 
This is the first step for creating a generalized kinetic Monte Carlo code 
for surface growth phenomena targeting chemical vapor and atomic layer deposition processes.  
It is still under development but please feel free to contact me if you have something in mind! 

Compile - qmake
--------------------------------------------------------------------------------------------------------------
In the Apothesis directory run 
```
qmake
```

and then 

```
make
```
From `src/input.kmc` copy the `input.kmc` file.

Paste the `input.kmc` file in the `Apothesis` folder.

And then in terminal type
```
./apothesis

```
to run the simulation.

This should do the trick.

Having QtCreator will make things a lot easier on (mostly) windows and linux OS. 
If you are more comfortable using Qt GUI then follow,

Open Qt creater and go to `Open Projects` and open the `Apothesis.pro` file.
Then hit `Configure Project`.

After this move up a directory from the `Apothesis` directory and you should see the Qt build folder `build-Apothesis-Desktop-Debug`.
Paste the `input.kmc` file there.

Now you can run the application.

Since the project is not based on Qt framework you can use any IDE of your choice. 

Compile - cmake
--------------------------------------------------------------------------------------------------------------
In the Apothesis directory create build directory.
Inside the build directory type 

``` 
cmake ..
``` 

and then 

``` 
make -j
```
From `src/input.kmc` copy the `input.kmc` file and put it inside the build folder and then in terminal type
```
./Apothesis

```

Contact information:
Nikolaos (Nikos) Cheimarios: 
nixeimar@chemeng.ntua.gr
