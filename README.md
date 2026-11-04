<p align="center"> <img width="200" src="https://github.com/nixeimar/Apothesis/blob/master/ApothesisLogo_small.png" alt="Material Bread logo" > </p>
  

Apothesis 
--------------------------------------------------------------------------------------------------------------
Apothesis is an open source software for simulating deposition processes via the kinetic Monte Carlo method. 
This is the first step for creating a generalized kinetic Monte Carlo code 
for surface growth phenomena targeting chemical vapor and atomic layer deposition processes.  
It is still under development but please feel free to contact me if you have something in mind! 

## Simulation input

Use the versioned JSON format for new simulations. It supports named parameters
on the system, lattice, species, processes and kinetic models, including custom
nested values. Start with [`input.json`](input.json) and the
[input format guide](doc/input-format.md).

```sh
./apothesis --validate-input input.json
./apothesis --input input.json
```

Running `./apothesis` without arguments selects `input.json` in the working
directory. Input files must use the versioned JSON format.

Compile - qmake
--------------------------------------------------------------------------------------------------------------
To compile Apothesis using qmake, follow these steps:

Step 1: Compile the Code:
In the Apothesis directory run 
```
qmake
```

and then 

```
make
```
Step 2: Run with the JSON input supplied in the project root:

```sh
./apothesis --input input.json
```

This should do the trick.

Having QtCreator will make things a lot easier on (mostly) windows and linux OS. 
If you are more comfortable using Qt GUI then follow,

Open Qt creater and go to `Open Projects` and open the `Apothesis.pro` file.
Then hit `Configure Project`.

After this move up a directory from the `Apothesis` directory and you should see the Qt build folder `build-Apothesis-Desktop-Debug`.
Copy `input.json` from the project root into that build folder and set the
application arguments to `--input input.json`.

Now you can run the application.

Since the project is not based on Qt framework you can use any IDE of your choice. 

Compile - cmake
--------------------------------------------------------------------------------------------------------------
In the Apothesis directory create a "build" directory.
Inside the "build" directory type 

``` 
cmake ..
``` 

and then 

``` 
make -j
```
Copy `input.json` from the project root into the build folder, then run:
```
./apothesis --input input.json

```

## Build on Windows (MSYS2 + MinGW)

1. Install MSYS2
2. Install dependencies:
   ```pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-cmake```

3. Open "MSYS2 MinGW64" terminal

4. Build:
   ```
   mkdir build
   cd build
   cmake -G "MinGW Makefiles" -DCMAKE_CXX_STANDARD=17 ..
   mingw32-make
   ```

Contact information:
Nikolaos (Nikos) Cheimarios: 
nixeimar@chemeng.ntua.gr
