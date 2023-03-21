# flippy demos

each sub-folder contains a simulation based on flippy.
In each folder, you will (at least) find:
 - `main.cpp` file containing the `c++` code of the simulation
 - `CMakeLists.txt` file that contains build instructions for the `c++` code
 - `Readme.md` file that explains what the simulation is about

## how to compile
In order to compile the simulations, you can use the CMake file. From the shell, this can be done in a few steps:
1. open the command line in the demonstration folder or navigate to it.
2. create a build directory
```shell
mkdir cmake-build
```
3. change to the build directory
```shell
cd cmake-build
```
4. create build files with CMake
```shell
cmake ..
```
5. build the executable 
```shell
cmake --build . --config Release
```
The executable file will now be located in the `cmake-build` folder and have the same name as specified in `CMakeLists.txt` in the `project(...)` statement. One caveat for `MSVC` users, the executable will be located in the `Release` sub-folder.

This will work if you cloned the whole git repository and did not change the folder structure. `CMakeLists.txt` relies on the relative position of the `demo` and `flippy` folders. 
```
.
├── assets
├── demo
│        ├── just_a_wobbly_sphere_MC
│        └── simplest_MC
├── flippy
│        ├── external
│        └── utilities
├── single_header_flippy
└── tests
         ├─ ...
         ...
```
Alternatively, you can use the single header flippy if you want to compile from the command line directly. Just copy and paste `flippy.hpp` from the `single_header_flippy` folder to the same folder as the `main.cpp` file of the demo that you are compiling and run:
```shell
g++ --std=c++2a -O3 main.cpp -o demo
```
This will create an executable named demo in the same folder.
There are two caveats with the above line of code:
1. `g++` might need to be replaced with another compiler name, depending on what is installed on your system. For example, `g++-11` if you separately installed a newer `gcc` compiler or `clang++` if you use the compilers' clang family.
2. The name of the executable will need the `.exe` suffix on Windows.
3. `--std=c++2a` and `-O3` are important flags for the compilation, the first one sets the language standard to `c++20` and is strictly necessary for the code to compile. The second one sets the optimization level to `O3`, the highest possible value. This is not strictly necessary, but it will make the code run significantly faster.
## how to use shell
### MacOS
Open the `Terminal` app, which comes preinstalled. If you have never used this app before, then you will need to first install the command line tools by executing
```shell
xcode-select —install
```
a dialog window will appear, which you should confirm by clicking install.
This step will take a bit and will install command line tools and a `c++` compiler on your mac.

### Windows
The easiest way to get a shell on windows (that I am aware of) is to install the [`PowerShell`](https://www.microsoft.com/en-us/p/powershell/9mz1snwt0n5d#activetab=pivot:overviewtab).
### Linux
If you are a Linux user, you know what to do ;)