### Organizing the Source Code
The `src` folder contains the files `NavierStokes.hpp`, `NavierStokes.cpp`, and some example files demonstrating their usage.

Binary files must not be uploaded to the repository (including executables).

The `mesh` folder contains `.geo` files for `gmsh`. To generate the meshes, navigate to the `mesh` folder and run:
```bash
$ cd mesh
$ gmsh domain2D.geo -2
$ gmsh domain3D.geo -3
```

### Compilation
To compile the executable, make sure you have loaded the necessary modules with:
```bash
$ module load gcc-glibc dealii
```
Then run the following commands:
```bash
$ mkdir build
$ cd build
$ cmake ..
$ make -j2 NavierStokes
$ make -j2 postprocess
```
The executable will be created in the `build` folder and can be executed with:
```bash
$ ./executable-name
```

### Running Tests
To run the tests, navigate to the appropriate test folder:
#### 2D Tests
```bash
$ cd tests/2D/test_0X
$ mkdir build
$ cd build
$ cmake ..
$ make -j2
```
#### 3D Tests
```bash
$ cd tests/3D/test_0X
$ mkdir build
$ cd build
$ cmake ..
$ make -j2
```

