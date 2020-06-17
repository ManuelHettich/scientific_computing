Curve Approximation
===================

We want to approximate a set of given constraints `(x_i, y_i)` by a
polynomial curve. The constraints can be defined by clicking with the
left mouse button. Use the key `c` to clear the points. Use
`1` to `7` to try the seven test cases (a simple quadratic function and a sine curve with
increasing density). Press `a` and `i` to switch to approximation and
interpolation, respectively. You can use the `spacebar` to switch between
your own LU solver and one provided by [Eigen](http://eigen.tuxfamily.org/),
or between your own Cholesky solver, one provided by Eigen, and a QR solver
provided by Eigen -- depending on whether approximation or interpolation
is active.
You can find some typedefs for Eigen's matrices and vectors in `types.h`.
Use `Scalar` for floating point variables and change it's definition in
`types.h` if you want to use `float` instead of `double`.

Building under Linux/MacOS
--------------------------

Inside the exercise's top-level directory, execute the following commands:

    mkdir build
    cd build
    cmake ..
    make

Using Xcode under MacOS
-----------------------

If you prefer to use Xcode on Mac, to the following:

    mkdir build
    cd build
    cmake -G Xcode ..

and then open the generated Xcode project.

Building and running under Windows
----------------------------------

* Download and install [Visual Studio Community](https://www.visualstudio.com/vs/community/)
* Make sure that you check "Desktop development with C++" during installation
* Download [cmake](https://cmake.org/download/) via the platform windows .zip version and extract it somewhere
* Create an empty build folder inside the project's top-level directory
* Start cmake-gui.exe (located in cmake's bin folder)
* Specify the top-level directory as source directory (button Browse source...)
* Specify the build folder as build directory (button Browse build...)
* Select Configure using your Visual Studio Version as option.
* When configuration is finished, select Generate.
* Start Visual Studio Community
* Open the project via File -> open -> project -> .sln in build folder
* In the project explorer window on the right, right-click the project (approximation) and set it as startup-project
* Switch to release mode
* Hit CTRL + F5 to build and run (or CTRL + SHIFT + B to build)


Documentation
-------------

Open the pre-built documentation from within the `build` directory via

    firefox ../doc/index.html

or using any other web browser.

You can also build the HTML documentation on your own as long as you have [Doxygen](www.doxygen.org/) installed. To do so, still inside the directory `build`, execute the following command:

    make doc

View the documentation by opening the file `html/index.html` with any web browser / HTML viewer. If you are into LaTeX, navigate into the directory `latex` and execute the command `make` to create a printable version of the documentation.


Running
-------

After building (i.e. execution of `make`), still in the directory `build`, launch the program by executing the following command:

    ./approximation

Todo
----

* In the file `src/InterpolationViewer.cpp`, implement the function `fit_curve()`. Set up the linear system `A * x = b` and use either your or Eigen's solver to solve the system. Note, that your solver is not implemented yet.

* In the file `src/ApproximationViewer.cpp`, implement the function `approximate()`. Set up the overdetermined system and the normal equations, and use either your solver, Eigen's Cholesky solver or Eigen's QR solver to compute the least squares solution. Note, that your solver is not implemented yet.

* In the file `src/lu.cpp`, implement the `solve()` function that solves `A * x = b` after the factorization `A = L * U` has been computed. Here you have to solve two triangular systems involving `L` and `U`. Finish the function `factorize()` to actually compute the factorization.

* In the file `src/ApproximationViewer.cpp`, implement the function `add_noise()`. Use it to add some noise to the constraints' `y`-coordinates.

* In the file `src/cholesky.cpp`, implement the `factorize()`. Then, make this solver available as an alternative in `src/ApproximationViewer.cpp`.

License
-------

The code examples are copyright Computer Graphics Group, Bielefeld University.


Have fun!

