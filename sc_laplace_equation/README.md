Diffusion -- Gradient Descent and Conjugate Gradients
=====================================================

In this exercise we solve the Laplace equation, i.e., the steady-state
equilibrium of the Heat Equation. This solution was already computed in the
last exercise, but the dense Cholesky solver was only applicable to small
grid resolutions. This time, we use either the gradient descent solver or
the conjugate gradients solver to solve the Laplace equation on 
a 100x100 grid (or higher).

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
* In the project explorer window on the right, right-click the project (gradients) and set it as startup-project
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

    ./gradients

Controls
--------

* `Esc`: quit application
* `←`/`→`: rotate the grid.
* `s`: do a single time step
* `a`: turn animtion on/off
* `space`: compute the equilibrium by solving the Laplace equation
* `1`-`3`: setup some test cases
* Left mouse button: create heat at the mouse cursor
* Use the GUI to decrease/increase both grid resolution and time step, toggle wireframe rendering on/off, and choose the solver to compute the equilibrium

Todo
----

* Copy your results from the previous exercise to complete `src/HeatEquationViewer.cpp`.

* In the files `src/gradient_descent.h` and `src/conjugate_gradients.h`, implement the respective solvers. They are called from `src/HeatEquationViewer.cpp`, where the sparse matrix is set up. Note that we do not use Eigen's implementation of sparse matrices, but a custom variety.

* For both solvers, implement some way of storing the error after each iteration and visualize it afterwards. To plot the errors curves, you can use any software you like, e.g. Gnuplot, Octave, Matlab, or even Office.

License
-------

The code examples are copyright Computer Graphics Group, Bielefeld University.


Have fun!
