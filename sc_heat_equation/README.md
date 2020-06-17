Diffusion
=========

We want to implement a dynamic simulation of heat diffusion on a regular 2-dimensional grid, visualized in 3D, with heat represented by height and color shifted from blue to red. For the actual simulation over time, we use explicit time integration (Euler). Then we also compute the steady-state equilibrium by solving the Laplace equation.

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
* In the project explorer window on the right, right-click the project (diffusion) and set it as startup-project
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

    ./diffusion

Controls
--------

* `Esc`: quit application
* `←`/`→`: rotate the grid.
* `s`: do a single time step
* `a`: turn animtion on/off
* `space`: compute the equilibrium by solving the Laplace equation
* `1`-`3`: setup some test cases
* left mouse button: create heat at the mouse cursor
* GUI: decrease/increase both grid resolution and time step, and toggle wireframe rendering on/off

Todo
----

* In the file `src/HeatEquationViewer.cpp`, implement the function `explicit_euler_step()`. This will allow you to execute a single time step and to run the animation.

* In the file `src/HeatEquationViewer.cpp`, implement the function `solve_equilibrium()`. Afterwards, you will be able to directly find a configuration's equilibrium. This will look most interesting for the test case accessible by pressing `3`.

License
-------

The code examples are copyright Computer Graphics Group, Bielefeld University.


Have fun!
