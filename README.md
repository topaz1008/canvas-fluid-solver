HTML5 Canvas Fluid Solver
==========================

A Simple fluid solver implementation in javascript.

A live demo can be found [HERE](https://topaz1008.github.io/canvas-fluid-solver).

The demo uses [dat.gui](https://github.com/dataarts/dat.gui) for the GUI (included in the `vendors` folder).

Simulates the [Navier–Stokes](https://en.wikipedia.org/wiki/Navier-Stokes_equations) equations for incompressible fluids.

Largely based on Jos Stam's paper [Real-Time Fluid Dynamics for Games](https://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf).

The density field is rendered on an off-screen bitmap and blit'ed to the canvas for performance reasons.

The code is commented as best I could where ever I felt it was necessary.

If you want to understand this code I suggest you read the original paper (You will need a good background in math and physics).

The best reference is the original paper which is still online and contains a very clear explanation of the solver and methods it uses. It also contains a complete C implementation of the solver.

If the original paper even goes offline I have a copy of it [here](https://github.com/topaz1008/canvas-fluid-solver/blob/master/doc/GDC03.pdf).

Other implementations I've looked at while making this. (these links are very old and are dead and/or insecure so proceed at you own risk)

* [http://www.multires.caltech.edu/teaching/demos/java/stablefluids.htm](http://www.multires.caltech.edu/teaching/demos/java/stablefluids.htm)
* ~~[http://blog.inspirit.ru/fluidsolver-as3-port-of-msafluid/](http://blog.inspirit.ru/fluidsolver-as3-port-of-msafluid/)~~
* ~~[http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html](http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html)~~

Running
---------

Just open `index.html` in your browser or see the live demo [HERE](https://topaz1008.github.io/canvas-fluid-solver).

Wishlist
---------
* Add internal grid boundaries.
* Maybe think about considering thinking about a 3D implementation ;)
