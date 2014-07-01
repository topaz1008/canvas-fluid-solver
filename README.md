HTML5 Canvas Fluid Solver
==========================

A Simple fluid solver implementation in javascript. A live demo can be found [here](http://topaz1008.github.io/canvas-fluid-solver).

The demo uses [dat.gui](https://github.com/dataarts/dat.gui) for the GUI (included in the `vendors` folder).

Simulates the [Navier–Stokes](http://en.wikipedia.org/wiki/Navier-Stokes_equations) equations for incompressible fluids.

Largely based on Jos Stam's paper [Real-Time Fluid Dynamics for Games](http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf).

The density field is rendered on an off-screen bitmap and blit'ed to the canvas for performance reasons.

The code is commented as best I could where ever I felt it was necessary.

If you want to understand this code I suggest you read the original paper (You will need a good background in math and physics).

Other implementations I've looked at while making this.

* [http://www.multires.caltech.edu/teaching/demos/java/stablefluids.htm](http://www.multires.caltech.edu/teaching/demos/java/stablefluids.htm)
* [http://blog.inspirit.ru/fluidsolver-as3-port-of-msafluid/](http://blog.inspirit.ru/fluidsolver-as3-port-of-msafluid/)
* [http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html](http://http.developer.nvidia.com/GPUGems/gpugems_ch38.html)

Running
---------

Just open `index.html` in your browser.
