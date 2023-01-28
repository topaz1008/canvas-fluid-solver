/**
 * A Simple fluid solver implementation in javascript.
 *
 * Largely based on Jos Stam's paper "Real-Time Fluid Dynamics for Games".
 * @link https://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
 *
 * Simulates the Navier–Stokes equations for incompressible fluids.
 * @link https://en.wikipedia.org/wiki/Navier-Stokes_equations
 *
 * @author Topaz Bar <topaz1008@gmail.com>
 */
export class FluidSolver {

    // Boundaries enumeration.
    static BOUNDARY_NONE = 0;
    static BOUNDARY_LEFT_RIGHT = 1;
    static BOUNDARY_TOP_BOTTOM = 2;

    /**
     * @param n {Number} Number of fluid cells for the simulation grid in each dimension (NxN)
     * @constructor
     */
    constructor(n) {
        this.n = n;

        this.dt = 0.23; // The simulation time-step
        this.diffusion = 0.0004; // The amount of diffusion
        this.viscosity = 0; // The fluid's viscosity

        // Number of iterations to use in the Gauss-Seidel method in linearSolve()
        this.iterations = 15;

        this.doVorticityConfinement = true;
        this.doBuoyancy = true;

        // Two extra cells in each dimension for the boundaries
        this.numOfCells = (n + 2) * (n + 2);

        this.tmp = null; // Scratch space for references swapping

        // This might benefit from using typed arrays like Float32Array in some configuration.
        // But I haven't seen any significant improvement on Chrome because V8 probably does it on its own.

        // Values for current simulation step
        this.u = new Array(this.numOfCells); // Velocity x
        this.v = new Array(this.numOfCells); // Velocity y
        this.d = new Array(this.numOfCells); // Density

        // Values from the last simulation step
        this.uOld = new Array(this.numOfCells);
        this.vOld = new Array(this.numOfCells);
        this.dOld = new Array(this.numOfCells);

        this.curlData = new Array(this.numOfCells); // The cell's curl

        // Initialize everything to zero
        for (let i = 0; i < this.numOfCells; i++) {
            this.d[i] = this.u[i] = this.v[i] = 0;
            this.dOld[i] = this.uOld[i] = this.vOld[i] = 0;
            this.curlData[i] = 0;
        }
    }

    /**
     * Fluid cell indexing helper function.
     * (x | x) is a faster Math.floor(x)
     *
     * For public use.
     *
     * @return {number}
     * @public
     */
    I(i, j) {
        return (i | i) + (this.n + 2) * (j | j);
    }

    /**
     * Density step.
     */
    densityStep() {
        this.#addSource(this.d, this.dOld);

        this.#swapD();
        this.#diffuse(FluidSolver.BOUNDARY_NONE, this.d, this.dOld, this.diffusion);

        this.#swapD();
        this.#advect(FluidSolver.BOUNDARY_NONE, this.d, this.dOld, this.u, this.v);

        // Reset for next step
        for (let i = 0; i < this.numOfCells; i++) {
            this.dOld[i] = 0;
        }
    }

    /**
     * Velocity step.
     */
    velocityStep() {
        this.#addSource(this.u, this.uOld);
        this.#addSource(this.v, this.vOld);

        if (this.doVorticityConfinement) {
            this.#vorticityConfinement(this.uOld, this.vOld);
            this.#addSource(this.u, this.uOld);
            this.#addSource(this.v, this.vOld);
        }

        if (this.doBuoyancy) {
            this.#buoyancy(this.vOld);
            this.#addSource(this.v, this.vOld);
        }

        this.#swapU();
        this.#diffuse(FluidSolver.BOUNDARY_LEFT_RIGHT, this.u, this.uOld, this.viscosity);

        this.#swapV();
        this.#diffuse(FluidSolver.BOUNDARY_TOP_BOTTOM, this.v, this.vOld, this.viscosity);

        this.#project(this.u, this.v, this.uOld, this.vOld);
        this.#swapU();
        this.#swapV();

        this.#advect(FluidSolver.BOUNDARY_LEFT_RIGHT, this.u, this.uOld, this.uOld, this.vOld);
        this.#advect(FluidSolver.BOUNDARY_TOP_BOTTOM, this.v, this.vOld, this.uOld, this.vOld);

        this.#project(this.u, this.v, this.uOld, this.vOld);

        // Reset for next step
        for (let i = 0; i < this.numOfCells; i++) {
            this.uOld[i] = this.vOld[i] = 0;
        }
    }

    /**
     * Resets the density.
     */
    resetDensity() {
        for (let i = 0; i < this.numOfCells; i++) {
            this.d[i] = 0;
        }
    }

    /**
     * Resets the velocity.
     */
    resetVelocity() {
        for (let i = 0; i < this.numOfCells; i++) {
            // Set a small value, so we can render the velocity field
            this.v[i] = this.u[i] = 0.001;
        }
    }

    /**
     * Swap velocity x reference.
     * @private
     */
    #swapU() {
        this.tmp = this.u;
        this.u = this.uOld;
        this.uOld = this.tmp;
    }

    /**
     * Swap velocity y reference.
     * @private
     */
    #swapV() {
        this.tmp = this.v;
        this.v = this.vOld;
        this.vOld = this.tmp;
    }

    /**
     * Swap density reference.
     * @private
     */
    #swapD() {
        this.tmp = this.d;
        this.d = this.dOld;
        this.dOld = this.tmp;
    }

    /**
     * Integrate the density sources.
     *
     * @param x {Array<Number>}
     * @param s {Array<Number>}
     * @private
     */
    #addSource(x, s) {
        for (let i = 0; i < this.numOfCells; i++) {
            x[i] += s[i] * this.dt;
        }
    };

    /**
     * Calculate the curl at cell (i, j)
     * This represents the vortex strength at the cell.
     * Computed as: w = (del x U) where U is the velocity vector at (i, j).
     *
     * @param i Number
     * @param j {Number}
     * @return {Number}
     * @private
     */
    #curl(i, j) {
        const duDy = (this.u[this.I(i, j + 1)] - this.u[this.I(i, j - 1)]) * 0.5,
            dvDx = (this.v[this.I(i + 1, j)] - this.v[this.I(i - 1, j)]) * 0.5;

        return duDy - dvDx;
    }

    /**
     * Calculate the vorticity confinement force for each cell.
     * Fvc = (N x W) where W is the curl at (i, j) and N = del |W| / |del |W||.
     * N is the vector pointing to the vortex center, hence we
     * add force perpendicular to N.
     *
     * @param vcX {Array<Number>}
     * @param vcY {Array<Number>}
     * @private
     */
    #vorticityConfinement(vcX, vcY) {
        // Calculate magnitude of curl(i, j) for each cell
        for (let i = 1; i <= this.n; i++) {
            for (let j = 1; j <= this.n; j++) {
                this.curlData[this.I(i, j)] = Math.abs(this.#curl(i, j));
            }
        }

        for (let i = 1; i <= this.n; i++) {
            for (let j = 1; j <= this.n; j++) {
                // Calculate the derivative of the magnitude (n = del |w|)
                let dx = (this.curlData[this.I(i + 1, j)] - this.curlData[this.I(i - 1, j)]) * 0.5;
                let dy = (this.curlData[this.I(i, j + 1)] - this.curlData[this.I(i, j - 1)]) * 0.5;

                let norm = Math.sqrt((dx * dx) + (dy * dy));
                if (norm === 0) {
                    // Avoid divide by zero
                    norm = 1;
                }

                dx /= norm;
                dy /= norm;

                const v = this.#curl(i, j);

                // N x W
                vcX[this.I(i, j)] = dy * v * -1;
                vcY[this.I(i, j)] = dx * v;
            }
        }
    }

    /**
     * Calculate the buoyancy force for the grid.
     * Fbuoy = -a * d * Y + b * (T - Tamb) * Y where Y = (0,1)
     * The constants a and b are positive with physically meaningful quantities.
     * T is the temperature at the current cell, Tamb is the average temperature of the fluid grid
     *
     * In this simplified implementation we say that the temperature is synonymous with density
     * and because there are no other heat sources we can just use the density field instead of adding a new
     * temperature field.
     *
     * @param buoy {Array<Number>}
     * @private
     */
    #buoyancy(buoy) {
        const a = 0.000625,
            b = 0.025;

        // Sum all temperatures
//    for (let i = 1; i <= this.n; i++) {
//        for (let j = 1; j <= this.n; j++) {
//            tAmb += this.d[this.I(i, j)];
//        }
//    }

        // Sum all temperatures (faster)
        let tAmb = 0;
        for (let i = 0; i < this.d.length; i++) {
            tAmb += this.d[i];
        }

        // Calculate average temperature of the grid
        tAmb /= (this.n * this.n);

        // For each cell compute buoyancy force
        for (let i = 1; i <= this.n; i++) {
            for (let j = 1; j <= this.n; j++) {
                buoy[this.I(i, j)] = a * this.d[this.I(i, j)] + -b * (this.d[this.I(i, j)] - tAmb);
            }
        }
    }

    /**
     * Diffuse the density between neighbouring cells.
     *
     * @param b {Number}
     * @param x {Array<Number>}
     * @param x0 {Array<Number>}
     * @param diffusion {Number}
     * @private
     */
    #diffuse(b, x, x0, diffusion) {
        const a = this.dt * diffusion * this.n * this.n;

        this.#linearSolve(b, x, x0, a, 1 + 4 * a);
    }

    /**
     * The advection step moves the density through the static velocity field.
     * Instead of moving the cells forward in time, we treat the cell's center as a particle
     * and then trace it back in time to look for the 'particles' which end up at the cell's center.
     *
     * @param b {Number}
     * @param d {Array<Number>}
     * @param d0 {Array<Number>}
     * @param u {Array<Number>}
     * @param v {Array<Number>}
     * @private
     */
    #advect(b, d, d0, u, v) {
        const dt0 = this.dt * this.n;
        for (let i = 1; i <= this.n; i++) {
            for (let j = 1; j <= this.n; j++) {
                let x = i - dt0 * u[this.I(i, j)];
                let y = j - dt0 * v[this.I(i, j)];

                if (x < 0.5) x = 0.5;
                if (x > this.n + 0.5) x = this.n + 0.5;

                const i0 = (x | x);
                const i1 = i0 + 1;

                if (y < 0.5) y = 0.5;
                if (y > this.n + 0.5) y = this.n + 0.5;

                const j0 = (y | y);
                const j1 = j0 + 1;
                const s1 = x - i0;
                const s0 = 1 - s1;
                const t1 = y - j0;
                const t0 = 1 - t1;

                d[this.I(i, j)] = s0 * (t0 * d0[this.I(i0, j0)] + t1 * d0[this.I(i0, j1)]) +
                    s1 * (t0 * d0[this.I(i1, j0)] + t1 * d0[this.I(i1, j1)]);
            }
        }

        this.#setBoundary(b, d);
    }

    /**
     * Forces the velocity field to be mass conserving.
     * This step is what actually produces the nice looking swirly vortices.
     *
     * It uses a result called Hodge Decomposition which says that every velocity field is the sum
     * of a mass conserving field, and a gradient field. So we calculate the gradient field, and subtract
     * it from the velocity field to get a mass conserving one.
     * It solves a linear system of equations called Poisson Equation.
     *
     * @param u {Array<Number>}
     * @param v {Array<Number>}
     * @param p {Array<Number>}
     * @param div {Array<Number>}
     * @private
     */
    #project(u, v, p, div) {
        // Calculate the gradient field
        const h = 1.0 / this.n;
        for (let i = 1; i <= this.n; i++) {
            for (let j = 1; j <= this.n; j++) {
                div[this.I(i, j)] = -0.5 * h * (u[this.I(i + 1, j)] - u[this.I(i - 1, j)] +
                    v[this.I(i, j + 1)] - v[this.I(i, j - 1)]);

                p[this.I(i, j)] = 0;
            }
        }

        this.#setBoundary(FluidSolver.BOUNDARY_NONE, div);
        this.#setBoundary(FluidSolver.BOUNDARY_NONE, p);

        // Solve the Poisson equations
        this.#linearSolve(FluidSolver.BOUNDARY_NONE, p, div, 1, 4);

        // Subtract the gradient field from the velocity field to get a mass conserving velocity field.
        for (let i = 1; i <= this.n; i++) {
            for (let j = 1; j <= this.n; j++) {
                u[this.I(i, j)] -= 0.5 * (p[this.I(i + 1, j)] - p[this.I(i - 1, j)]) / h;
                v[this.I(i, j)] -= 0.5 * (p[this.I(i, j + 1)] - p[this.I(i, j - 1)]) / h;
            }
        }

        this.#setBoundary(FluidSolver.BOUNDARY_LEFT_RIGHT, u);
        this.#setBoundary(FluidSolver.BOUNDARY_TOP_BOTTOM, v);
    }

    /**
     * Solve a linear system of equations using Gauss-Seidel method.
     *
     * @param b {Number}
     * @param x {Array<Number>}
     * @param x0 {Array<Number>}
     * @param a {Number}
     * @param c {Number}
     * @private
     */
    #linearSolve(b, x, x0, a, c) {
        const invC = 1.0 / c;

        for (let k = 0; k < this.iterations; k++) {
            for (let i = 1; i <= this.n; i++) {
                for (let j = 1; j <= this.n; j++) {
                    x[this.I(i, j)] = (x0[this.I(i, j)] + a * (x[this.I(i - 1, j)] + x[this.I(i + 1, j)] +
                        x[this.I(i, j - 1)] + x[this.I(i, j + 1)])) * invC;
                }
            }

            this.#setBoundary(b, x);
        }
    }

    /**
     * Set boundary conditions.
     *
     * @param b {Number}
     * @param x {Array<Number>}
     * @private
     */
    #setBoundary(b, x) {
        for (let i = 1; i <= this.n; i++) {
            x[this.I(0, i)] = (b === FluidSolver.BOUNDARY_LEFT_RIGHT) ?
                -x[this.I(1, i)] : x[this.I(1, i)];

            x[this.I(this.n + 1, i)] = (b === FluidSolver.BOUNDARY_LEFT_RIGHT) ?
                -x[this.I(this.n, i)] : x[this.I(this.n, i)];

            x[this.I(i, 0)] = (b === FluidSolver.BOUNDARY_TOP_BOTTOM) ?
                -x[this.I(i, 1)] : x[this.I(i, 1)];

            x[this.I(i, this.n + 1)] = (b === FluidSolver.BOUNDARY_TOP_BOTTOM) ?
                -x[this.I(i, this.n)] : x[this.I(i, this.n)];
        }

        x[this.I(0, 0)] = 0.5 * (x[this.I(1, 0)] + x[this.I(0, 1)]);
        x[this.I(0, this.n + 1)] = 0.5 * (x[this.I(1, this.n + 1)] + x[this.I(0, this.n)]);
        x[this.I(this.n + 1, 0)] = 0.5 * (x[this.I(this.n, 0)] + x[this.I(this.n + 1, 1)]);
        x[this.I(this.n + 1, this.n + 1)] = 0.5 * (x[this.I(this.n, this.n + 1)] + x[this.I(this.n + 1, this.n)]);
    }
}
