/**
 * Demo usage of the FluidSolver class.
 *
 * @author Topaz Bar <topaz1008@gmail.com>
 */
(function (window, document) {
    'use strict';

    var NUM_OF_CELLS = 128, // Number of cells (not including the boundary)
        VIEW_SIZE = 640,    // View size (square)
        FPS = 60;           // Frames per second

    // Check if we're on a mobile device
    var isMobile = /mobile/i.test(window.navigator.userAgent);

    if (isMobile) {
        // Lower settings if we're on mobile
        NUM_OF_CELLS = 64;
        VIEW_SIZE = window.outerWidth - 2; // 1px border on each side
    }

    // setInterval still seems to be faster than this most of the time.
    window.requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame ||
        window.webkitRequestAnimationFrame || window.msRequestAnimationFrame;

    requestAnimationFrame(update);
//    setInterval(update, 1000 / FPS);

    var CELL_SIZE = VIEW_SIZE / NUM_OF_CELLS,  // Size of each cell in pixels
        CELL_SIZE_CEIL = Math.ceil(CELL_SIZE); // Size of each cell in pixels (ceiling)

    /**
     * A simple particle class.
     *
     * @param x {Number}
     * @param y {Number}
     * @constructor
     */
    function Particle(x, y) {
        this.x = x;
        this.y = y;
        this.vx = 0;
        this.vy = 0;
        this.age = 0;
        this.dead = false;
    }

    Particle.TIME_TO_LIVE = 5; // Time to live in seconds

    // Globals
    var canvas = document.getElementById('main'),
        context = canvas.getContext('2d');

    // Create the fluid solver
    var fs = new FluidSolver(NUM_OF_CELLS);
    fs.resetVelocity();

    // We draw the density on a bitmap for performance reasons
    var fdBuffer = context.createImageData(VIEW_SIZE, VIEW_SIZE);

    // Demo app variables
    var lastTime = Date.now(),
        isMouseDown = false,
        oldMouseX = 0,
        oldMouseY = 0,
        particles = [];

    var options = {
        drawVelocityField: false,
        drawDensityField: true,
        drawParticles: true,
        grayscale: false,
        resetParticles: function () { particles.length = 0; }
    };

    // Setup the gui
    var gui = new dat.GUI({
        width: ((isMobile) ? VIEW_SIZE : 350),
        autoPlace: !isMobile
    });

    gui.add(fs, 'dt', 0.05, 0.5).step(0.01).name('Time Step');
    gui.add(fs, 'iterations', 5, 40).step(1).name('Solver Iterations');
    gui.add(fs, 'diffusion', 0.0, 0.001).step(0.0001).name('Diffusion');
    gui.add(fs, 'viscosity', { None: 0.0, Low: 0.0002, High: 0.001 }).name('Viscosity');
    gui.add(fs, 'doVorticityConfinement').name('Vorticity Confinement');
    gui.add(fs, 'doBuoyancy').name('Buoyancy');

    gui.add(options, 'grayscale').name('Grayscale');
    gui.add(options, 'drawVelocityField').name('Draw Velocity Field');
    gui.add(options, 'drawDensityField').name('Draw Density Field');
    gui.add(options, 'drawParticles').name('Draw Particle Effect');

    gui.add(fs, 'resetVelocity').name('Reset Velocity');
    gui.add(fs, 'resetDensity').name('Reset Density');
    gui.add(options, 'resetParticles').name('Reset Particles');

    if (isMobile) {
        // Custom placement on mobile.
        document.getElementById('gui-container').appendChild(gui.domElement);
    }

    // Set render states
    canvas.width = canvas.height = VIEW_SIZE;       // View size
    context.lineWidth = 1;                          // Velocity field line width
    context.strokeStyle = 'rgb(192, 0, 0)';         // Velocity field color
    //context.globalCompositeOperation = 'screen';  // Blend mode

    // Disable smoothing when using floating point pixel values
    context.imageSmoothingEnabled = false;
    context.webkitImageSmoothingEnabled = false;
    context.mozImageSmoothingEnabled = false;

    /**
     * Add event listeners
     */
    document.addEventListener('mouseup', function () { isMouseDown = false; }, false);
    document.addEventListener('mousedown', function () { isMouseDown = true; }, false);

    // Mouse move listener (on the canvas element)
    canvas.addEventListener('mousemove', onMouseMove, false);

    // Touch listeners (on the canvas element)
    if (isMobile) {
        canvas.addEventListener('touchstart', onTouchStart, false);
        canvas.addEventListener('touchend', onTouchEnd, false);
        canvas.addEventListener('touchleave', onTouchEnd, false);
        canvas.addEventListener('touchcancel', onTouchEnd, false);
        canvas.addEventListener('touchmove', onTouchMove, false);
    }

    /**
     * Main mouse move listener
     *
     * @param event {MouseEvent|Object}
     */
    function onMouseMove(event) {
        var k, p;

        var mouseX = event.offsetX,
            mouseY = event.offsetY;

        // Find the cell below the mouse
        var i = (mouseX / VIEW_SIZE) * NUM_OF_CELLS + 1,
            j = (mouseY / VIEW_SIZE) * NUM_OF_CELLS + 1;

        // Dont overflow grid bounds
        if (i > NUM_OF_CELLS || i < 1 || j > NUM_OF_CELLS || j < 1) return;

        // Mouse velocity
        var du = (mouseX - oldMouseX) * 1.5,
            dv = (mouseY - oldMouseY) * 1.5;

        // Add the mouse velocity to cells above, below, to the left, and to the right as well.
        fs.uOld[fs.I(i, j)] = du;
        fs.vOld[fs.I(i, j)] = dv;

        fs.uOld[fs.I(i + 1, j)] = du;
        fs.vOld[fs.I(i + 1, j)] = dv;

        fs.uOld[fs.I(i - 1, j)] = du;
        fs.vOld[fs.I(i - 1, j)] = dv;

        fs.uOld[fs.I(i, j + 1)] = du;
        fs.vOld[fs.I(i, j + 1)] = dv;

        fs.uOld[fs.I(i, j - 1)] = du;
        fs.vOld[fs.I(i, j - 1)] = dv;

        if (isMouseDown) {
            // If holding down the mouse, add density to the cell below the mouse
            fs.dOld[fs.I(i, j)] = 50;
        }

        if (isMouseDown && options.drawParticles) {
            // Add particles
            for (k = 0; k < 5; k++) {
                p = new Particle(mouseX + getRandom(-50, 50), mouseY + getRandom(-50, 50));

                p.vx = du;
                p.vy = dv;
                particles.push(p);
            }
        }

        // Save current mouse position for next frame
        oldMouseX = mouseX;
        oldMouseY = mouseY;

    } // End onMouseMove()

    /**
     * Touch listeners.
     *
     * preventDefault() prevents the mouse events from being dispatched.
     */
    function onTouchStart(event) { event.preventDefault(); isMouseDown = true; }
    function onTouchEnd(event) { event.preventDefault(); isMouseDown = false; }

    /**
     * Touch move listener.
     * just passes the call to onMouseMove with the correct coordinates.
     *
     * @param event {*} The TouchEvent
     */
    function onTouchMove(event) {
        event.preventDefault();

        //noinspection JSUnresolvedVariable
        var touches = event.touches;
        if (touches && touches.length > 0) {
            onMouseMove({
                // Offset by the layer's (canvas) position.
                offsetX: touches[0].screenX + event.layerX,
                offsetY: touches[0].screenY + event.layerY - 50 // 50px top margin
            });
        }
    }

    /**
     * Draw the simulation grid.
     */
    function drawGrid() {
        var i;

        context.lineWidth = 1;
        context.strokeStyle = 'rgb(255, 255, 255)';
        context.beginPath();

        // Vertical
        for (i = 0; i <= VIEW_SIZE; i += CELL_SIZE_CEIL) {
            context.moveTo(i, 0);
            context.lineTo(i, VIEW_SIZE);
        }

        // Horizontal
        for (i = 0; i <= VIEW_SIZE; i += CELL_SIZE_CEIL) {
            context.moveTo(0, i);
            context.lineTo(VIEW_SIZE, i);
        }

        context.stroke();
    }

    /**
     * Update loop
     */
    function update(/*time*/) {
        var i, j, k, l, m, dx, dy, color, cellIndex, deltaTime,
            r, g, b, u, v, density, pxIdx, pxX, pxY,
            invMaxColor = 1.0 / 255;

        deltaTime = (Date.now() - lastTime) / 1000;

        // Step the fluid simulation
        fs.velocityStep();
        fs.densityStep();

        // Clear the canvas
        context.clearRect(0, 0, VIEW_SIZE, VIEW_SIZE);

        // Draw the last frame's buffer and clear for drawing the current.
        context.putImageData(fdBuffer, 0, 0);
        clearImageData(fdBuffer);

        //drawGrid();

        if (options.drawVelocityField) {
            // Call once per frame
            context.lineWidth = 1;
            context.strokeStyle = 'rgb(192, 0, 0)';
            context.beginPath();
        }

        // Render fluid
        for (i = 1; i <= NUM_OF_CELLS; i++) {
            // The x position of current cell
            dx = (i - 0.5) * CELL_SIZE;

            for (j = 1; j <= NUM_OF_CELLS; j++) {
                // The y position of current cell
                dy = (j - 0.5) * CELL_SIZE;

                cellIndex = i + (NUM_OF_CELLS + 2) * j;

                // Draw density
                density = fs.d[cellIndex];
                if (options.drawDensityField && density > 0) {
                    color = density * 255;

                    // fdBuffer.data is actually a Uint8ClampedArray so there is no need to manually clamp color values
                    //if (color < 0) color = 0;
                    //if (color > 255) color = 255;

                    r = color;
                    g = ((options.grayscale) ? color : color * dx * invMaxColor);
                    b = ((options.grayscale) ? color : color * dy * invMaxColor);

                    // Draw the cell on an image for performance reasons
                    for (l = 0; l < CELL_SIZE_CEIL; l++) {
                        for (m = 0; m < CELL_SIZE_CEIL; m++) {
                            pxX = (i - 1) * CELL_SIZE + l;
                            pxY = (j - 1) * CELL_SIZE + m;
                            pxIdx = ((pxX | pxX) + (pxY | pxY) * VIEW_SIZE) * 4;

                            fdBuffer.data[pxIdx    ] = r;
                            fdBuffer.data[pxIdx + 1] = g;
                            fdBuffer.data[pxIdx + 2] = b;
                            fdBuffer.data[pxIdx + 3] = 255;
                        }
                    }
                }

                // Draw velocity field ?
                if (options.drawVelocityField && (i % 2) === 0 && (j % 2) === 0) {
                    u = fs.u[cellIndex] * 50;
                    v = fs.v[cellIndex] * 50;

                    context.moveTo(dx, dy);
                    context.lineTo(dx + u, dy + v);
                }

            } // End for all cells in the y direction

        } // End for all cells in the x direction

        if (options.drawVelocityField) {
            // Call once per frame
            context.stroke();
        }

        // Update and render particles
        var x0, y0, p, alpha, lastAlpha = 0, particlesLength = particles.length;
        if (options.drawParticles) {
            context.lineWidth = 2;
            context.strokeStyle = 'rgb(255, 255, 255)';
            context.beginPath();

            for (k = 0; k < particlesLength; k++) {
                p = particles[k];

                p.age += deltaTime;
                alpha = (1 - (p.age / Particle.TIME_TO_LIVE));
                if ((alpha < 0.01) ||
                    (p.age >= Particle.TIME_TO_LIVE) ||
                    (p.x <= 0 || p.x >= VIEW_SIZE || p.y <= 0 || p.y >= VIEW_SIZE)) {
                    p.dead = true;

                } else {
                    x0 = (p.x / VIEW_SIZE) * NUM_OF_CELLS + 2;
                    y0 = (p.y / VIEW_SIZE) * NUM_OF_CELLS + 2;

                    cellIndex = fs.I(x0, y0);

                    p.vx = fs.u[cellIndex] * 50;
                    p.vy = fs.v[cellIndex] * 50;

                    p.x += p.vx;
                    p.y += p.vy;

                    if (Math.abs(alpha - lastAlpha) > 0.001) {
                        // Only change stroke style if the alpha changed to save on render state changes.
                        context.strokeStyle = 'rgba(255, 255, 255, ' + alpha + ')';
                        lastAlpha = alpha;
                    }

                    context.moveTo(p.x, p.y);
                    context.lineTo(p.x + p.vx, p.y + p.vy);
                }

                if (p.dead) {
                    // Remove dead particles, and update the length manually
                    particles.splice(k, 1);
                    particlesLength = particles.length;
                }

            } // End for all particles

            context.stroke();

        } // End if drawParticles

        // lastTime is now
        lastTime = Date.now();

        requestAnimationFrame(update);

    } // End update()

    /**
     * @param min {Number}
     * @param max {Number}
     * @returns {Number}
     */
    function getRandom(min, max) {
        return min + Math.random() * ((max + 1) - min);
    }

    /**
     * Clears all the pixels on the image data.
     *
     * @param image {ImageData}
     */
    function clearImageData(image) {
        var i, length = image.data.length;
        for (i = 0; i < length; i++) {
            image.data[i] = 0;
        }
    }

})(window, document);
