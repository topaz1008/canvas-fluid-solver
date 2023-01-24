/**
 * Demo usage of the FluidSolver class.
 *
 * @author Topaz Bar <topaz1008@gmail.com>
 */
import { FluidSolver } from './fluidsolver.js';

const NUM_OF_CELLS = 128, // Number of cells (not including the boundary)
    VIEW_SIZE = 640;    // View size (square)

const CELL_SIZE = VIEW_SIZE / NUM_OF_CELLS,  // Size of each cell in pixels
    CELL_SIZE_CEIL = Math.ceil(CELL_SIZE); // Size of each cell in pixels (ceiling)

const isMobile = window.navigator.userAgentData.mobile;
console.log(`isMobile: ${isMobile}`);

/**
 * A simple particle class.
 */
class Particle {
    static TIME_TO_LIVE = 5; // Time to live in seconds

    constructor(x, y) {
        this.x = x;
        this.y = y;
        this.vx = 0;
        this.vy = 0;
        this.age = 0;
        this.dead = false;
    }
}

// Globals
const canvas = document.getElementById('main-canvas'),
    context = canvas.getContext('2d');

// Create the fluid solver
const fs = new FluidSolver(NUM_OF_CELLS);
fs.resetVelocity();

// We draw the density on a bitmap for performance reasons
const fdBuffer = context.createImageData(VIEW_SIZE, VIEW_SIZE);

// Demo app variables
let lastTime = Date.now(),
    isMouseDown = true,
    oldMouseX = 0,
    oldMouseY = 0,
    particles = [];

const options = {
    drawVelocityField: false,
    drawDensityField: true,
    drawParticles: false,
    grayscale: true,
    resetParticles: function () { particles.length = 0; }
};

// Set up the gui
const gui = new dat.GUI({
    width: 400,
    autoPlace: false
});

gui.add(fs, 'dt', 0.05, 0.5).step(0.01).name('Time Step');
gui.add(fs, 'iterations', 5, 40).step(1).name('Solver Iterations');
gui.add(fs, 'diffusion', 0.0, 0.001).step(0.0001).name('Diffusion');
gui.add(fs, 'viscosity', { None: 0.0, 'Very Low': 0.0003, Low: 0.0002, High: 0.001 }).name('Viscosity');
gui.add(fs, 'doVorticityConfinement').name('Vorticity Confinement');
gui.add(fs, 'doBuoyancy').name('Buoyancy');

gui.add(options, 'grayscale').name('Grayscale');
gui.add(options, 'drawVelocityField').name('Draw Velocity Field');
gui.add(options, 'drawDensityField').name('Draw Density Field');
gui.add(options, 'drawParticles').name('Draw Particle Effect');

gui.add(fs, 'resetVelocity').name('Reset Velocity');
gui.add(fs, 'resetDensity').name('Reset Density');
gui.add(options, 'resetParticles').name('Reset Particles');

// Attach gui to dom
document.getElementById('gui-container').appendChild(gui.domElement);

// Set render states
canvas.width = canvas.height = VIEW_SIZE;       // View size
context.lineWidth = 1;                          // Velocity field line width
context.strokeStyle = 'rgb(192, 0, 0)';         // Velocity field color
// context.globalCompositeOperation = 'screen';  // Blend mode

// Disable smoothing when using floating point pixel values
context.imageSmoothingEnabled = false;

/**
 * Add event listeners
 */
document.addEventListener('mouseup', () => { isMouseDown = false; }, false);
document.addEventListener('mousedown', () => { isMouseDown = true; }, false);

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
 * @param e {MouseEvent|Object}
 */
function onMouseMove(e) {
    const mouseX = e.offsetX,
        mouseY = e.offsetY;

    // console.log(`onMouseMove: x: ${mouseX} y: ${mouseY}`);

    // Find the cell below the mouse
    const i = (mouseX / VIEW_SIZE) * NUM_OF_CELLS + 1,
        j = (mouseY / VIEW_SIZE) * NUM_OF_CELLS + 1;

    // Don't overflow grid bounds
    if (i > NUM_OF_CELLS || i < 1 || j > NUM_OF_CELLS || j < 1) return;

    // Mouse velocity
    const du = (mouseX - oldMouseX) * 0.5,
        dv = (mouseY - oldMouseY) * 0.5;

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
        fs.dOld[fs.I(i, j)] = 150;
    }

    if (isMouseDown && options.drawParticles) {
        // Add particles
        for (let k = 0; k < 5; k++) {
            const p = new Particle(mouseX + getRandom(-50, 50), mouseY + getRandom(-50, 50));

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
function onTouchStart(e) { e.preventDefault(); isMouseDown = true; }
function onTouchEnd(e) { e.preventDefault(); isMouseDown = false; }

/**
 * Touch move listener.
 * just passes the call to onMouseMove with the correct coordinates.
 *
 * @param e {*} The TouchEvent
 */
function onTouchMove(e) {
    e.preventDefault();

    //noinspection JSUnresolvedVariable
    const touches = e.touches;
    if (touches && touches.length > 0) {
        // Offset by the layer's (canvas) position.
        const x = touches[0].screenX;
        const y = touches[0].screenY;
        console.log(`onTouchMove: x: ${x} y: ${y}`);
        console.log(touches[0]);
        onMouseMove({ offsetX: x, offsetY: y });
    }
}

//// End event listeners

/**
 * Update loop
 */
function update(/*time*/) {
    const invMaxColor = 1.0 / 255;

    const deltaTime = (Date.now() - lastTime) / 1000;

    // Step the fluid simulation
    fs.velocityStep();
    fs.densityStep();

    // Clear the canvas
    context.clearRect(0, 0, VIEW_SIZE, VIEW_SIZE);

    // Draw the last frame's buffer and clear for drawing the current.
    context.putImageData(fdBuffer, 0, 0);
    clearImageData(fdBuffer);

    // drawGrid();

    if (options.drawVelocityField) {
        // Call once per frame
        context.lineWidth = 1;
        context.strokeStyle = 'rgb(192, 0, 0)';
        context.beginPath();
    }

    // Render fluid
    for (let i = 1; i <= NUM_OF_CELLS; i++) {
        // The x position of current cell
        const dx = (i - 0.5) * CELL_SIZE;

        for (let j = 1; j <= NUM_OF_CELLS; j++) {
            // The y position of current cell
            const dy = (j - 0.5) * CELL_SIZE;

            const cellIndex = i + (NUM_OF_CELLS + 2) * j;

            // Draw density
            const density = fs.d[cellIndex];
            if (options.drawDensityField && density > 0) {
                const color = density * 255;

                // fdBuffer.data is actually a Uint8ClampedArray so there is no need to manually clamp color values
                //if (color < 0) color = 0;
                //if (color > 255) color = 255;

                const r = color;
                const g = ((options.grayscale) ? color : color * dx * invMaxColor);
                const b = ((options.grayscale) ? color : color * dy * invMaxColor);

                // Draw the cell on an image for performance reasons
                for (let l = 0; l < CELL_SIZE_CEIL; l++) {
                    for (let m = 0; m < CELL_SIZE_CEIL; m++) {
                        const pxX = (i - 1) * CELL_SIZE + l;
                        const pxY = (j - 1) * CELL_SIZE + m;
                        const pxIdx = ((pxX | pxX) + (pxY | pxY) * VIEW_SIZE) * 4;

                        fdBuffer.data[pxIdx    ] = r;
                        fdBuffer.data[pxIdx + 1] = g;
                        fdBuffer.data[pxIdx + 2] = b;
                        fdBuffer.data[pxIdx + 3] = 255;
                    }
                }
            }

            // Draw velocity field ?
            if (options.drawVelocityField && (i % 2) === 0 && (j % 2) === 0) {
                const u = fs.u[cellIndex] * 50;
                const v = fs.v[cellIndex] * 50;

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
    let lastAlpha = 0,
        particlesLength = particles.length;
    if (options.drawParticles) {
        context.lineWidth = 2;
        context.strokeStyle = 'rgb(255, 255, 255)';
        context.beginPath();

        for (let k = 0; k < particlesLength; k++) {
            const p = particles[k];

            p.age += deltaTime;
            const alpha = (1 - (p.age / Particle.TIME_TO_LIVE));
            if ((alpha < 0.01) ||
                (p.age >= Particle.TIME_TO_LIVE) ||
                (p.x <= 0 || p.x >= VIEW_SIZE || p.y <= 0 || p.y >= VIEW_SIZE)) {
                p.dead = true;

            } else {
                const x0 = (p.x / VIEW_SIZE) * NUM_OF_CELLS + 2;
                const y0 = (p.y / VIEW_SIZE) * NUM_OF_CELLS + 2;

                const cellIndex = fs.I(x0, y0);

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

// Start app
update();

/**
 * @param min {Number}
 * @param max {Number}
 * @returns {Number}
 */
function getRandom(min, max) {
    return min + Math.random() * (max - min);
}

/**
 * Clears all the pixels on the image data.
 *
 * @param image {ImageData}
 */
function clearImageData(image) {
    for (let i = 0; i < image.data.length; i++) {
        if ((i % 4) === 0) {
            image.data[i] = 255;

        } else {
            image.data[i] = 0;
        }
    }
}

/**
 * Draw the simulation grid.
 */
function drawGrid() {
    context.lineWidth = 1;
    context.strokeStyle = 'rgb(255, 255, 255)';
    context.beginPath();

    // Vertical
    for (let i = 0; i <= VIEW_SIZE; i += CELL_SIZE_CEIL) {
        context.moveTo(i, 0);
        context.lineTo(i, VIEW_SIZE);
    }

    // Horizontal
    for (let i = 0; i <= VIEW_SIZE; i += CELL_SIZE_CEIL) {
        context.moveTo(0, i);
        context.lineTo(VIEW_SIZE, i);
    }

    context.stroke();
}
