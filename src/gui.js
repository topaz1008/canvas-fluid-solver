/**
 * A simple wrapper around lil-gui
 * just to keep all this boilerplate stuff out of the main file.
 *
 * @author Topaz Bar <topaz1008@gmail.com>
 */
export class AppGUI {
    static CONTAINER_ELEMENT_ID = 'gui-container';

    #gui = null;
    #appOptions = null;
    #fluidSolver = null;

    constructor(GUI, guiOptions, appOptions) {
        this.#gui = new GUI(guiOptions);
        this.#appOptions = appOptions;
        this.#fluidSolver = appOptions.fluidSolver;
    }

    init() {
        this.#gui.add(this.#fluidSolver, 'dt', 0.05, 0.5).step(0.01).name('Time Step');
        this.#gui.add(this.#fluidSolver, 'iterations', 5, 40).step(1).name('Solver Iterations');
        this.#gui.add(this.#fluidSolver, 'diffusion', 0.0, 0.001).step(0.0001).name('Diffusion');

        const viscosities = {
            None: 0,
            'Very Low': 1 / 100000,
            Low: 1 / 5000,
            High: 1 / 1000
        };
        this.#gui.add(this.#fluidSolver, 'viscosity', viscosities).name('Viscosity');

        this.#gui.add(this.#fluidSolver, 'doVorticityConfinement').name('Vorticity Confinement');
        this.#gui.add(this.#fluidSolver, 'doBuoyancy').name('Buoyancy');

        this.#gui.add(this.#appOptions, 'grayscale').name('Grayscale');
        this.#gui.add(this.#appOptions, 'drawVelocityField').name('Draw Velocity Field');
        this.#gui.add(this.#appOptions, 'drawDensityField').name('Draw Density Field');
        this.#gui.add(this.#appOptions, 'drawParticles').name('Draw Particle Effect');

        this.#gui.add(this.#fluidSolver, 'resetVelocity').name('Reset Velocity');
        this.#gui.add(this.#fluidSolver, 'resetDensity').name('Reset Density');
        this.#gui.add(this.#appOptions, 'resetParticles').name('Reset Particles');

        document.getElementById(AppGUI.CONTAINER_ELEMENT_ID)
            .appendChild(this.#gui.domElement);
    }
}
