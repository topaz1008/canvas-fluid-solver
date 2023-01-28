/**
 * A simple particle class.
 *
 * @author Topaz Bar <topaz1008@gmail.com>
 */
export class Particle {
    static TIME_TO_LIVE = 5; // Time to live in seconds

    constructor(x, y) {
        this.x = x; this.y = y;   // Position
        this.vx = 0; this.vy = 0; // Velocity
        this.age = 0;
        this.dead = false;
    }
}
