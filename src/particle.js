/**
 * A simple particle class.
 */
export class Particle {
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