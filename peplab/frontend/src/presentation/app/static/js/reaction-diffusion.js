/**
 * Implements GPU-accelerated Gray-Scott reaction-diffusion algorithm
 */
class ReactionDiffusion {
    constructor(canvas) {
        this.canvas = canvas;
        this.ctx = canvas.getContext('2d');
        this.useCPU = true;

        // Add scale factor to reduce resolution
        this.scale = 2; // Each simulation cell will be 4x4 pixels on screen

        // Reaction-diffusion parameters
        this.dA = 0.485;    // Diffusion rate of A
        this.dB = 0.65 - 0.47;    // Diffusion rate of B
        this.feed = 0.055 + 0.0 ; // Feed rate
        this.kill = 0.062 + 0.00; // Kill rate
        this.dt = 0.49;    // Time step

        // Create two buffers for the simulation at reduced resolution
        this.width = Math.floor(canvas.width / this.scale);
        this.height = Math.floor(canvas.height / this.scale);
        this.current = this.createBuffer();
        this.next = this.createBuffer();
        
        // Initialize the simulation state
        this.init();

        // Add steps per frame parameter
        this.stepsPerFrame = 8; // Number of simulation steps per animation frame

        // Add binary mode option
        this.binaryMode = false; // false for gradient, true for binary black/white
        
        // Add invert colors option
        this.invertColors = true; // false for black pattern, true for white pattern
    }

    createBuffer() {
        return {
            a: new Float32Array(this.width * this.height),
            b: new Float32Array(this.width * this.height)
        };
    }

    init() {
        // Fill with chemical A
        for (let i = 0; i < this.width * this.height; i++) {
            this.current.a[i] = 1.0;
            this.current.b[i] = 0.0;
        }

        // Add random seeds of chemical B
        const numSeeds = 10;
        const seedSize = Math.min(21, Math.floor(Math.min(this.width, this.height) * 1));
        
        for (let seed = 0; seed < numSeeds; seed++) {
            const startX = Math.floor(Math.random() * (this.width - seedSize));
            const startY = Math.floor(Math.random() * (this.height - seedSize));
            
            for (let y = 0; y < seedSize; y++) {
                for (let x = 0; x < seedSize; x++) {
                    const i = (startY + y) * this.width + (startX + x);
                    this.current.b[i] = 1.0;
                }
            }
        }
    }

    getIndex(x, y) {
        // Handle wrapping around edges
        x = (x + this.width) % this.width;
        y = (y + this.height) % this.height;
        return y * this.width + x;
    }

    update() {
        // For each pixel, compute the new values
        for (let y = 0; y < this.height; y++) {
            for (let x = 0; x < this.width; x++) {
                const i = this.getIndex(x, y);
                
                // Get the current values
                const a = this.current.a[i];
                const b = this.current.b[i];
                
                // Compute Laplacian
                let laplaceA = 0;
                let laplaceB = 0;
                
                // Check neighboring cells
                const neighbors = [
                    this.getIndex(x-1, y),
                    this.getIndex(x+1, y),
                    this.getIndex(x, y-1),
                    this.getIndex(x, y+1)
                ];
                
                for (const ni of neighbors) {
                    laplaceA += this.current.a[ni];
                    laplaceB += this.current.b[ni];
                }
                
                laplaceA -= 4 * a;
                laplaceB -= 4 * b;
                
                // Reaction-diffusion formula
                const da = this.dA * laplaceA - a * b * b + this.feed * (1 - a);
                const db = this.dB * laplaceB + a * b * b - (this.kill + this.feed) * b;
                
                // Update next state
                this.next.a[i] = Math.max(0, Math.min(1, a + da * this.dt));
                this.next.b[i] = Math.max(0, Math.min(1, b + db * this.dt));
            }
        }
        
        // Swap buffers
        [this.current, this.next] = [this.next, this.current];
    }

    draw() {
        const imageData = this.ctx.createImageData(this.canvas.width, this.canvas.height);
        const data = imageData.data;
        
        // Background color matching card container (rgba(245, 246, 250, 0.9))
        const bgColor = {r: 245, g: 246, b: 250};
        
        // Scale up the simulation to fill the canvas
        for (let y = 0; y < this.canvas.height; y++) {
            for (let x = 0; x < this.canvas.width; x++) {
                // Map canvas coordinates to simulation coordinates
                const simX = Math.floor(x / this.scale);
                const simY = Math.floor(y / this.scale);
                const simIdx = simY * this.width + simX;
                
                // Get the value and apply inversion if needed
                let intensity = this.current.b[simIdx];
                if (!this.invertColors) {
                    intensity = 1 - intensity; // Invert for black pattern
                }
                
                // Apply threshold if in binary mode
                if (this.binaryMode) {
                    intensity = intensity > 0.111 ? 1 : 0;
                }
                
                // Interpolate between background color and black
                const canvasIdx = (y * this.canvas.width + x) * 4;
                data[canvasIdx] = bgColor.r * (1 - intensity);     // Red
                data[canvasIdx + 1] = bgColor.g * (1 - intensity); // Green
                data[canvasIdx + 2] = bgColor.b * (1 - intensity); // Blue
                data[canvasIdx + 3] = 255;                         // Alpha
            }
        }
        
        this.ctx.putImageData(imageData, 0, 0);
    }

    animate() {
        // Run multiple simulation steps per frame
        for (let i = 0; i < this.stepsPerFrame; i++) {
            this.update();
        }
        this.draw();
        requestAnimationFrame(() => this.animate());
    }

    resize(width, height) {
        this.canvas.width = width;
        this.canvas.height = height;
        
        // Update simulation dimensions at reduced resolution
        this.width = Math.floor(width / this.scale);
        this.height = Math.floor(height / this.scale);
        
        // Create new buffers at new size
        this.current = this.createBuffer();
        this.next = this.createBuffer();
        
        // Reinitialize
        this.init();
    }
}

// Initialize when document is loaded
document.addEventListener('DOMContentLoaded', () => {
    const canvas = document.getElementById('reaction-canvas');
    if (!canvas) {
        console.error('Reaction diffusion canvas not found');
        return;
    }

    // Set canvas size
    const setCanvasSize = () => {
        const width = window.innerWidth;
        const height = window.innerHeight;
        canvas.style.width = width + 'px';
        canvas.style.height = height + 'px';
        canvas.width = width;
        canvas.height = height;
        console.log(`Canvas size set to ${width}x${height}`);
    };
    
    setCanvasSize();
    
    const rd = new ReactionDiffusion(canvas);
    rd.animate();
    
    // Debounce resize events
    let resizeTimeout;
    window.addEventListener('resize', () => {
        clearTimeout(resizeTimeout);
        resizeTimeout = setTimeout(() => {
            console.log('Window resized - updating canvas');
            setCanvasSize();
            rd.resize(canvas.width, canvas.height);
        }, 250);
    });
}); 