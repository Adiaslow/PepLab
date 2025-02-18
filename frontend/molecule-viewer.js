class MoleculeViewer {
    constructor(numMolecules = 5) {
        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(
            75,
            window.innerWidth / window.innerHeight,
            1,
            200,
        );

        // Enhanced renderer settings
        this.renderer = new THREE.WebGLRenderer({
            alpha: true,
            antialias: true,
            logarithmicDepthBuffer: true,
            precision: "highp",
            powerPreference: "high-performance",
            stencil: false,
        });

        // Enable proper pixel ratio handling
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

        this.atomColors = {
            C: new THREE.Color(0x222222),
            N: new THREE.Color(0x2e50f0),
            O: new THREE.Color(0xf02e2e),
            S: new THREE.Color(0xf0f02e),
            H: new THREE.Color(0xfcfcfc),
            P: new THREE.Color(0xffa500),
        };

        // Van der Waals radii (Å)
        this.atomicRadii = {
            C: 0.77,
            N: 0.75,
            O: 0.73,
            S: 1.02,
            H: 0.37,
            P: 1.05,
        };

        // Silhouette settings
        this.silhouetteScale = 1.05;
        this.silhouetteColor = 0x000000;
        this.silhouetteOpacity = 1.0;

        this.molecules = [];
        this.numMolecules = numMolecules;
        this.conformers = [];

        // Background and particle settings
        this.particles = [];
        this.numParticles = 200;
        this.particleGroups = {
            near: { z: [-50, 0], speed: 0.005, size: 0.15 },
            mid: { z: [-100, -50], speed: 0.003, size: 0.1 },
            far: { z: [-150, -100], speed: 0.001, size: 0.05 },
        };

        this.waterMolecules = [];
        this.numWaterMolecules = 30; // Adjust this number as needed

        // Water molecule settings
        this.waterScale = 0.8; // Makes water molecules slightly smaller than main molecules
        this.waterVelocityFactor = 0.01; // Makes water molecules move a bit faster

        // Enhanced collision settings
        this.collisionElasticity = 0.8; // How bouncy collisions are (0-1)
        this.boundingSphereRadius = 3; // Radius for collision detection
        this.collisionDamping = 0.98; // Reduces velocity after collision
        this.minVelocity = 0.0001; // Minimum velocity threshold
        this.maxVelocity = 0.1; // Maximum velocity threshold
        this.rotationTransfer = 0.0; // How much collision affects rotation

        // Collection of all molecular objects for collision detection
        this.allMolecules = [];

        this.init();
        this.createBackgroundParticles();
    }

    init() {
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        document.body.insertBefore(
            this.renderer.domElement,
            document.body.firstChild,
        );
        this.renderer.domElement.style.position = "fixed";
        this.renderer.domElement.style.top = "0";
        this.renderer.domElement.style.left = "0";
        this.renderer.domElement.style.zIndex = "-1";

        // Adjusted camera position
        this.camera.position.z = 30;

        // Enhanced lighting setup
        const ambientLight = new THREE.AmbientLight(0xffffff, 1.5);
        this.scene.add(ambientLight);

        const frontLight = new THREE.DirectionalLight(0xffffff, 1.2);
        frontLight.position.set(2, 2, 2);
        this.scene.add(frontLight);

        const backLight = new THREE.DirectionalLight(0xffffff, 0.3);
        backLight.position.set(-2, -2, -2);
        this.scene.add(backLight);

        // Add soft point lights for better depth perception
        const pointLight1 = new THREE.PointLight(0xffffff, 0.5);
        pointLight1.position.set(10, 10, 10);
        this.scene.add(pointLight1);

        const pointLight2 = new THREE.PointLight(0xffffff, 0.3);
        pointLight2.position.set(-10, -10, -10);
        this.scene.add(pointLight2);

        // Adjusted fog
        this.scene.fog = new THREE.FogExp2(0x777777, 0.005);

        window.addEventListener("resize", () => this.onWindowResize(), false);
    }

    createBackgroundParticles() {
        const particleGeometry = new THREE.BufferGeometry();
        const positions = [];
        const velocities = [];
        const sizes = [];

        Object.entries(this.particleGroups).forEach(([depth, settings]) => {
            const particlesInGroup = Math.floor(this.numParticles / 3);

            for (let i = 0; i < particlesInGroup; i++) {
                // Random position
                positions.push(
                    (Math.random() - 0.5) * 200, // x
                    (Math.random() - 0.5) * 200, // y
                    Math.random() * (settings.z[1] - settings.z[0]) +
                        settings.z[0], // z
                );

                // Velocity (mainly upward with slight variation)
                velocities.push(
                    (Math.random() - 0.5) * settings.speed,
                    (Math.random() - 0.5) * settings.speed,
                    0,
                );

                // Size
                sizes.push(settings.size);
            }
        });

        particleGeometry.setAttribute(
            "position",
            new THREE.Float32BufferAttribute(positions, 3),
        );
        particleGeometry.setAttribute(
            "velocity",
            new THREE.Float32BufferAttribute(velocities, 3),
        );
        particleGeometry.setAttribute(
            "size",
            new THREE.Float32BufferAttribute(sizes, 1),
        );

        const particleMaterial = new THREE.PointsMaterial({
            color: 0xaaaaaa,
            size: 0.1,
            transparent: true,
            opacity: 0.5,
            sizeAttenuation: true,
            // blending: THREE.AdditiveBlending,
        });

        this.particleSystem = new THREE.Points(
            particleGeometry,
            particleMaterial,
        );
        this.scene.add(this.particleSystem);
    }

    updateParticles() {
        const positions =
            this.particleSystem.geometry.attributes.position.array;
        const velocities =
            this.particleSystem.geometry.attributes.velocity.array;

        for (let i = 0; i < positions.length; i += 3) {
            // Update positions
            positions[i] += velocities[i]; // x
            positions[i + 1] += velocities[i + 1]; // y

            // Reset particles that move too far
            if (positions[i + 1] > 100) {
                positions[i + 1] = -100;
                positions[i] = (Math.random() - 0.5) * 200;
            }
        }

        this.particleSystem.geometry.attributes.position.needsUpdate = true;
    }

    parsePDB(pdbContent) {
        const atoms = [];
        const bonds = new Set(); // To store unique bonds
        let currentModel = null;

        // Split into lines and process each line
        const lines = pdbContent.split("\n");

        for (const line of lines) {
            if (line.startsWith("MODEL")) {
                currentModel = parseInt(line.slice(10, 14).trim());
                if (!this.conformers[currentModel]) {
                    this.conformers[currentModel] = {
                        atoms: [],
                        bonds: new Set(),
                    };
                }
            } else if (line.startsWith("ATOM") || line.startsWith("HETATM")) {
                const atom = {
                    serial: parseInt(line.slice(6, 11)),
                    name: line.slice(12, 16).trim(),
                    resName: line.slice(17, 20).trim(),
                    chainID: line.slice(21, 22),
                    resSeq: parseInt(line.slice(22, 26)),
                    x: parseFloat(line.slice(30, 38)),
                    y: parseFloat(line.slice(38, 46)),
                    z: parseFloat(line.slice(46, 54)),
                    element: line.slice(76, 78).trim(),
                };

                if (currentModel !== null) {
                    this.conformers[currentModel].atoms.push(atom);
                } else {
                    atoms.push(atom);
                }
            } else if (line.startsWith("CONECT")) {
                const numbers = line
                    .slice(6)
                    .trim()
                    .split(/\s+/)
                    .map((n) => parseInt(n));
                const atomIndex = numbers[0];

                for (let i = 1; i < numbers.length; i++) {
                    const bond = [
                        Math.min(atomIndex, numbers[i]),
                        Math.max(atomIndex, numbers[i]),
                    ];
                    const bondString = bond.join("-");

                    if (currentModel !== null) {
                        this.conformers[currentModel].bonds.add(bondString);
                    } else {
                        bonds.add(bondString);
                    }
                }
            }
        }

        return { atoms, bonds: Array.from(bonds) };
    }

    async loadPDBFile(url) {
        try {
            const response = await fetch(url);
            const pdbContent = await response.text();
            const structure = this.parsePDB(pdbContent);
            await this.createMolecules(structure);
        } catch (error) {
            console.error("Error loading PDB file:", error);
        }
    }

    getDefaultStructure() {
        return {
            atoms: [
                // Backbone atoms
                { serial: 1, element: "N", x: 0, y: 0, z: 0 },
                { serial: 2, element: "C", x: 1.47, y: 0, z: 0 },
                { serial: 3, element: "C", x: 1.47, y: 1.54, z: 0 },
                { serial: 4, element: "O", x: 0.74, y: 2.31, z: 0 },
                { serial: 5, element: "C", x: 2.94, y: -0.3, z: 0 },
                { serial: 6, element: "S", x: 3.68, y: -1.92, z: 0 },
                // Additional backbone atoms
                { serial: 7, element: "N", x: 2.7, y: 1.85, z: 0 },
                { serial: 8, element: "C", x: 2.94, y: 3.24, z: 0 },
                { serial: 9, element: "C", x: 4.41, y: 3.54, z: 0 },
                { serial: 10, element: "O", x: 4.89, y: 4.62, z: 0 },
            ],
            bonds: [
                // Define connections between atoms
                "1-2",
                "2-3",
                "3-4",
                "2-5",
                "5-6",
                "3-7",
                "7-8",
                "8-9",
                "9-10",
            ],
        };
    }

    createAtom(atom, moleculeGroup, scale = 1.0) {
        const element = atom.element || "C";
        const radius = (this.atomicRadii[element] || 0.77) * scale;

        // Higher quality geometry
        const geometry = new THREE.SphereGeometry(radius, 32, 32);
        const material = new THREE.MeshPhysicalMaterial({
            color: this.atomColors[element] || 0x808080,
            metalness: 0.1,
            roughness: 0.8,
            clearcoat: 0.3,
            clearcoatRoughness: 0.25,
            envMapIntensity: 1.0,
            dithering: true,
        });

        const mesh = new THREE.Mesh(geometry, material);
        mesh.position.set(atom.x * scale, atom.y * scale, atom.z * scale);
        mesh.userData.atom = atom;

        // Improved silhouette with better depth handling
        const silhouetteGeometry = new THREE.SphereGeometry(
            radius * this.silhouetteScale,
            32,
            32,
        );
        const silhouetteMaterial = new THREE.MeshBasicMaterial({
            color: this.silhouetteColor,
            transparent: true,
            opacity: this.silhouetteOpacity,
            side: THREE.BackSide,
            depthWrite: false,
        });

        const silhouette = new THREE.Mesh(
            silhouetteGeometry,
            silhouetteMaterial,
        );
        silhouette.position.copy(mesh.position);
        silhouette.renderOrder = -1;

        moleculeGroup.add(silhouette);
        moleculeGroup.add(mesh);

        return mesh;
    }

    createBond(atom1, atom2, moleculeGroup, scale = 1) {
        const start = new THREE.Vector3(atom1.x, atom1.y, atom1.z);
        const end = new THREE.Vector3(atom2.x, atom2.y, atom2.z);

        const direction = new THREE.Vector3()
            .subVectors(end, start)
            .normalize();
        const distance = start.distanceTo(end);

        // Create main bond
        const bondGeometry = new THREE.CylinderGeometry(
            0.1 * scale,
            0.1 * scale,
            distance * scale,
            0,
        );
        const bondMaterial = new THREE.MeshPhongMaterial({
            color: 0x808080,
            shininess: 30,
        });

        // Create silhouette for bond
        const silhouetteGeometry = new THREE.CylinderGeometry(
            0.15 * scale,
            0.15 * scale,
            distance * scale,
            8,
        );
        const silhouetteMaterial = new THREE.MeshBasicMaterial({
            color: this.silhouetteColor,
            transparent: true,
            opacity: this.silhouetteOpacity,
            side: THREE.BackSide,
        });

        const bond = new THREE.Mesh(bondGeometry, bondMaterial);
        const silhouette = new THREE.Mesh(
            silhouetteGeometry,
            silhouetteMaterial,
        );

        // Position and rotate both meshes
        [bond, silhouette].forEach((mesh) => {
            mesh.position.copy(start.multiplyScalar(scale));
            mesh.position.add(direction.multiplyScalar((distance * scale) / 2));

            const quaternion = new THREE.Quaternion();
            quaternion.setFromUnitVectors(
                new THREE.Vector3(0, 1, 0),
                direction,
            );
            mesh.setRotationFromQuaternion(quaternion);

            moleculeGroup.add(mesh);
        });
    }

    getWaterStructure() {
        return {
            atoms: [
                // Oxygen atom at center
                { serial: 1, element: "O", x: 0, y: 0, z: 0 },
                // Hydrogen atoms with tetrahedral angle (104.5°)
                { serial: 2, element: "H", x: 0.757, y: 0.586, z: 0 },
                { serial: 3, element: "H", x: -0.757, y: 0.586, z: 0 },
            ],
            bonds: [
                "1-2", // O-H bond
                "1-3", // O-H bond
            ],
        };
    }

    createWaterMolecules() {
        const waterStructure = this.getWaterStructure();
        for (let i = 0; i < this.numWaterMolecules; i++) {
            const waterMolecule = this.createMoleculeFromStructure(
                waterStructure,
                this.waterScale,
            );

            // Customize water molecule properties
            waterMolecule.position.set(
                (Math.random() - 0.5) * 100,
                (Math.random() - 0.5) * 50,
                (Math.random() - 0.5) * 100,
            );

            waterMolecule.velocity = new THREE.Vector3(
                (Math.random() - 0.5) * this.waterVelocityFactor,
                (Math.random() - 0.5) * this.waterVelocityFactor,
                (Math.random() - 0.5) * this.waterVelocityFactor,
            );

            // Add some wobble to the rotation
            waterMolecule.rotationVelocity = new THREE.Vector3(
                (Math.random() - 0.5) * 0.002,
                (Math.random() - 0.5) * 0.002,
                (Math.random() - 0.5) * 0.002,
            );

            this.waterMolecules.push(waterMolecule);
            this.scene.add(waterMolecule);
        }
    }

    createMoleculeFromStructure(structure, scale = 1, isWater = false) {
        const moleculeGroup = new THREE.Group();
        const atomMap = new Map();

        // Create atoms and bonds as before
        structure.atoms.forEach((atom) => {
            const mesh = this.createAtom(atom, moleculeGroup, scale);
            atomMap.set(atom.serial, atom);
        });

        structure.bonds.forEach((bondString) => {
            const [atom1Serial, atom2Serial] = bondString
                .split("-")
                .map((n) => parseInt(n));
            const atom1 = atomMap.get(atom1Serial);
            const atom2 = atomMap.get(atom2Serial);
            /*
            if (atom1 && atom2) {
                this.createBond(atom1, atom2, moleculeGroup, scale);
            }
             */
        });

        // Add bounding sphere for collision detection
        const boundingSphere = new THREE.Sphere(
            new THREE.Vector3(),
            this.boundingSphereRadius * scale,
        );
        moleculeGroup.boundingSphere = boundingSphere;

        // Set initial position and velocity
        moleculeGroup.position.set(
            (Math.random() - 0.5) * 128,
            (Math.random() - 0.5) * 64,
            (Math.random() - 0.5) * 128,
        );

        const velocityFactor = isWater ? this.waterVelocityFactor : 0.001;
        moleculeGroup.velocity = new THREE.Vector3(
            (Math.random() - 0.5) * velocityFactor,
            (Math.random() - 0.5) * velocityFactor,
            (Math.random() - 0.5) * velocityFactor,
        );

        moleculeGroup.rotation.set(
            Math.random() * Math.PI,
            Math.random() * Math.PI,
            Math.random() * Math.PI,
        );

        moleculeGroup.rotationVelocity = new THREE.Vector3(
            (Math.random() - 0.5) * 0.001,
            (Math.random() - 0.5) * 0.001,
            (Math.random() - 0.5) * 0.001,
        );

        moleculeGroup.mass = isWater ? 1 : 2; // Water molecules are lighter
        moleculeGroup.isWater = isWater;

        this.scene.add(moleculeGroup);
        this.allMolecules.push(moleculeGroup);
        return moleculeGroup;
    }

    checkCollisions() {
        for (let i = 0; i < this.allMolecules.length; i++) {
            const mol1 = this.allMolecules[i];

            // Update bounding sphere position
            mol1.boundingSphere.center.copy(mol1.position);

            for (let j = i + 1; j < this.allMolecules.length; j++) {
                const mol2 = this.allMolecules[j];
                mol2.boundingSphere.center.copy(mol2.position);

                // Check if bounding spheres intersect
                const distance = mol1.position.distanceTo(mol2.position);
                const minDistance =
                    mol1.boundingSphere.radius + mol2.boundingSphere.radius;

                if (distance < minDistance) {
                    this.handleCollision(mol1, mol2);
                }
            }
        }
    }

    handleCollision(mol1, mol2) {
        // Calculate normal vector between molecules
        const normal = new THREE.Vector3()
            .subVectors(mol2.position, mol1.position)
            .normalize();

        // Calculate relative velocity
        const relativeVelocity = new THREE.Vector3().subVectors(
            mol2.velocity,
            mol1.velocity,
        );

        // Calculate relative velocity along normal
        const velocityAlongNormal = relativeVelocity.dot(normal);

        // If molecules are moving apart, skip collision response
        if (velocityAlongNormal > 0) return;

        // Calculate impulse scalar
        const totalMass = mol1.mass + mol2.mass;
        const j = -(1 + this.collisionElasticity) * velocityAlongNormal;
        const impulse = j / totalMass;

        // Apply impulse to velocities
        const impulseVector = normal.multiplyScalar(impulse);

        mol1.velocity.sub(impulseVector.multiplyScalar(mol2.mass));
        mol2.velocity.add(impulseVector.multiplyScalar(mol1.mass));

        // Transfer some linear momentum to angular momentum
        const rotationImpulse = new THREE.Vector3(
            (Math.random() - 0.5) * this.rotationTransfer,
            (Math.random() - 0.5) * this.rotationTransfer,
            (Math.random() - 0.5) * this.rotationTransfer,
        );

        mol1.rotationVelocity.add(rotationImpulse);
        mol2.rotationVelocity.add(rotationImpulse.negate());

        // Apply velocity constraints
        [mol1, mol2].forEach((mol) => {
            mol.velocity.multiplyScalar(this.collisionDamping);
            const speed = mol.velocity.length();
            if (speed > this.maxVelocity) {
                mol.velocity.multiplyScalar(this.maxVelocity / speed);
            } else if (speed < this.minVelocity) {
                mol.velocity.multiplyScalar(this.minVelocity / speed);
            }
        });

        // Separate molecules to prevent sticking
        const overlap =
            mol1.boundingSphere.radius +
            mol2.boundingSphere.radius -
            mol1.position.distanceTo(mol2.position);
        if (overlap > 0) {
            const separationVector = normal.multiplyScalar(overlap * 0.5);
            mol1.position.sub(separationVector);
            mol2.position.add(separationVector);
        }
    }

    async createMolecules(structure = null) {
        const moleculeData = structure || this.getDefaultStructure();
        this.allMolecules = []; // Clear existing molecules array

        // Create main molecules
        for (let i = 0; i < this.numMolecules; i++) {
            const molecule = this.createMoleculeFromStructure(
                moleculeData,
                1,
                false,
            );
            this.molecules.push(molecule);
        }

        // Create water molecules
        const waterStructure = this.getWaterStructure();
        for (let i = 0; i < this.numWaterMolecules; i++) {
            const waterMolecule = this.createMoleculeFromStructure(
                waterStructure,
                this.waterScale,
                true,
            );
            this.waterMolecules.push(waterMolecule);
        }

        this.animate();
    }

    animate() {
        requestAnimationFrame(() => this.animate());

        this.updateParticles();
        this.checkCollisions();

        // Update molecules
        this.allMolecules.forEach((molecule) => {
            molecule.position.add(molecule.velocity);
            molecule.rotation.x += molecule.rotationVelocity.x;
            molecule.rotation.y += molecule.rotationVelocity.y;
            molecule.rotation.z += molecule.rotationVelocity.z;

            ["x", "y", "z"].forEach((axis) => {
                const limit =
                    axis === "z"
                        ? molecule.isWater
                            ? 15
                            : 20
                        : molecule.isWater
                          ? 100
                          : 128;
                if (Math.abs(molecule.position[axis]) > limit) {
                    molecule.velocity[axis] *= -0.95;
                    if (molecule.isWater) {
                        molecule.velocity[axis] +=
                            (Math.random() - 0.5) * 0.0005;
                    }
                }
            });
        });

        this.renderer.render(this.scene, this.camera);
    }

    onWindowResize() {
        const width = window.innerWidth;
        const height = window.innerHeight;

        this.camera.aspect = width / height;
        this.camera.updateProjectionMatrix();

        this.renderer.setSize(width, height);
        this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    }

    switchToConformer(modelNum) {
        if (this.conformers[modelNum]) {
            const structure = {
                atoms: this.conformers[modelNum].atoms,
                bonds: Array.from(this.conformers[modelNum].bonds),
            };
            // Clear existing molecules and water
            this.molecules.forEach((m) => this.scene.remove(m));
            this.waterMolecules.forEach((w) => this.scene.remove(w));
            this.molecules = [];
            this.waterMolecules = [];
            // Create new molecules with the selected conformer
            this.createMolecules(structure);
        }
    }
}
