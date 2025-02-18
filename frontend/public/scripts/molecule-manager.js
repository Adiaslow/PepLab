// molecule-manager.js
export class MoleculeManager {
    constructor(viewer) {
        this.viewer = viewer;
        this.setupEventListeners();
    }

    setupEventListeners() {
        document.addEventListener("pdbFileLoaded", async (event) => {
            await this.handleNewPDBFile(event.detail.content);
        });
    }

    async handleNewPDBFile(content) {
        this.viewer.molecules = [];
        await this.viewer.loadPDBFile(content);
        await this.createRandomMolecules(5);
    }

    async createRandomMolecules(numMolecules) {
        for (let i = 0; i < numMolecules; i++) {
            const randomConformerIndex = Math.floor(
                Math.random() * this.viewer.conformers.length,
            );
            const conformer = this.viewer.conformers[randomConformerIndex];
            if (conformer) {
                const structure = {
                    atoms: conformer.atoms,
                    bonds: Array.from(conformer.bonds),
                };
                const molecule =
                    this.viewer.createMoleculeFromStructure(structure);
                this.viewer.molecules.push(molecule);
            }
        }
    }

    async initialize() {
        await this.viewer.loadPDBFile("oxytocin_conformers.pdb");
        await this.createRandomMolecules(6);
        this.viewer.animate();
    }
}
