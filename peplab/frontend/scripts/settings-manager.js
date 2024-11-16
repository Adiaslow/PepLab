export class SettingsManager {
    constructor() {
        this.toggles = new Map();
        this.settingsContainer = null;
    }

    initialize(containerId) {
        this.settingsContainer = document.getElementById(containerId);
        if (!this.settingsContainer) {
            console.error(
                `Settings container with ID '${containerId}' not found`,
            );
            return;
        }

        // Initialize toggles
        this.createThemeToggle();
        this.createWaterToggle();
        this.createPeptideToggle();
        this.createParticleToggle();
    }

    createThemeToggle() {
        const themeToggle = new SettingsToggle({
            id: "themeToggle",
            label: "Dark Theme",
            storageKey: "theme",
            initialState: window.matchMedia("(prefers-color-scheme: dark)")
                .matches,
            onChange: (state) => {
                document.documentElement.setAttribute(
                    "data-theme",
                    state ? "dark" : "light",
                );
            },
        });
        this.addToggle("theme", themeToggle);
    }

    createWaterToggle() {
        const waterToggle = new SettingsToggle({
            id: "waterToggle",
            label: "Show Water",
            storageKey: "show-water",
            initialState: true,
            onChange: (state) => {
                // Add your water visibility logic here
            },
        });
        this.addToggle("water", waterToggle);
    }

    createPeptideToggle() {
        const peptideToggle = new SettingsToggle({
            id: "peptideToggle",
            label: "Show Peptides",
            storageKey: "show-peptides",
            initialState: true,
            onChange: (state) => {
                // Add your peptide visibility logic here
            },
        });
        this.addToggle("peptides", peptideToggle);
    }

    createParticleToggle() {
        const particleToggle = new SettingsToggle({
            id: "particleToggle",
            label: "Show Particles",
            storageKey: "show-particles",
            initialState: true,
            onChange: (state) => {
                // Add your particle visibility logic here
            },
        });
        this.addToggle("particles", particleToggle);
    }

    addToggle(key, toggle) {
        this.toggles.set(key, toggle);
        if (this.settingsContainer) {
            this.settingsContainer.appendChild(toggle.getElement());
        }
    }

    getToggle(key) {
        return this.toggles.get(key);
    }
}
