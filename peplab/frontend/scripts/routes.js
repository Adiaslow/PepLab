import { SettingsManager } from "./settings-manager.js";
import { SettingsToggle } from "./settings-toggle.js";
import { Dropdown } from "./dropdown.js";

export const routes = {
    "/home": {
        template: `
            <div class="route-content">
                <input
                    type="file"
                    id="fileInput"
                    accept=".csv"
                    style="display: none"
                />
                <h1>Home</h1>
                <button id="loadLibrary">Load Library</button>
                <button id="designLibrary">Design Library</button>
            </div>
        `,
        init: function () {
            document
                .getElementById("loadLibrary")
                .addEventListener("click", () => {
                    document.getElementById("fileInput").click();
                });

            // Add event listener for file input change
            document
                .getElementById("fileInput")
                .addEventListener("change", (event) => {
                    const file = event.target.files[0];
                    if (file) {
                        // Wait for file to be processed
                        const reader = new FileReader();
                        reader.onload = () => {
                            // After file is read, navigate to library page
                            window.location.hash = "/library";
                        };
                        reader.readAsText(file);
                    }
                });
            document
                .getElementById("designLibrary")
                .addEventListener("click", () => {
                    window.location.hash = "/design";
                });
        },
    },
    "/library": {
        template: `
            <div class="route-content">
                <h1>Library</h1>
                </div>
                `,
    },
    "/design": {
        template: `
                <div class="route-content">
                    <h1>Library Design</h1>
                    <!-- <div id="libraryOptions" class="library-options"></div> -->
                    <button id="combinatoric">Combinatoric</button>
                    <button id="generative">Generative</button>
                    <button id="genetic">Genetic</button>
                    <button id="mcmc">Markov Chain Monte Carlo</button>
                </div>
            `,
        init: function () {
            // Add event listeners for all buttons
            document
                .getElementById("combinatoric")
                .addEventListener("click", () => {
                    window.location.hash = "/combinatoric";
                });

            document.getElementById("mcmc").addEventListener("click", () => {
                window.location.hash = "/mcmc";
            });

            document.getElementById("genetic").addEventListener("click", () => {
                window.location.hash = "/genetic";
            });

            document
                .getElementById("generative")
                .addEventListener("click", () => {
                    window.location.hash = "/generative";
                });
        },
        /*
        init: function () {
            const dropdown = new Dropdown({
                items: [
                    "Combinatorial",
                    "Genetic Algorithm",
                    "Markov Chain Monte Carlo",
                    "Generative Neural Network",
                ],
                placeholder: "Design Method",
                onChange: (value, label) => {
                    console.log(`Selected design method: ${label}`);
                },
            });

            document
                .getElementById("libraryOptions")
                .appendChild(dropdown.element);
        },
         */
    },
    "/combinatoric": {
        template: `
            <div class="route-content">
                <h1>Combinatoric Design</h1>
                <button id="analyze">Cartesian Product</button>
                <button id="analyze">Combinative</button>
                <button id="analyze">Permutative</button>
                <button id="analyze">Group Theoretic</button>
            </div>
        `,
    },
    "/mcmc": {
        template: `
            <div class="route-content">
                <h1>MCMC Design</h1>
            </div>
        `,
    },
    "/genetic": {
        template: `
            <div class="route-content">
                <h1>Genetic Design</h1>
            </div>
        `,
    },
    "/generative": {
        template: `
            <div class="route-content">
                <h1>Generative Design</h1>
            </div>
        `,
    },
    "/analyze": {
        template: `
            <div class="route-content">
                <h1>Analyze</h1>
            </div>
        `,
    },
    "/visualize": {
        template: `
            <div class="route-content">
                <h1>Visualize</h1>
            </div>
        `,
    },
    "/settings": {
        template: `
            <div class="route-content">
                <h1>Settings</h1>
                <div id="settingsContainer" class="settings-container">
                    <!-- Toggles will be inserted here by SettingsManager -->
                </div>
            </div>
        `,
        init: function () {
            console.log("Settings route initialized");
            const settingsManager = new SettingsManager();
            settingsManager.initialize("settingsContainer");
        },
    },
};
