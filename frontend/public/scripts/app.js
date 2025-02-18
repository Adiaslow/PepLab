// Update app.js
import { Router } from "./router.js";
import { routes } from "./routes.js";
import { MoleculeManager } from "./molecule-manager.js";
import { Menu } from "./menu.js";
import { initThemeToggle } from "./theme-toggle.js";

// app.js
document.addEventListener("DOMContentLoaded", async function () {
    try {
        const viewer = new MoleculeViewer(5);
        const moleculeManager = new MoleculeManager(viewer);
        await moleculeManager.initialize();

        // Wait for initial route template to be inserted
        const setupRouter = () => {
            const router = new Router(routes);
            window.router = router;
        };

        // Delay router setup slightly
        setTimeout(setupRouter, 200);

        new Menu();
        new FileHandler();
        initThemeToggle();
    } catch (error) {
        console.error("Initialization error:", error);
    }
});
