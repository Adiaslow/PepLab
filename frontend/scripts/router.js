export class Router {
    constructor(routes) {
        this.routes = routes;
        this.app = document.getElementById("app");
        this.container = this.app.querySelector(".container");
        this.currentContent = null;

        if (!this.app || !this.container) {
            console.error("Required elements not found");
            return;
        }

        window.addEventListener("hashchange", () => this.handleRoute());
        this.handleRoute();
    }

    async handleRoute() {
        const hash = window.location.hash || "#/";
        const path = hash.slice(1);
        const route = this.routes[path] || this.routes["/"];

        if (!route) {
            console.error(`No route found for path: ${path}`);
            return;
        }

        // Create a wrapper for the new content
        const contentWrapper = document.createElement("div");
        contentWrapper.className = "route-content";
        contentWrapper.innerHTML = route.template;

        // Create a measurement div
        const measureDiv = contentWrapper.cloneNode(true);
        measureDiv.style.position = "absolute";
        measureDiv.style.visibility = "hidden";
        measureDiv.style.width = "auto";
        measureDiv.style.height = "auto";
        document.body.appendChild(measureDiv);

        // Measure the content
        const newWidth = measureDiv.offsetWidth;
        const newHeight = measureDiv.offsetHeight;

        // Clean up measurement div
        document.body.removeChild(measureDiv);

        // Start transition
        this.container.classList.add("morphing");

        // If there's current content, fade it out
        if (this.currentContent) {
            this.currentContent.classList.add("fade-out");
            await new Promise((resolve) => setTimeout(resolve, 300)); // Match transition duration
        }

        // Clear container and add new content
        this.container.innerHTML = "";
        this.container.appendChild(contentWrapper);
        this.currentContent = contentWrapper;

        // Set new dimensions with a small buffer
        this.container.style.width = `${newWidth + 60}px`; // Add padding
        this.container.style.height = `${newHeight + 60}px`; // Add padding

        // Initialize route if needed
        if (route.init) {
            route.init();
        }

        // Remove transition class after animation
        setTimeout(() => {
            this.container.classList.remove("morphing");
        }, 300);
    }
}
