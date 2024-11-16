// menu.js
export class Menu {
    constructor() {
        this.hamburger = document.querySelector(".hamburger");
        this.nav = document.querySelector("nav");
        this.isOpen = false;

        this.init();
    }

    init() {
        // Add the navigation links to the hamburger container
        const navContent = `
            <div class="menu-content">
                <div class="line"></div>
                <div class="line"></div>
                <div class="line"></div>
                <nav class="nav-links">
                    <a href="#/home" data-page="home">Home</a>
                    <a href="#/library" data-page="about">Library</a>
                    <a href="#/design" data-page="about">Design</a>
                    <a href="#/analyze" data-page="about">Analyze</a>
                    <a href="#/visualize" data-page="about">Visualize</a>
                    <a href="#/settings" data-page="settings">Settings</a>
                    <a href="https://github.com/Adiaslow/PepLab" data-page="github" target="_blank">About</a>
                </nav>
            </div>
        `;

        this.hamburger.innerHTML = navContent;

        this.hamburger.addEventListener("click", () => this.toggleMenu());

        // Close menu when clicking outside
        document.addEventListener("click", (e) => {
            if (!this.hamburger.contains(e.target)) {
                this.closeMenu();
            }
        });
    }

    toggleMenu() {
        this.isOpen = !this.isOpen;
        this.hamburger.classList.toggle("active");
    }

    closeMenu() {
        this.isOpen = false;
        this.hamburger.classList.remove("active");
    }
}
