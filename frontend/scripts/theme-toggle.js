export function initThemeToggle() {
    // Wait for DOM to be ready
    const themeToggle = document.getElementById("themeToggle");
    if (!themeToggle) {
        console.error("Theme toggle element not found");
        return;
    }

    // Get saved theme or use system preference as fallback
    const savedTheme =
        localStorage.getItem("theme") ||
        (window.matchMedia("(prefers-color-scheme: dark)").matches
            ? "dark"
            : "light");

    // Set initial state
    updateThemeState(savedTheme);

    function updateThemeState(theme) {
        // Update DOM
        document.documentElement.setAttribute("data-theme", theme);
        // Update toggle state (active = dark theme now)
        themeToggle.classList.toggle("active", theme === "dark");
        themeToggle.setAttribute("aria-checked", theme === "dark");
        // Save preference
        localStorage.setItem("theme", theme);
    }

    function toggleTheme() {
        const currentTheme =
            document.documentElement.getAttribute("data-theme");
        const newTheme = currentTheme === "light" ? "dark" : "light";
        updateThemeState(newTheme);
    }

    // Add event listeners
    themeToggle.addEventListener("click", toggleTheme);
    themeToggle.addEventListener("keydown", (e) => {
        if (e.key === "Enter" || e.key === " ") {
            e.preventDefault();
            toggleTheme();
        }
    });
}
