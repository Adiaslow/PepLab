// settings-toggle.js
export class SettingsToggle {
    constructor(options) {
        this.options = {
            id: "",
            label: "",
            initialState: false,
            storageKey: "",
            onChange: () => {},
            ...options,
        };

        this.element = this.createToggle();
        this.init();
    }

    createToggle() {
        const wrapper = document.createElement("div");
        wrapper.className = "settings-toggle-wrapper";

        const label = document.createElement("label");
        label.className = "settings-label";
        label.textContent = this.options.label;

        const toggle = document.createElement("div");
        toggle.className = "settings-toggle";
        toggle.id = this.options.id;
        toggle.setAttribute("role", "switch");
        toggle.setAttribute("tabindex", "0");
        toggle.setAttribute("aria-checked", "false");

        wrapper.appendChild(label);
        wrapper.appendChild(toggle);

        return wrapper;
    }

    init() {
        const toggle = this.element.querySelector(".settings-toggle");

        // Get saved state
        const savedState = localStorage.getItem(this.options.storageKey);
        const initialState =
            savedState !== null
                ? savedState === "true"
                : this.options.initialState;

        this.updateState(initialState);

        // Add event listeners
        toggle.addEventListener("click", () => this.toggleState());
        toggle.addEventListener("keydown", (e) => {
            if (e.key === "Enter" || e.key === " ") {
                e.preventDefault();
                this.toggleState();
            }
        });
    }

    updateState(state) {
        const toggle = this.element.querySelector(".settings-toggle");
        toggle.classList.toggle("active", state);
        toggle.setAttribute("aria-checked", state);
        localStorage.setItem(this.options.storageKey, state);
        this.options.onChange(state);
    }

    toggleState() {
        const toggle = this.element.querySelector(".settings-toggle");
        const currentState = toggle.classList.contains("active");
        this.updateState(!currentState);
    }

    getElement() {
        return this.element;
    }
}
