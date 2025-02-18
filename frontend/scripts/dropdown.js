// dropdown.js
export class Dropdown {
    constructor(options = {}) {
        this.options = {
            items: options.items || [],
            placeholder: options.placeholder || "Select an option",
            onChange: options.onChange || (() => {}),
            // width: options.width || "280px",
        };

        this.isOpen = false;
        this.selectedItem = null;
        this.element = this.createElement();
        this.setupEventListeners();
    }

    createElement() {
        const wrapper = document.createElement("div");
        wrapper.className = "dropdown";
        wrapper.innerHTML = `
            <button class="dropdown-trigger">
                <span class="dropdown-selected">${this.options.placeholder}</span>
            </button>
            <div class="dropdown-content">
                ${this.options.items
                    .map(
                        (item, index) => `
                    <a class="dropdown-item" data-value="${item.value || item}">
                        ${item.label || item}
                    </a>
                `,
                    )
                    .join("")}
            </div>
        `;
        return wrapper;
    }

    setupEventListeners() {
        const trigger = this.element.querySelector(".dropdown-trigger");
        const content = this.element.querySelector(".dropdown-content");
        const items = this.element.querySelectorAll(".dropdown-item");

        trigger.addEventListener("click", (e) => {
            e.stopPropagation();
            this.isOpen = !this.isOpen;
            trigger.classList.toggle("active");
            content.classList.toggle("active");
        });

        items.forEach((item) => {
            item.addEventListener("click", (e) => {
                const value = e.target.dataset.value;
                const label = e.target.textContent.trim();
                this.selectItem(value, label);
                this.closeDropdown();
            });
        });

        document.addEventListener("click", (e) => {
            if (!this.element.contains(e.target)) {
                this.closeDropdown();
            }
        });
    }

    selectItem(value, label) {
        this.selectedItem = { value, label };
        this.element.querySelector(".dropdown-selected").textContent = label;
        this.options.onChange(value, label);

        // Update selected state
        this.element.querySelectorAll(".dropdown-item").forEach((item) => {
            item.classList.toggle("selected", item.dataset.value === value);
        });
    }

    closeDropdown() {
        this.isOpen = false;
        this.element
            .querySelector(".dropdown-trigger")
            .classList.remove("active");
        this.element
            .querySelector(".dropdown-content")
            .classList.remove("active");
    }

    getValue() {
        return this.selectedItem?.value;
    }

    setValue(value) {
        const item = this.element.querySelector(
            `.dropdown-item[data-value="${value}"]`,
        );
        if (item) {
            this.selectItem(value, item.textContent.trim());
        }
    }
}
