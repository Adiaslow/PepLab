document.addEventListener("DOMContentLoaded", function () {
    console.log("🚀 design.js loaded successfully!");

    // Get elements
    const openModal = document.getElementById("open-parameters");
    const closeModal = document.getElementById("close-modal");
    const submitBtn = document.getElementById("submit-parameters");
    const modal = document.getElementById("parameters-modal");
    const backButton = document.getElementById("back-button");
    const runButton = document.getElementById("run-button");

    // Debugging
    console.log("🔍 Checking for modal elements...");
    console.log("Open Parameters Button:", openModal);
    console.log("Close Button:", closeModal);
    console.log("Submit Button:", submitBtn);
    console.log("Modal:", modal);

    if (!openModal) console.error("🚨 Open Parameters button is missing!");
    if (!closeModal) console.error("🚨 Close button is missing!");
    if (!submitBtn) console.error("🚨 Submit button is missing!");
    if (!modal) console.error("🚨 Modal is missing!");

    // Open Modal
    if (openModal) {
        openModal.addEventListener("click", function () {
            console.log("✅ Open Parameters button clicked!");
            modal.style.display = "block";
        });
    }

    // Close Modal
    if (closeModal) {
        closeModal.addEventListener("click", function () {
            console.log("✅ Close button clicked!");
            modal.style.display = "none";
        });
    }

    // Handle Submit Button
    if (submitBtn) {
        submitBtn.addEventListener("click", async function () {
            console.log("✅ Submit button clicked!");

            const sequence = document.getElementById("sequence").value;
            const numPermutations = document.getElementById("num_permutations").value;

            if (!sequence || !numPermutations) {
                alert("Please enter both values.");
                return;
            }

            try {
                const urlParts = window.location.pathname.split("/");
                const strategyValue = urlParts[urlParts.length - 1];
                
                let parsedInput = sequence.includes(",") ? sequence.split(",") : sequence.split(/\s+/);
                parsedInput = parsedInput.map(s => s.trim()).filter(s => s.length > 0);

                const response = await fetch("/api/design/generate", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({
                        strategy: strategyValue,
                        input_data: parsedInput,
                        r: parseInt(numPermutations)
                    }),
                });

                const result = await response.json();
                console.log("Generated Sequences:", result);
                displayResults(result);
                modal.style.display = "none";
            } catch (error) {
                console.error("Error:", error);
                alert("Failed to generate peptides.");
            }
        });
    }

    // Handle Back Button
    if (backButton) {
        backButton.addEventListener("click", function () {
            history.back();
        });
    }

    // Handle Run Button
    if (runButton) {
        runButton.addEventListener("click", function () {
            console.log("🚀 Running design method:", window.location.pathname);
        });
    }
});

// Function to display generated sequences
function displayResults(data) {
    let outputDiv = document.getElementById("output");

    if (!outputDiv) {
        outputDiv = document.createElement("div");
        outputDiv.id = "output";
        outputDiv.className = "results-container";
        document.body.appendChild(outputDiv);
    }

    outputDiv.innerHTML = "<h3>Generated Sequences:</h3>";

    if (data.result && data.result.length > 0) {
        const list = document.createElement("ul");
        data.result.forEach(seq => {
            const item = document.createElement("li");
            item.textContent = seq.join(" ");
            list.appendChild(item);
        });
        outputDiv.appendChild(list);
        
        // Add Save to Database button
        const saveBtn = document.createElement("button");
        saveBtn.className = "btn-save";
        saveBtn.textContent = "Save to Database";
        saveBtn.onclick = function() {
            openSaveModal(data.result);
        };
        outputDiv.appendChild(saveBtn);
    } else {
        outputDiv.innerHTML += "<p>No sequences generated.</p>";
    }
}

// Open the save-name modal and handle the API call on confirm
function openSaveModal(sequences) {
    const saveModal = document.getElementById("save-name-modal");
    const nameInput = document.getElementById("library-name-input");
    const confirmBtn = document.getElementById("confirm-save");
    const cancelBtn = document.getElementById("cancel-save");

    if (!saveModal) return;

    nameInput.value = "";
    saveModal.style.display = "flex";
    setTimeout(() => nameInput.focus(), 50);

    // Remove any previous listeners to avoid stacking
    const newConfirm = confirmBtn.cloneNode(true);
    const newCancel = cancelBtn.cloneNode(true);
    confirmBtn.parentNode.replaceChild(newConfirm, confirmBtn);
    cancelBtn.parentNode.replaceChild(newCancel, cancelBtn);

    newCancel.addEventListener("click", () => {
        saveModal.style.display = "none";
    });

    newConfirm.addEventListener("click", async () => {
        const libraryName = nameInput.value.trim();
        if (!libraryName) {
            nameInput.focus();
            return;
        }

        newConfirm.disabled = true;
        newConfirm.textContent = "Saving…";

        try {
            const response = await fetch("/api/library/save", {
                method: "POST",
                headers: { "Content-Type": "application/json" },
                body: JSON.stringify({ name: libraryName, sequences })
            });
            const resData = await response.json();
            saveModal.style.display = "none";
            if (response.ok) {
                showToast(`Saved ${resData.peptide_count} peptides to "${libraryName}"`, "success");
            } else {
                showToast("Save failed: " + (resData.error || "Unknown error"), "error");
            }
        } catch (err) {
            saveModal.style.display = "none";
            showToast("Could not connect to API.", "error");
        } finally {
            newConfirm.disabled = false;
            newConfirm.textContent = "Save";
        }
    });

    // Also allow Enter key to confirm
    nameInput.addEventListener("keydown", function handler(e) {
        if (e.key === "Enter") { newConfirm.click(); nameInput.removeEventListener("keydown", handler); }
        if (e.key === "Escape") { newCancel.click(); nameInput.removeEventListener("keydown", handler); }
    });
}

// Lightweight toast notification — no alert() required
function showToast(message, type = "success") {
    const toast = document.createElement("div");
    toast.textContent = message;
    Object.assign(toast.style, {
        position: "fixed",
        bottom: "2rem",
        right: "2rem",
        padding: "0.85rem 1.4rem",
        borderRadius: "10px",
        background: type === "success" ? "#2ecc71" : "#e74c3c",
        color: "#fff",
        fontFamily: "'Cartograph CF', sans-serif",
        fontSize: "0.9rem",
        fontWeight: "500",
        zIndex: "9999",
        boxShadow: "0 4px 12px rgba(0,0,0,0.2)",
        opacity: "1",
        transition: "opacity 0.4s ease"
    });
    document.body.appendChild(toast);
    setTimeout(() => { toast.style.opacity = "0"; }, 2600);
    setTimeout(() => toast.remove(), 3100);
}

// Handle navigation between design methods
function handleDesignMethod(method) {
    console.log("Navigating to method:", method);
    window.location.href = `/design/${method}`;
}

// Handle combinatorial method selection
function handleCombinatoric(method) {
    console.log("Navigating to combinatoric method:", method);
    window.location.href = `/design/combinatoric/${method}`;
}