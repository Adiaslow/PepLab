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
        saveBtn.className = "btn btn-primary mt-3";
        saveBtn.textContent = "Save to Database";
        saveBtn.onclick = async function() {
            const libraryName = prompt("Enter a name for this Library:", "Combinatorial Generation");
            if (!libraryName) return;
            
            saveBtn.disabled = true;
            saveBtn.textContent = "Saving...";
            try {
                const response = await fetch("/api/library/save", {
                    method: "POST",
                    headers: { "Content-Type": "application/json" },
                    body: JSON.stringify({
                        name: libraryName,
                        sequences: data.result
                    })
                });
                const resData = await response.json();
                if (response.ok) {
                    alert(`Successfully saved ${resData.peptide_count} combinations to Library: ${libraryName}`);
                } else {
                    alert("Failure: " + (resData.error || "Unknown Error"));
                }
            } catch (err) {
                alert("Failed to connect to API.");
            } finally {
                saveBtn.disabled = false;
                saveBtn.textContent = "Save to Database";
            }
        };
        outputDiv.appendChild(saveBtn);
    } else {
        outputDiv.innerHTML += "<p>No sequences generated.</p>";
    }
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