(() => {
    const apiUrl = "http://127.0.0.1:5000/generate"; // Backend API URL

    // Function to send data to the backend
    async function generateComposition(strategy, inputData, r = null) {
        const payload = {
            strategy: strategy,
            input_data: inputData,
            r: r
        };

        console.log("Sending payload:", payload); // Debugging output

        try {
            showLoading(true); // Show loading indicator while processing

            const response = await fetch(apiUrl, {
                method: "POST",
                headers: {
                    "Content-Type": "application/json"
                },
                body: JSON.stringify(payload)
            });

            const result = await response.json();
            showLoading(false); // Hide loading indicator when done

            if (response.ok) {
                console.log("Response received:", result); // Debugging output
                displayResults(result.result); // Display results on the webpage
            } else {
                console.error("Error from server:", result.error);
                displayError(result.error || "An unknown error occurred.");
            }
        } catch (error) {
            showLoading(false); // Ensure the loading indicator is hidden on error
            console.error("Network error:", error);
            displayError("Failed to connect to the server. Please try again later.");
        }
    }

    // Function to display results on the webpage
    function displayResults(results) {
        const outputDiv = document.getElementById("output");
        if (!outputDiv) {
            console.error("Output div not found"); // Debugging message
            return;
        }
        outputDiv.innerHTML = results.map(res => `<div>${JSON.stringify(res)}</div>`).join("");
    }

    // Function to display an error message
    function displayError(errorMessage) {
        const outputDiv = document.getElementById("output");
        if (!outputDiv) {
            console.error("Output div not found"); // Debugging message
            return;
        }
        outputDiv.innerHTML = `<div class="error">${errorMessage}</div>`;
    }

    // Show or hide a loading indicator
    function showLoading(isLoading) {
        const loadingIndicator = document.getElementById("loading");
        if (loadingIndicator) {
            loadingIndicator.style.display = isLoading ? "block" : "none";
        }
    }
    function displayResults(results) {
        const outputDiv = document.getElementById("output");
        if (outputDiv) {
            outputDiv.innerHTML = results.map(result => `<div>${JSON.stringify(result)}</div>`).join("");
        }
    }
    

    // Add event listeners to buttons
    document.getElementById("combinativeButton").addEventListener("click", () => {
        const inputData = [1, 2, 3]; // Input data for combinative composition
        generateComposition("combinative", inputData, 2);
    });

    document.getElementById("cyclicButton").addEventListener("click", () => {
        const inputData = ["A", "B", "C"]; // Input data for cyclic composition
        generateComposition("cyclic", inputData);
    });

    // Add other buttons here as needed
})();
