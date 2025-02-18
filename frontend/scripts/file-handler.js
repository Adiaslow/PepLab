// file-handler.js
class FileHandler {
    constructor() {
        this.fileInput = document.getElementById("fileInput");
        this.loadButton = document.getElementById("loadLibrary");
        this.setupEventListeners();
    }

    setupEventListeners() {
        this.loadButton.addEventListener("click", () =>
            this.openFileExplorer(),
        );
        this.fileInput.addEventListener("change", (event) =>
            this.handleFileSelect(event),
        );
    }

    openFileExplorer() {
        this.fileInput.click();
    }

    async handleFileSelect(event) {
        const file = event.target.files[0];
        if (file) {
            try {
                const content = await file.text();
                // Parse CSV content
                const rows = content.split("\n");
                const headers = rows[0].split(",");
                const data = rows.slice(1).map((row) => {
                    const values = row.split(",");
                    const rowData = {};
                    headers.forEach((header, index) => {
                        rowData[header.trim()] = values[index]?.trim();
                    });
                    return rowData;
                });

                console.log("Loaded CSV data:", data); // For debugging

                // Dispatch event with parsed data
                const fileLoadEvent = new CustomEvent("csvLoaded", {
                    detail: {
                        headers: headers,
                        data: data,
                    },
                });
                document.dispatchEvent(fileLoadEvent);
            } catch (error) {
                console.error("Error loading file:", error);
                alert("Error loading CSV file. Please try again.");
            }
        }
    }
}
