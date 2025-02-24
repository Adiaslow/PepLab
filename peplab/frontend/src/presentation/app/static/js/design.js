// Handle initial design method selection
function handleDesignMethod(method) {
    window.location.href = `/design/${method}`;
}

// Handle combinatoric method selection
function handleCombinatoric(method) {
    window.location.href = `/design/combinatoric/${method}`;
}

// Handle method parameters and execution
function handleRun() {
    // Get parameters from form
    const parameters = {
        // Add parameter collection logic here
    };

    // TODO: Implement API call to run the design method
    console.log('Running design with parameters:', parameters);
}

// Handle building blocks file upload
function handleBlocksUpload(input) {
    if (input.files.length > 0) {
        const formData = new FormData();
        formData.append('blocks', input.files[0]);
        
        fetch('/design/upload-blocks', {
            method: 'POST',
            body: formData
        }).then(response => {
            window.location.reload();
        }).catch(error => {
            console.error('Upload failed:', error);
        });
    }
} 