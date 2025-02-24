function saveThemeSettings(event) {
    event.preventDefault();
    
    const themeData = {
        theme: document.querySelector('#theme-select').value,
        accent_color: document.querySelector('#accent-select').value
    };

    fetch('/save_appearance', {
        method: 'POST',
        headers: {
            'Content-Type': 'application/json',
            'X-CSRFToken': document.querySelector('meta[name="csrf-token"]').content
        },
        body: JSON.stringify(themeData)
    })
    .then(response => {
        if (!response.ok) {
            throw new Error('Network response was not ok');
        }
        return response.json();
    })
    .then(data => {
        if (data.success) {
            window.location.reload();
        } else {
            console.error('Error:', data.error);
        }
    })
    .catch(error => {
        console.error('Error:', error);
    });
} 