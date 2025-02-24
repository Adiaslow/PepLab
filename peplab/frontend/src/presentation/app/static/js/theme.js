function setTheme(theme) {
    if (theme === 'system') {
        // Check system preference
        if (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches) {
            document.documentElement.setAttribute('data-theme', 'dark');
        } else {
            document.documentElement.setAttribute('data-theme', 'light');
        }
    } else {
        document.documentElement.setAttribute('data-theme', theme);
    }
    
    // Save theme preference
    localStorage.setItem('theme', theme);
}

function initTheme() {
    // Get saved theme or default to system
    const savedTheme = localStorage.getItem('theme') || 'system';
    setTheme(savedTheme);
    
    // Update theme selector if it exists
    const themeSelect = document.getElementById('theme');
    if (themeSelect) {
        themeSelect.value = savedTheme;
    }
    
    // Listen for system theme changes
    if (window.matchMedia) {
        window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', e => {
            if (localStorage.getItem('theme') === 'system') {
                setTheme('system');
            }
        });
    }
}

// Initialize theme on page load
document.addEventListener('DOMContentLoaded', function() {
    const themeOptions = document.querySelectorAll('.theme-option');
    const saveButton = document.querySelector('.settings-save');
    let selectedTheme = localStorage.getItem('theme') || 'light';

    // Set initial theme
    document.documentElement.setAttribute('data-theme', selectedTheme);

    // Update selected theme option
    themeOptions.forEach(option => {
        if (option.dataset.theme === selectedTheme) {
            option.classList.add('selected');
        }

        option.addEventListener('click', () => {
            // Update selection UI
            themeOptions.forEach(opt => opt.classList.remove('selected'));
            option.classList.add('selected');
            
            // Update theme immediately
            selectedTheme = option.dataset.theme;
            document.documentElement.setAttribute('data-theme', selectedTheme);
        });
    });

    // Save theme preference
    saveButton?.addEventListener('click', () => {
        localStorage.setItem('theme', selectedTheme);
    });
}); 