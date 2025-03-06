class BuildingBlockExplorer {
    constructor(containerId) {
        this.container = document.getElementById(containerId);
        this.buildingBlocks = [];
        this.initializeUI();
    }

    initializeUI() {
        // Create explorer container
        const explorerContainer = document.createElement('div');
        explorerContainer.className = 'building-block-explorer';
        
        // Create load button
        const loadButton = document.createElement('button');
        loadButton.textContent = 'Load Building Blocks';
        loadButton.className = 'button primary';
        loadButton.addEventListener('click', () => this.loadBuildingBlocks());
        
        // Create grid container for building blocks
        this.gridContainer = document.createElement('div');
        this.gridContainer.className = 'building-block-grid';
        
        // Append elements
        explorerContainer.appendChild(loadButton);
        explorerContainer.appendChild(this.gridContainer);
        this.container.appendChild(explorerContainer);
    }

    async loadBuildingBlocks() {
        try {
            const response = await fetch('/api/building-blocks');
            const data = await response.json();
            
            if (data.error) {
                throw new Error(data.error);
            }
            
            this.buildingBlocks = data.building_blocks;
            this.displayBuildingBlocks();
        } catch (error) {
            console.error('Error loading building blocks:', error);
            // Show error message to user
            alert('Error loading building blocks: ' + error.message);
        }
    }

    displayBuildingBlocks() {
        this.gridContainer.innerHTML = '';
        
        this.buildingBlocks.forEach(block => {
            const card = document.createElement('div');
            card.className = 'building-block-card';
            
            card.innerHTML = `
                <div class="molecule-image">
                    <img src="${block.image}" alt="${block.name}">
                </div>
                <div class="molecule-info">
                    <h3>${block.name}</h3>
                    <p>Aliases: ${block.alt_name1} (${block.alt_name2})</p>
                    <p>SMILES: ${block.smiles}</p>
                    <p>Position: ${block.position}</p>
                </div>
            `;
            
            this.gridContainer.appendChild(card);
        });
    }
} 