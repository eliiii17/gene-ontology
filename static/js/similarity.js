// static/js/similarity.js

// Set up autocompletion for term inputs
setupTermAutocomplete('term_a_query', 'term_a_id', 'results_a');
setupTermAutocomplete('term_b_query', 'term_b_id', 'results_b');

// Set up autocompletion for gene inputs
setupGeneAutocomplete('gene_a_query', 'results_gene_a');
setupGeneAutocomplete('gene_b_query', 'results_gene_b');

// --- Dynamic Form Logic ---
document.addEventListener('DOMContentLoaded', function () {
    const modeRadios = document.querySelectorAll('input[name="mode"]');
    const termInputsContainer = document.getElementById('term-inputs');
    const geneInputsContainer = document.getElementById('gene-inputs');
    // Check if elements exist to avoid errors on pages that might include this script partially
    if (!termInputsContainer || !geneInputsContainer) return;

    const termQueryInputs = termInputsContainer.querySelectorAll('input[type="text"]');
    const geneQueryInputs = geneInputsContainer.querySelectorAll('input[type="text"]');

    function setRequired(inputs, isRequired) {
        inputs.forEach(input => {
            if (isRequired) {
                input.setAttribute('required', '');
            } else {
                input.removeAttribute('required');
            }
        });
    }

    function updateFormVisibility() {
        const selectedMode = document.querySelector('input[name="mode"]:checked').value;

        termInputsContainer.classList.add('hidden');
        geneInputsContainer.classList.add('hidden');
        document.getElementById('matrix-inputs').classList.add('hidden');

        setRequired(termQueryInputs, false);
        setRequired(geneQueryInputs, false);
        const geneListInput = document.querySelector('input[name="gene_list_query"]');
        if (geneListInput) geneListInput.removeAttribute('required');

        if (selectedMode === 'term') {
            termInputsContainer.classList.remove('hidden');
            setRequired(termQueryInputs, true);
        } else if (selectedMode === 'gene') {
            geneInputsContainer.classList.remove('hidden');
            setRequired(geneQueryInputs, true);
        } else if (selectedMode === 'matrix') {
            document.getElementById('matrix-inputs').classList.remove('hidden');
            if (geneListInput) geneListInput.setAttribute('required', '');
        }
    }

    // Add event listeners to radio buttons
    modeRadios.forEach(radio => radio.addEventListener('change', updateFormVisibility));

    // Set initial state on page load
    updateFormVisibility();
});

function loadTopGenes() {
    const btn = document.querySelector('button[onclick="loadTopGenes()"]');
    const originalText = btn.innerText;
    btn.innerText = "Loading...";
    btn.disabled = true;

    fetch('/api/top-genes')
        .then(response => response.json())
        .then(data => {
            const input = document.getElementById('gene_list_query');
            input.value = data.genes;
            btn.innerText = originalText;
            btn.disabled = false;
        })
        .catch(err => {
            console.error(err);
            btn.innerText = "Error";
            btn.disabled = false;
        });
}

// --- Form Submission Fallback ---
const form = document.querySelector('.similarity-form');
if (form) {
    form.addEventListener('submit', function (e) {
        const selectedModeInput = document.querySelector('input[name="mode"]:checked');
        if (!selectedModeInput) return;
        
        const selectedMode = selectedModeInput.value;
        if (selectedMode === 'term') {
            const termAIdInput = document.getElementById('term_a_id');
            const termBIdInput = document.getElementById('term_b_id');
            if (!termAIdInput.value) {
                termAIdInput.value = document.getElementById('term_a_query').value;
            }
            if (!termBIdInput.value) {
                termBIdInput.value = document.getElementById('term_b_query').value;
            }
        }
    });
}
