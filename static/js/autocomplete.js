// static/js/autocomplete.js

function setupTermAutocomplete(queryInputId, idInputId, resultsContainerId) {
    const queryInput = document.getElementById(queryInputId);
    const idInput = document.getElementById(idInputId);
    const resultsContainer = document.getElementById(resultsContainerId);

    // If any element is missing on the page, silently do nothing — makes the
    // function safe to call from any template without conditionals.
    if (!queryInput || !idInput || !resultsContainer) return;

    let debounceTimer;

    queryInput.addEventListener('input', () => {
        clearTimeout(debounceTimer);
        const query = queryInput.value;
        idInput.value = ''; // Clear the hidden ID when the user types

        if (query.length < 3) {
            resultsContainer.innerHTML = '';
            return;
        }

        // Debounce to avoid sending too many requests while typing
        debounceTimer = setTimeout(() => {
            fetch(`/api/search-terms?q=${encodeURIComponent(query)}`)
                .then(response => response.json())
                .then(data => {
                    resultsContainer.innerHTML = '';
                    if (data.results.length > 0) {
                        data.results.forEach(term => {
                            const item = document.createElement('a');
                            item.href = '#';

                            item.innerHTML = `<strong>${term.id}</strong> - ${term.name} <span class="muted-small">${term.namespace}</span>`;
                            item.addEventListener('click', (e) => {
                                e.preventDefault();
                                queryInput.value = `${term.id} - ${term.name}`;
                                idInput.value = term.id;
                                resultsContainer.innerHTML = '';
                            });
                            resultsContainer.appendChild(item);
                        });
                    } else {
                        resultsContainer.innerHTML = '<div class="no-results">No matches found</div>';
                    }
                });
        }, 250); // Wait for 250ms after typing stops
    });

    // Hide results when clicking outside
    document.addEventListener('click', (e) => {
        if (!resultsContainer.contains(e.target) && e.target !== queryInput) {
            resultsContainer.innerHTML = '';
        }
    });
}


function setupGeneAutocomplete(queryInputId, resultsContainerId) {
    const queryInput = document.getElementById(queryInputId);
    const resultsContainer = document.getElementById(resultsContainerId);

    if (!queryInput || !resultsContainer) return;

    let debounceTimer;

    queryInput.addEventListener('input', () => {
        clearTimeout(debounceTimer);
        const query = queryInput.value;

        if (query.length < 2) {
            resultsContainer.innerHTML = '';
            return;
        }

        debounceTimer = setTimeout(() => {
            fetch(`/api/search-genes?q=${encodeURIComponent(query)}`)
                .then(response => response.json())
                .then(data => {
                    resultsContainer.innerHTML = '';
                    if (data.results.length > 0) {
                        data.results.forEach(gene => {
                            const item = document.createElement('a');
                            item.href = '#';
                            item.innerHTML = `<strong>${gene.symbol}</strong> - ${gene.name}`;
                            item.addEventListener('click', (e) => {
                                e.preventDefault();
                                queryInput.value = gene.symbol;
                                resultsContainer.innerHTML = '';
                            });
                            resultsContainer.appendChild(item);
                        });
                    } else {
                        resultsContainer.innerHTML = '<div class="no-results">No matches found</div>';
                    }
                });
        }, 250);
    });

    document.addEventListener('click', (e) => {
        if (!resultsContainer.contains(e.target) && e.target !== queryInput) {
            resultsContainer.innerHTML = '';
        }
    });
}
