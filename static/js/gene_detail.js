// static/js/gene_detail.js

function filterAnnotations() {
    const aspect = document.getElementById('aspectFilter').value;
    const evidenceType = document.getElementById('evidenceFilter').value;
    const rows = document.querySelectorAll('.annotations-table tbody tr');
    let visibleCount = 0;

    // Classification of Evidence Codes
    const experimental = ['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'];
    const computational = ['ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'RCA'];
    const electronic = ['IEA'];

    rows.forEach(row => {
        const rowAspect = row.getAttribute('data-aspect');
        const rowEvidence = row.getAttribute('data-evidence');

        let aspectMatch = (aspect === 'all' || rowAspect === aspect);
        let evidenceMatch = false;

        if (evidenceType === 'all') {
            evidenceMatch = true;
        } else if (evidenceType === 'EXP') {
            evidenceMatch = experimental.includes(rowEvidence);
        } else if (evidenceType === 'COMP') {
            evidenceMatch = computational.includes(rowEvidence);
        } else if (evidenceType === 'IEA') {
            evidenceMatch = electronic.includes(rowEvidence);
        } else if (evidenceType === 'OTHER') {
            evidenceMatch = !experimental.includes(rowEvidence) &&
                           !computational.includes(rowEvidence) &&
                           !electronic.includes(rowEvidence);
        }

        if (aspectMatch && evidenceMatch) {
            row.style.display = '';
            visibleCount++;
        } else {
            row.style.display = 'none';
        }
    });

    document.getElementById('filteredCountBadge').innerText = `Showing ${visibleCount} terms`;
}

function resetFilters() {
    document.getElementById('aspectFilter').value = 'all';
    document.getElementById('evidenceFilter').value = 'all';
    filterAnnotations();
}
