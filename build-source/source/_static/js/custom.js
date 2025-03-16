document.addEventListener('DOMContentLoaded', function () {
    const links = document.querySelectorAll('.reference.internal[href*="#equation-"]');
    links.forEach(function(link) {
        // Create a div to be inserted under the equation link
        const tooltip = document.createElement('div');
        tooltip.className = 'eq-tooltip';
        link.classList.add('equation-link');

        const eqnUrl = link.getAttribute('href');
        const eqnId = eqnUrl.split('#')[1];
        fetch(eqnUrl)
            .then(response => response.text())
            .then(html => {
                const parser = new DOMParser();
                const doc = parser.parseFromString(html, 'text/html');
                const eqnElem = doc.getElementById(eqnId);
                if (eqnElem) {
                    // Set `eqno` invisible if there is any
                    const eqnoElem = eqnElem.querySelector('.eqno');
                    if (eqnoElem) {
                        eqnoElem.style.display = 'none';
                    }

                    tooltip.innerHTML = eqnElem.innerHTML;

                    if (typeof MathJax !== 'undefined') {
                        MathJax.typeset([tooltip]);
                    }

                    link.appendChild(tooltip);
                }
            })
            .catch(error => {
                console.log('Error fetching equation: ', error);
            });
    });
});
