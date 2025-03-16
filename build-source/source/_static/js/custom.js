document.addEventListener('DOMContentLoaded', function () {
    const links = document.querySelectorAll('.reference.internal[href*="#equation-"]');
    links.forEach(function(link) {
        const tooltip = document.createElement('div');
        tooltip.className = 'eq-tooltip';
        link.classList.add('equation-link');

        const eqUrl = link.getAttribute('href');
        const eqId = eqUrl.split('#')[1];
        fetch(eqUrl)
            .then(response => response.text())
            .then(html => {
                const parser = new DOMParser();
                const doc = parser.parseFromString(html, 'text/html');
                const eqElem = doc.getElementById(eqId);
                if (eqElem) {
                    // Set eqno invisible if there is any
                    let eqnoElem = eqElem.querySelector('.eqno');
                    if (eqnoElem) {
                        eqnoElem.style.display = 'none';
                    }

                    tooltip.innerHTML = eqElem.innerHTML;

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
