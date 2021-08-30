This repo contains the source code to generate the [mandy website](https://yhuang85.github.io/mandy/).

To make a similar website of your own, either start from scratch and follow [Sphinx tutorial](https://www.sphinx-doc.org/en/master/tutorial/index.html) or do the following

1. Clone the repo to your local machine by `git clone https://github.com/yhuang85/mandy_source.git`. Rename the folder **mandy_source** to whatever you like, but remember to carry your favorite name along in what follows.
2. Make sure you have installed a reasonably new [Python](https://www.python.org/) and [pipenv](https://github.com/pypa/pipenv).
3. Inside the folder **mandy_source**, run `pipenv sync` to install all the dependencies.
4. Modify the configurations such as project information, extensions etc in **docs/source/conf.py** (open it with your favorite text editor) to suit your needs. See more about Sphinx configurations [here](https://www.sphinx-doc.org/en/master/usage/configuration.html).
5. Run `make html` in the folder **mandy_source/docs** (inside the pipenv). This should create a folder **mandy_source/docs/build** which contains everything necessary to actually render a webpage.
6. Create a new Github Page repo and copy the folder **mandy_source/docs/build/html** to the repo. Rename the folder as **docs**. This is the content of [mandy](https://github.com/yhuang85/mandy.git).
7. In the new repo, go to *Settings* -> *Pages* and set the **Source** to be **main/docs**. Click **save** and you should see the url of your webpage. Of course, it's always a good idea to check out the [official documentation](https://docs.github.com/en/pages).

---
**Notes and Tips**

* The **build** folder is automatically created by `make`, and it's not included in the repo (see **.gitignore**).

* The `make html` command will overwrite the **build** folder but won't remove files in it which are no longer in use, e.g. outdated posts or images. Run `make clean` to remove all the content in **build**, or `make clean html` to clean up and rebuild in one go.

* Check out **conf.py** for examples of how Sphinx (>=4.0) supports [MathJax](http://docs.mathjax.org/en/latest) in an awesome way.
