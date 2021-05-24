
.. _editing_this_user_guide:

Editing this User Guide
=======================

This user guide is generated with `Sphinx <https://www.sphinx-doc.org/>`_. 
Sphinx is an open-source Python project designed to make writing software documentation easier. 
The documentation is written in a reStructuredText (it's similar to markdown), which Sphinx extends for software documentation.
The source for the documentation is the :file:`docs/source` directory in top-level of the source code.

Quick start
-----------

To build this user guide on your local machine, you need to install Sphinx. Sphinx is a Python 3 package and
it is available via :program:`pip`. This user guide uses the Read The Docs theme, so you will also need to 
install :literal:`sphinx-rtd-theme`. It also uses the `sphinxcontrib-bibtex <https://pypi.org/project/sphinxcontrib-bibtex/>`_
and `recommonmark <https://recommonmark.readthedocs.io/>`_ extensions, which you'll need to install.

.. code-block:: console

   $ pip install sphinx sphinx-rtd-theme sphinxcontrib-bibtex recommonmark


To build this user guide locally, navigate to the :file:`docs/` directory and make the :literal:`html` target.

.. code-block:: console

   gcuser:~$ cd gcpy/docs
   gcuser:~/gcpy/docs$ make html


This will build the user guide in :file:`docs/build/html`, and you can open :file:`index.html` in your 
web-browser. The source files for the user guide are found in :file:`docs/source`.  

.. note::

   You can clean the documentation with :code:`make clean`.


Learning reST
-------------

Writing reST can be tricky at first. Whitespace matters, and some directives
can be easily miswritten. Two important things you should know right away are:

* Indents are 3-spaces
* "Things" are separated by 1 blank line. For example, a list or code-block following a paragraph should be separated from the paragraph by 1 blank line.

You should keep these in mind when you're first getting started. Dedicating an hour to learning reST
will save you time in the long-run. Below are some good resources for learning reST.

* `reStructuredText primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_: (single best resource; however, it's better read than skimmed)
* Official `reStructuredText reference <https://docutils.sourceforge.io/docs/user/rst/quickref.html>`_ (there is *a lot* of information here)
* `Presentation by Eric Holscher <https://www.youtube.com/watch?v=eWNiwMwMcr4>`_ (co-founder of Read The Docs) at DjangoCon US 2015 (the entire presentation is good, but reST is described from 9:03 to 21:04)
* `YouTube tutorial by Audrey Tavares's <https://www.youtube.com/watch?v=DSIuLnoKLd8>`_

A good starting point would be Eric Holscher's presentations followed by the reStructuredText primer.

Style guidelines
----------------

.. important::

   This user guide is written in semantic markup. This is important so that the user guide remains
   maintainable. Before contributing to this documentation, please review our style guidelines
   (below). When editing the source, please refrain from using elements with the wrong semantic
   meaning for aesthetic reasons. Aesthetic issues can be addressed by changes to the theme.

For **titles and headers**:

* Section headers should be underlined by :literal:`#` characters
* Subsection headers should be underlined by :literal:`-` characters
* Subsubsection headers should be underlined by :literal:`^` characters
* Subsubsubsection headers should be avoided, but if necessary, they should be underlined by :literal:`"` characters

**File paths** (including directories) occuring in the text should use the :literal:`:file:` role.

**Program names** (e.g. :program:`cmake`) occuring in the text should use the :literal:`:program:` role.

**OS-level commands** (e.g. :command:`rm`) occuring in the text should use the :literal:`:command:` role.

**Environment variables** occuring in the text should use the :literal:`:envvar:` role.

**Inline code** or code variables occuring in the text should use the :literal:`:code:` role.

**Code snippets** should use :literal:`.. code-block:: <language>` directive like so

.. code-block:: none

   .. code-block:: python

      import gcpy
      print("hello world")

The language can be "none" to omit syntax highlighting. 

For command line instructions, the "console" language should be used. The :literal:`$` should be used
to denote the console's prompt. If the current working directory is relevant to the instructions,
a prompt like :literal:`gcuser:~/path1/path2$` should be used.

**Inline literals** (e.g. the :literal:`$` above) should use the :literal:`:literal:` role.