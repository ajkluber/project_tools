Things To Do and Keep in Mind
=============================

Things To Keep in Mind
----------------------

Our overall coding goals (in roughly decreasing importance):

1. Write code that is easy to understand and extend: style, design, and documentation.

  * Please force yourself to adhere to `PEP 8 style guide. <http://legacy.python.org/dev/peps/pep-0008>`_

    Standardizing our coding style promotes readability. 

    20 useful adages: `Zen of Python <http://legacy.python.org/dev/peps/pep-0020/>`_

  * Program design and development is nontrivial, ideally we would follow the 
    new school philosophy of `"agile development" <http://en.wikipedia.org/wiki/Agile_software_development>`_
    where development process involves many small iterations. 

2. Write code that works: unit testing.

  * Code is broken until tested.

3. Eventually write code that is efficient with time and memory, but only after #2.


Things To Do
------------

Things to be improved, tested, or implemented for the first time.

General Objectives
^^^^^^^^^^^^^^^^^^

1. Generally, we want to reduce our dependence on outside packages (e.g.
   MODELLER, Gromacs commands) and calling command line commands when these
   tasks can be done in pure Python. For example, some commands may not 
   work with different versions of Gromacs.

2. Descriptive and concise docstrings in all modules. Decide on essentials
   for docstrings.

3. Figure out how to grab source code docstrings into Sphinx documentation.

3. Write a simple example Document example.

Parameter Fitting
^^^^^^^^^^^^^^^^^

1. ``ddG_MC2004``

    1. Storage of entire matrices for fij is not necessary because they are
       very sparse (only ~10 nonzero entries). Instead store just the indices
       and values of nonzero entries to file reading overhead.
    2. Is there a way to determine the noise in the f_ij calculation?
    3. Is there

2. ``FRET``

    1. Decide on procedure and format for target data. 
    2. Integrate into the procedure for solving for the solutions.

3. ``RMSF``

    1. Write ``compute_Jacobian.py``
    2. Decide on procedure and format for target data. 
    3. Integrate into the procedure for solving for the solutions.

4. ``solver.py``

    1. Integrate the "TSVD" and "Cplex" options into the procedure.

Analysis
^^^^^^^^

1. Integrate ``bootstrap.py`` to calculate errors on WHAM free energy curves
   by bootstrapping.

    1. 
