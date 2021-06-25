================
Developer Guide
================


Coding Style
==============

Please stick to the following coding style for easier code readability:
    - class variables start with an underscore and lowercase letters e.g. ``_foo``
    - functions start with a capital letter e.g. ``GetSettings()``
    - any variable/function names have capital letters at each individual word e.g. ``GetAllCellsAdjacentTo(Cell i)``

Please also use the provided ``code/.clang-format`` style format to format your code before pushing your latest commits.
Some editors offer to automatically apply the style format upon saving a file (e.g. ``Qtcreator``).


Continuous Integration (CI)
============================

Every commit on the master branch will trigger a build test and unit tests.
If either of the tests fail you will get an email, telling you the 'pipeline' has failed. If this happens, visit the 'CI/CD' tab and check what went wrong and try to fix the error as soon as possible.
Creating a merge request from a branch into the master branch will (not tested yet) also trigger the tests and your merge will be rejected in case any test fails.

If you add additional libraries to the code, these also need to be added to the test environment, i.e. the respective docker container.

Little guide:

Test the library installation in the docker container

.. code-block:: bash 

    docker run -it --rm rtsn/test:latest bash

Note the steps required and add them to the `Dockerfile` in the `scripts` folder.
Build the new container (takes some time)

.. code-block:: bash 

    cd docker build -t rtsn/test:latest .

or commit your changes to the image (google that procedure).
Push the new image to `hub.docker.com`

.. code-block:: bash 
     
   docker push rtsn/test:latest

This last step requires a preceeding `docker login`. Ask Jannick for login credentials.