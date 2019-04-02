=====================
Chimera scipion plugin
=====================

Plugin to use chimera programs within the Scipion framework: rigit fit, open chimera, model using an atomic struct reference and find contacts between atoms.


- **TUTORIAL:** https://github.com/I2PC/scipion/wiki/User-Documentation#tutorials

=====
Setup
=====

- **Install this plugin:**

.. code-block::

    scipion installp -p scipion-em-chimera

OR

  - through the plugin manager GUI by launching Scipion and following **Configuration** >> **Plugins**

Alternatively, in devel mode:

.. code-block::

    scipion installp -p local/path/to/scipion-em-chimera --devel

- **TESTS:**
    - ./scipion test chimera.tests.test_protocol_chimera_operate
    - ./scipion test chimera.tests.test_protocol_chimera_fit
    - ./scipion test chimera.tests.test_protocol_modeller_search
    - ./scipion test chimera.tests.test_protocol_contact


![build status](http://arquimedes.cnb.csic.es:9980/badges/chimera_devel.svg)