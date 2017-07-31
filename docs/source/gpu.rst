.. _gpu:

*****
GPGPU
*****

The most computationally intensive parts of gprMax, which are the FDTD solver loops, can optionally be executed using General-purpose computing on graphics processing units (GPGPU). This has been achieved through use of the NVIDIA CUDA programming environment, therefore a `NVIDIA CUDA-Enabled GPU <https://developer.nvidia.com/cuda-gpus>`_ is required to take advantage of the GPU-based solver.

Extra installation steps for GPU usage
======================================

The following steps provide guidance on how to install the extra components to allow gprMax to run on your GPU:

1. Install the `NVIDIA CUDA Toolkit <https://developer.nvidia.com/cuda-toolkit>`_. You can follow the Installation Guides in the `NVIDIA CUDA Toolkit Documentation <http://docs.nvidia.com/cuda/index.html#installation-guides>`_
2. Install the pycuda Python module. Open a Terminal (Linux/macOS) or Command Prompt (Windows), navigate into the top-level gprMax directory, and if it is not already active, activate the gprMax conda environment :code:`source activate gprMax` (Linux/macOS) or :code:`activate gprMax` (Windows). Run :code:`pip install pycuda`


Running gprMax using GPU(s)
===========================

Open a Terminal (Linux/macOS) or Command Prompt (Windows), navigate into the top-level gprMax directory, and if it is not already active, activate the gprMax conda environment :code:`source activate gprMax` (Linux/macOS) or :code:`activate gprMax` (Windows)

Run one of the test models:

.. code-block:: none

    (gprMax)$ python -m gprMax user_models/cylinder_Ascan_2D.in -gpu


Combining MPI and GPU usage
---------------------------

Message Passing Interface (MPI) has been utilised to implement a simple task farm that can be used to distribute a series of models as independent tasks. This is described in more detail in the :ref:`OpenMP, MPI, HPC section <openmp-mpi>`. MPI can be combined with the GPU functionality to allow a series models to be distributed to multiple GPUs on the same machine (node). For example, to run a B-scan that contains 60 A-scans (traces) on a system with 4 GPUs:

.. code-block:: none

    (gprMax)$ python -m gprMax user_models/cylinder_Bscan_2D.in -n 60 -mpi 5 -gpu

.. note::

    The argument given with `-mpi` is number of MPI tasks, i.e. master + workers, for MPI task farm. So in this case, 1 master (CPU) and 4 workers (GPU cards).

