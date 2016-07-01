# nek_pipe
[jcanton](mailto:jcanton@mech.kth.se)

This is a simple Poiseuille flow setup for nek5000.
The objective is to use it to test the non-linear adjoint algorithm developed by Oana Marin and Michel Schanen.

## Checkpointing
A two-layer checkpointing scheme, distinguishing between disc and memory checkpointing, is implemented in order to avoid the storage of the velocity field at every time-step.
For the memory checkpointing, the revolve algorithm is implemented.
A fixed threshold on the number of checkpoints is assumed, and the writing and reading of checkpoints is assumed to be at no cost.
The disc checkpointing implements linear checkpointing and assumes an infinite amount of space, but a limited amount of bandwidth.
