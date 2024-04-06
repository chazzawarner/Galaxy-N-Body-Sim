# Low-Compute Galaxy Merger Simulations
---
Generates a galaxy and simulates its movement using an N-body sim that utilises the Barnes-Hut algorithm. Galactic potential-density pairs are used to initalise the galaxy. Example galaxies (including the Milky Way and Andromeda) are available to play with.

Python acts as a backend, generating a .csv file of positions that the frontend can read to visualise the sim using the HTML canvas element. Currently, frontend support is broken as it does not handle 3D simulations. However, animations and plots can be generated by the backend for visualisation that do have support for handling 3D.

This was made as part of a project for a university module. The code requires some tweaks for more accurate simulations but this is a working prototype.