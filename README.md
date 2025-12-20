# Trochoidal Paths for a Multi-Robot Swarm

This repository contains MATLAB simulation code associated with the paper:  
[**Design of planar collision-free trochoidal paths for a multi-robot swarm**](https://www.sciencedirect.com/science/article/abs/pii/S0947358024002036)

## Overview

This work extends an existing distributed consensus protocol that enables a connected swarm of agents to generate trochoidal motion patterns, while additionally satisfying a set of geometric and temporal constraints.  
These constraints include bounds on inter-agent distance, distance to trochoid centers, and velocity limits to ensure the generated trajectories are physically trackable.

The design methodology for enforcing these constraints is discussed in detail in the accompanying paper. The proposed trajectories are validated both in simulation and on an indoor mobile robot platform.

## Repository Contents

- MATLAB scripts for trajectory generation and simulation  
- Visualization scripts for solution space and trajectory analysis  

## Real-World Implementation

The following videos demonstrate collision-free trochoidal trajectories executed by a multi-robot indoor mobile robot platform.
<video src="videos/Troch_16x_PVs.mp4" controls width="80%">
</video>
