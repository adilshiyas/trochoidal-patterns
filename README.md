# Trochoidal Paths for a Multi-Robot Swarm

This work extends an existing distributed consensus protocol that enables a connected swarm of agents to generate trochoidal motion patterns, while additionally satisfying a set of geometric and temporal constraints.  
These constraints include bounds on inter-agent distance, distance to trochoid centers, and velocity limits to ensure the generated trajectories are physically trackable.

This repository contains the MATLAB implementation accompanying the paper [**Design of planar collision-free trochoidal paths for a multi-robot swarm**](https://www.sciencedirect.com/science/article/abs/pii/S0947358024002036). The proposed method extends an existing distributed consensus protocol by incorporating geometric and temporal constraints that guarantee physically trackable, collision-free trochoidal trajectories.

The design methodology for enforcing these constraints is discussed in detail in the accompanying paper. The proposed trajectories are validated both in simulation and on an indoor mobile robot platform.

![demo](videos/trochoids.gif)

<p align="center">
  Trochoidal trajectories executed by a 3-robot indoor swarm (16× speed)
</p>

## Highlights

- Decentralized multi-robot trajectory generation
- Consensus-based control with geometric and temporal constraints
- Collision-free trochoidal trajectories
- Validated on an indoor mobile robot platform


## Repository Contents

- MATLAB scripts for trajectory generation and simulation  
- Visualization scripts for solution space and trajectory analysis  



