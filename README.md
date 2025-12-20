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

<video width="720" height="405" controls autoplay muted loop playsinline>
<source src='https://github.com/media/adilshiyas/adilshiyas.github.io/blob/master/Troch_16x_PV (1).mp4?raw=true' type="video/mp4">
 Your browser does not support the video tag.
</video>

<video controls="" width="800" height="500" muted="" loop="" autoplay="">
<source src="https://github.com/adilshiyas/adilshiyas.github.io/raw/main/Troch_16x_PV (1).mp4" type="video/mp4">
</video>

<video controls="" width="800" height="500" muted="" loop="" autoplay="">
<source src="https://github.com/adilshiyas/adilshiyas.github.io/Troch_16x_PV (1).mp4" type="video/mp4">
</video>

<video controls="" width="800" height="500" muted="" loop="" autoplay="">
<source src="https://github.com/adilshiyas/adilshiyas.github.io/blob/main/Troch_16x_PVs.mp4" type="video/mp4">
</video>

https://github.com/adilshiyas/adilshiyas.github.io/assets/153742460/f4bb2ff9-1873-4609-8347-361bcf31748e

<p align="center">
  Trochoidal trajectories executed by a 3-robot indoor swarm (16× speed)
</p>

<p align="center">
  <video width="80%" controls muted loop playsinline>
    <source src="https://github.com/adilshiyas/trochoidal-patterns/raw/main/videos/Troch_16x_PVs.mp4" type="video/mp4">
  </video>
  <br>
  Trochoidal trajectories executed by a 3-robot indoor swarm (16× speed)
</p>

