# FastReso 

FastReso calculates the irreducible spectral components for resonance decays for the given list of particles and decay chains.

Copyright (c) 2018-2021 Aleksas Mazeliauskas, Stefan Floerchinger, Eduardo Grossi, and Derek Teaney. 

All rights reserved.

FastReso is distributed under MIT license; see the LICENSE file available at: https://github.com/amazeliauskas/FastReso/

Steps for running on virgo3:
1) Clone this repository 
2) In your lustre directory, create a build folder
```
cd  /lustre/alice/users/<username>/FastReso/src
mkdir build
cd build
```
3) In the build directory, compile FastReso with
   ```
   cmake ..
   make
   ```
4) Enter the `analysis_FastReso.yaml` and change `<username>` with your personal user name
5) In the FastReso directory, run FastReso with `python3 submit_FastReso.py analysis_FastReso.yaml`

# Citing
  When using FastReso please refer to/cite the accompanying publication "Fast resonance decays in nuclear collisions"  Eur. Phys. J. C79 (2019) [arXiv:1809.11049] by Aleksas Mazeliauskas, Stefan Floerchinger, Eduardo Grossi, Derek Teaney
