# Percolation on spatial networks with triadic regulations

Codes for triadic percolation on spatial networks generated by the Waxman model. Written by Hanlin Sun (<hanlinsun.work@gmail.com>) and Ana P Millán (<ana@onsager.ugr.es>).

This repository contains codes for generating the time series of spatial patterns and the order parameter (the size of the giant component) implemented in Julia and Python. The codes have been tested on networks with up to 10,000 nodes. 

Link to the paper: https://doi.org/10.1093/pnasnexus/pgae270

 Julia codes: 

- `spatial.jl`: The codes include all the functions used for producing orbit diagrams and time series of spatial patterns from Monte Carlo simulation.

Python codes:
- `xxx.py`: xxxx



# How to use
The Julia codes are implemented in Julia 1.10. Packages used are `Plots`, `StatsBase`, `DelimitedFiles` and `LaTexStrings`.

The Python codes requires packages `xxx`.

# Citing
If you find the codes useful in your research, please cite the following paper:

```latex

@article{10.1093/pnasnexus/pgae270,
    author = {Millán, Ana P and Sun, Hanlin and Torres, Joaquín J and Bianconi, Ginestra},
    title = "{Triadic percolation induces dynamical topological patterns in higher-order networks}",
    journal = {PNAS Nexus},
    pages = {pgae270},
    year = {2024},
    month = {07},
    abstract = "{Triadic interactions are higher-order interactions which occur when a set of nodes affects the interaction between two other nodes. Examples of triadic interactions are present in the brain when glia modulate the synaptic signals among neuron pairs or when interneuron axo-axonic synapses enable presynaptic inhibition and facilitation, and in ecosystems when one or more species can affect the interaction among two other species. On random graphs, triadic percolation has been recently shown to turn percolation into a fully-fledged dynamical process in which the size of the giant component undergoes a route to chaos. However, in many real cases, triadic interactions are local and occur on spatially embedded networks. Here we show that triadic interactions in spatial networks induce a very complex spatio-temporal modulation of the giant component which gives rise to triadic percolation patterns with significantly different topology. We classify the observed patterns (stripes, octopus and small clusters) with topological data analysis and we assess their information content (entropy and complexity). Moreover we illustrate the multistability of the dynamics of the triadic percolation patterns and we provide a comprehensive phase diagram of the model. These results open new perspectives in percolation as they demonstrate that in presence of spatial triadic interactions, the giant component can acquire a time-varying topology. Hence, this work provides a theoretical framework that can be applied to model realistic scenarios in which the giant component is time-dependent as in neuroscience.}",
    issn = {2752-6542},
    doi = {10.1093/pnasnexus/pgae270},
    url = {https://doi.org/10.1093/pnasnexus/pgae270},
    eprint = {https://academic.oup.com/pnasnexus/advance-article-pdf/doi/10.1093/pnasnexus/pgae270/58483079/pgae270.pdf},
}

```
# License
This code can be redistributed and/or modified under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  
This program is distributed by the authors in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

