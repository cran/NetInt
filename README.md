# NetInt
This repository contains different network integration methods that can be classified into:

**Unweighted approaches**
These methods perform network integration without considering the "predictiveness" (i.e. 
the informativeness of each network) with respect to a prediction task. In particular,
the following integrations are implemented:
- *Unweighted Average (UA)*
- *Per-edge Unweighted Average (PUA)*
- *Maximum (MAX)*
- *Minimum (MIN)*
- *At least K (ATLEASTK)*


**Weighted approaches**
These methods require to provide as input a weight for each network, which are usually 
learned considering an appropriate learning algorithm:
- *Weighted Average Per-class (WAP)*
- *Weighted Average (WA)*


**Citation** - 
These methods were presented in the paper:

Valentini, Giorgio, et al. "An extensive analysis of disease-gene associations using network integration and 
fast kernel-based gene prioritization methods." Artificial Intelligence in Medicine 61.2 (2014): 63-78.

Corresponding bib entry:
```
@article{valentini2014extensive,
  title={An extensive analysis of disease-gene associations using network integration and fast kernel-based gene prioritization methods},
  author={Valentini, Giorgio and Paccanaro, Alberto and Caniza, Horacio and Romero, Alfonso E and Re, Matteo},
  journal={Artificial Intelligence in Medicine},
  volume={61},
  number={2},
  pages={63--78},
  year={2014},
  publisher={Elsevier}
  
}
```
