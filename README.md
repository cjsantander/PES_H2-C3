# PES_H2-C3

This repository contains the source code for the four dimensional potential of our paper:

[Deexcitation rate coefficients of C3 by collision with H2 at low temperatures](https://doi.org/10.1051/0004-6361/202142434)<br/>
Carlos Santander, Otoniel Denis-Alpizar and Carlos Cárdenas<br/>

```
@article{Santander2021,
  title = {Deexcitation rate coefficients of C3 by collision with H2 at low temperatures},
  author = {C. Santander and O. Denis-Alpizar and Carlos C{\'{a}}rdenas},
  doi = {10.1051/0004-6361/202142434},
  url = {https://doi.org/10.1051/0004-6361/202142434},
  year = {2021},
  month = nov,
  publisher = {{EDP} Sciences},
  journal = {Astronomy {\&} Astrophysics}
}
```

## Requirements

The code has been tested with GNU Fortran (Homebrew GCC 11.2.0), GNU Fortran (Ubuntu 9.3.0-17ubuntu1~20.04) and f2py to compile for Python

## Compile

To compile to use in Python, run

```Shell
f2py -m pot -c potential.f flegan.f
```

If it works, run the following command to print the global minimum of the PES

```Shell
python -c "import pot; print(pot.potential(4.76, 0, 0, 0))"
```
