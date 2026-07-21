# Collinear Lagrange Points in the Circular Restricted Three-Body Problem with Radiation Pressure

[![super-linter](../../actions/workflows/super-linter.yml/badge.svg)](../../actions/workflows/super-linter.yml) ![human-only code](https://img.shields.io/badge/human--only-code-white)

A Fortran collinear Lagrange points calculator.

---

## Table of Contents

- [Background](#background)
- [Key Files](#key-files)
- [Software Requirements](#software-requirements)
- [Quality Assurance](#quality-assurance)
- [Getting Started](#getting-started)
- [References](#references)

## Background

ForCLaP (Fortran Collinear Lagrange Points Calculator) is a Fortran command line
application demonstrating Newton-Raphson convergence on collinear solutions to
the circular restricted three-body problem. Optional displacement of the
collinear Lagrange points due to radiation pressure, relevant to small bodies
such as asteroidal dust particles, is incorporated.

ForCLaP assumes a traditional three-body system with a large, central primary
(m<sub>1</sub>), a smaller secondary (m<sub>2</sub>) and an infinitesimal
tertiary mass (m<sub>3</sub>). ForCLaP accepts a user-input value for the
m<sub>1</sub>/m<sub>2</sub> mass ratio of the system. Use of a mass ratio within
the Routh value limit, is enforced however.

ForCLaP also accepts a user-input value for a ratio of solar radiation pressure
to gravitational force experienced by m<sub>3</sub>, <i>&beta;</i>. This allows
simulation of a luminous primary such as the Sun, but ForCLaP can just as easily
be used for simulation of a non-luminous primary, such as planet, with an input
ratio, <i>&beta;</i>, of zero.

Collinear Lagrange point solution calculation is done with a Newton-Raphson
algorithm. The program therefore inherits limitations of that algorithm. The
algorithm can suffer runaway solution divergence, for example, but ForCLaP
explicitly flags this and provides graceful algorithm termination. For algorithm
convergence, iteration continues until the solution and the Newton derivative
are both stable to four significant figures.

### Locations of the Collinear Lagrange Points in the Circular Restricted Three-Body Problem with Radiation Pressure

Formulas for the locations of the collinear Lagrange points
<i>L</i><sub>1</sub>, <i>L</i><sub>2</sub> and <i>L</i><sub>3</sub> are given
below. Symbol definitions and explanatory theory are given in Stenborg
[[1](#references)].

#### <i>L</i><sub>1</sub>

$$
\begin{eqnarray}
f(r_2) &=& \frac{r_2^5 - 3r_2^4 + 3r_2^3 - \beta r_2^2}{-r_2^5 + 2r_2^4 - r_2^3 + r_2^2 - 2r_2 + 1} - \frac{m_2}{m_1}\\
f'(r_2) &=& \frac{-r_2^8 + 4r_2^7 - 3\beta r_2^6 - (7 - 2\beta)2r_2^5 + (26 - \beta)r_2^4 - 24r_2^3 + (9 - 2\beta)r_2^2 - 2\beta r_2}{r_2^{10} - 4r_2^9 + 6r_2^8 - 6r_2^7 + 9r_2^6 - 12r_2^5 + 9r_2^4 - 6r_2^3 + 6r_2^2 - 4r_2 + 1}\nonumber \\
\end{eqnarray}
$$

#### <i>L</i><sub>2</sub>

$$
\begin{eqnarray}
f(r_2) &=& \frac{r_2^5 + 3r_2^4 + 3r_2^3 + \beta r_2^2}{-r_2^5 - 2r_2^4 - r_2^3 + r_2^2 + 2r_2 + 1} - \frac{m_2}{m_1}\\
f'(r_2) &=& \frac{r_2^8 + 4r_2^7 + (6 + 3\beta)r_2^6 + (14 + 4\beta)r_2^5 + (26 + \beta)r_2^4 + 24r_2^3 + (9 + 2\beta)r_2^2 + 2\beta r_2}{r_2^{10} + 4r_2^9 + 6r_2^8 + 2r_2^7 - 7r_2^6 - 12r_2^5 - 7r_2^4 + 2r_2^3 + 6r_2^2 + 4r_2 + 1}\nonumber \\
\end{eqnarray}
$$

#### <i>L</i><sub>3</sub>

$$
\begin{eqnarray}
f(r_1) &=& \frac{-r_1^5 - 2r_1^4 - r_1^3 + (1 - \beta)r_1^2 + (1 - \beta)2r_1 - \beta + 1}{r_1^5 + 3r_1^4 - 3r_1^3} - \frac{m_2}{m_1}\\
f'(r_1) &=& \frac{-r_1^6 - 4r_1^5 - (2 - \beta)3r_1^4 - (1 + \beta)14r_1^3 - (1 - \beta)(26r_1^2 + 24r_1 + 9)}{r_1^8 + 6r_1^7 + 15r_1^6 + 18r_1^5 + 9r_1^4}\nonumber \\
\end{eqnarray}
$$

## Key Files

| File          | Notes            |
| :------------ | :--------------- |
| `forclap.f90` | Fortran program. |

## Software Requirements

| Software | Notes                                                                                |
| :------- | :----------------------------------------------------------------------------------- |
| Fortran  | [Details here](https://fortran-lang.org/). Free and proprietary compilers available. |

## Quality Assurance

ForCLaP has been tested with the GNU Fortran (gfortran) compiler, version
15.1.0.

## Getting Started

Program source code requires compilation on the target system with a Fortran
compiler. E.g., for gfortran:

    gfortran forclap.f90 -o forclap

## References

1. Stenborg, TN 2008, "Collinear Lagrange Point Solutions in the Circular
   Restricted Three-Body Problem with Radiation Pressure using Fortran", in RW
   Argyle, PS Bunclark and JR Lewis (eds), Astronomical Data Analysis Software
   and Systems XVII, vol. 394, Astronomical Society of the Pacific, Provo, UT,
   pp. 734&ndash;737.\
   [View PDF](https://aspbooks.org/publications/394/734.pdf)
   &nbsp; [View at publisher](https://aspbooks.org/custom/publications/paper/394-0734.html)
   &nbsp; [SciX](https://scixplorer.org/abs/2008ASPC..394..734S/abstract)
