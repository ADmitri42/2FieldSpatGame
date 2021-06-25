# 2FieldSpatGame
[![DOI](https://zenodo.org/badge/378931215.svg)](https://zenodo.org/badge/latestdoi/378931215)

This repository containts generalization of the code used in our [research](https://github.com/ADmitri42/spatial-evolutionary-game)
## spatial-evolutionary-game

In this game $2L^2$ individual "*players*" placed in a two two-dimensional spatial array.

Every individual can play one of two tactics: cooperate($\mathcal{C}$) or defect($\mathcal{D}$). At the beginning of the "game", each player has the probability of being assigned $\mathcal{C}$ with probability $p_c$ and tactics $\mathcal{D}$ otherwise.

In each round individuals "*play the game*" with 8 its neighbors and the *average cooperator* from this and another field. Tables of payoffs shown below

Table 1: Payoffs for games with neighbors.

| payoffs | $\mathcal{D}$ | $\mathcal{C}$ |
| --------- |:-------------:|----:|
| $\mathcal{D}$ | 0 | 0 |
| $\mathcal{C}$ | $b_j$ | 1 |

Table 2: Payoffs for games with *average cooperator*

| payoffs | $\mathcal{D}$ | $\mathcal{C}$ |
| --------- |:-------------:|----:|
| $\mathcal{D}$ | 0 | 0 |
| $\mathcal{C}$ | $b_j f_{ci}$ | $f_{ci}$ |

Where $f_ci$ - fraction of cooperators on the field $i$ and $b_i$ payoff on the field $i$
After playing all games the site occupied either by its original owner or by one of the neighbors who scores the highest total payoff in that round.

In general the total payoff $P(x)$ of the agent $\sigma$ with strategy $s(x)$ placed on field $j$ can be calculated as

$$P(x) = \sum_{y\in n.n.}s^T(x)\hat{H}^js(y) + \sum_{i=1}^2\lambda_{ji}s^T(x)\hat{H}^{ji}_{mf}s(y)$$ 

Where
$s(x) = (1,0)^T$ if $x$ - defector and $s(x) = (0,1)^T$ if $x$ - cooperator

Matrix $\hat{H}^j$ have a simmilar form to the Table 1

Matrix $\hat{H}^{ji}$ has a form

| $b_j f_{ci}$ | 0 |
| --------- |:-------------:|
| 0 | $f_{ci}$ |

## Setup
To compile Cython file
```bash
pip install -r requirements.txt
make
```
If you need to test 
To test
```bash
make pytest
```

## How to use it?
One can find example of the usage in notebook [example.ipynb](example.ipynb).

Two most important thing in compiled file is class
```python
TwoFieldSpatialGame
```
and function
```python
color_field_change(old_field, new_field)
```
which accept two 2-D array and returns 3-D array where
RED     - 0$\to$0
YELLOW  - 1$\to$0
GREEN   - 0$\to$1
BLUE    - 1$\to$1
