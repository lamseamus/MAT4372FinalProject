## An Exploration of Convergence in the Least Squares Monte Carlo Algorithm for American Put Options.
The complete report can be found [here](https://lamseamus.github.io/Exploration-of-Convergence-in-the-Least-Squares-Monte-Carlo-Algorithm/final_report.pdf).

## Abstract
The report presents an analysis of how varying parameters of the Least Squares Monte Carlo
(LSMC) algorithm affect the absolute convergence of the value of an estimated American put
option with respect to the underlying value derived from the binomial model. Parameters for the
models are obtained from the S&P 500 index and United States treasury bill rate data. The
number of simulation paths, the choice of polynomials, and the degree of the polynomials were
all considered as variables within the algorithm. Particularly, the effectiveness of four functions
as bases for the regression was explored: standard, Laguerre, and Hermite polynomials, in
addition to custom functions using European option payoffs. Further, additional analysis
regarding the number of possible exercise dates and the degree to which the put is In-The-Money are briefly touched upon. An At-The-Money American put was selected as a base case,
with a possibility of M=50 early exercise dates; however, analysis suggested that puts deeper
In-The-Money and with fewer possible early exercise dates would have contributed to lower
absolute errors. Results indicated that the choice of polynomial had a negligible effect on overall
convergence, while second-degree and higher degree polynomials displayed slightly reduced
absolute errors. The custom functions also outperformed both the standard and orthogonal
polynomials. A question manifesting itself throughout the study concerns how one can increase
the efficiency of simulation whilst still minimizing errors within the results. 

