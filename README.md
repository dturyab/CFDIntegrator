<img width="65" height="21" alt="image" src="https://github.com/user-attachments/assets/03bfba56-104e-4e13-9268-636fdc627a8c" /># CFDIntegrator
There is a C++ code for numerical integrating of hyperbolic systems of the form: $\frac{\partial \vec U}{\partial t}+\frac{\partial \vec F}{\partial x}=\vec S(\vec U, x)$.
## Numerical method
The developed program implements the Godunov method with TVD modification. The Godunov method may be formulatted in the form:
$$
\vec U_i^{n+1}=\vec U_i^n+\frac{\Delta t}{\Delta x}\left ( \vec F_{i-\frac{1}{2}}-\vec F_{i+\frac{1}{2}} \right ),
$$
where $dt$ - time step, $dx$ - space step, fractional indices $i\pm \frac{1}{2}$ mean edges of the cell with center in $x_i$.
So $F_{i\pm\frac{1}{2}}$ - solutions of the Riemann problem on the left and right edges of the cell respectively, these may be found with approximate Riemann solver. 
Current program uses HLLC[1] Riemann solver on the edges.

To increase converges rate up to 2nd TVD scheme is implemented:
1. Piece-wise linear reconstructon of the values inside each cell:
   $$
   U(t^k,x)=U_j^k+Q_j^k(x-x_j),x\in [x_j-\frac{1}{2}\Delta x,x_j+\frac{1}{2}\Delta x];
   $$
3. Predictor step:
   $$
   \frac{\hat{U}_j^{k+1}-U_j^k}{\frac{1}{2} \Delta t}+\frac{F(U_j^k+\frac{1}{2}\Delta xQ_j^k)-F(U_j^k-\frac{1}{2}\Delta xQ_j^k)}{\Delta x}=0;
   $$
4. Corrector step:
   $$
   \frac{\vec U_j^{k+1}-\vec U_j^k}{\Delta t}+\frac{\vec F_{j+\frac{1}{2}}-\vec F_{j-\frac{1}{2}}}{\Delta x}=0.
   $$

## Comprasion on exact solutions
### Sod shock tube comprasion:
![Shock tube](https://github.com/dturyab/CFDIntegrator/blob/main/img/SodExactNumericCompr.png)
Convergence to numerical solution(N=1600):
|N|dx|Dens. $L_1$|Dens. q|Vel. $L_1$|Vel. q|Pres. $L_1$|Pres. q|
|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
|1600|6.25E-04|||||||
|800|1.25E-03|0.00046769|1.45|0.000635|1.42|0.000342|1.52|
|400|2.5E-03|0.0012862|1.18|0.001701|1.24|0.000987|1.24|
|200|5E-03|0.002933|1.25|0.0040336|1.31|0.002152|1.23|
|100|1E-02|0.0070236||0.0100611||0.005071||

### Sedovproblem comprasion:
![Sedov problem](https://github.com/dturyab/CFDIntegrator/blob/main/img/SedovExactNumericCompr.png)
Convergence to exact solution:
|N|dx|Dens. $L_1$|Dens. q|Pres. $L_1$|Pres. q|
|:--:|:--:|:--:|:--:|:--:|:--:|
|5000|1.98E-05|7.12E-05||1.54E-05||
|2000|4.95E-05|2.8E-04|1.49|4.31E-05|1.12|
|1000|9.9E-05|7.2E-04|1.44|9.86E-05|1.15|
|500|1.98E-04|1.36E-03|1.28|1.71E-04|1.04|

1. Eleuterio F. Toro Riemann Solvers and Numerical Methods for Fluid Dynamics. - New York: Springer. - 2009, - 724 p.
