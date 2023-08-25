# Hess-Smith Panel Method
The Smith-Hess panel method is an incompressible potential flow method that allows for the determination of the surface pressure over a surface. The method divides the surface of the object (e.g., an aircraft wing) into a series of small panels or flat segments. Each panel is associated with both a source and a vortex. A source panel represents a distribution of sources that induce velocity at the panel's centroid, creating a normal velocity on the panel. A vortex panel represents a distribution of vortices along the panel, which induces circulation and tangential velocity.

Essentially we wish to solves Laplace's equation $\nabla^2\phi=0$ for the flow potential $\phi$ over the surface we are interested in. The flow potential is two-dimensions is
$$
\phi=\phi_{\infty}+\int_{S}\left[\frac{q\left(s\right)}{2\pi}\ln r-\frac{\gamma\left(s\right)}{2\pi}\theta\right]ds,
$$
where $\phi_{\infty}$ is the far-field potential, $q$ is the 2D source strength, $\gamma$ is the vortex singularity strength, $\theta=\arctan(y/x)$ and $r$ is blah.

If we then discretize our surface into $N$ panels, then we may approximate the potential by 
$$
\phi=\phi_{\infty}+\sum_{j=1}^{N}\int_{P_{j}}\left[\frac{q_{j}}{2\pi}\ln r_{j}-\frac{\gamma}{2\pi}\theta_{j}\right]ds_{j},
$$
where the integral in the sum is now a surface integral of the panel $P_j$. Laplace's equation can then be solved with the included condition that we equate velocity components tangential to the panels adjacent to the trailing edge on the upper and lower surface. This condition is known as the Kutta condition. Essentially we have a linear system of $N+1$ variables and $N+1$ unknowns. The unknowns are the source and vortex strengths which we place inside a vector labelled $\mathbf{x}$ where $x_i=q_i$ for $i=1,\dots,N$ and $x_{N+1}=\gamma$. The coefficients to this linear system are the elements of the $(N+1)\times (N+1)$ matrix $\mathbf{A}$.

The inner $N\times N$ matrix elements are given by [[1]](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://archive.aoe.vt.edu/mason/Mason_f/CAtxtChap4.pdf)
$$A_{i,j}=\frac{1}{2\pi}\sin\left(\theta_{i}-\theta_{j}\right)\ln\left(\frac{r_{i,j+1}}{r_{i,j}}\right)+\frac{1}{2\pi}\cos\left(\theta_{i}-\theta_{j}\right)\beta_{i,j}.$$

The $(N+1)^{\text{th}}$ column entries are given by 

$$
A_{i,N+1}=\frac{1}{2\pi}\sum_{j=1}^{N}\left[\cos\left(\theta_{i}-\theta_{j}\right)\ln\left(\frac{r_{i,j+1}}{r_{i,j}}\right)-\sin\left(\theta_{i}-\theta_{j}\right)\beta_{i,j}\right].
$$ 

The first $N$ elements of the vector $x$ are given by 
$$
b_{i}=V_{\infty}\sin\left(\theta_{i}-\alpha\right).
$$

We then apply the Kutta condition by equating velocity components tangential to the panels adjacent to the trailing edge on the upper and lower surface. This leads to 

$$ 
A_{N+1,j}=\frac{1}{2\pi}\bigg[\sin\theta_{i,j}\beta_{1,j}+\sin\theta_{N,J}\beta_{N,j}
-\cos\theta_{i,j}\ln\left(\frac{r_{1,j+1}}{r_{1,j}}\right)-\cos\theta_{N,j}\ln\left(\frac{r_{N,j+1}}{r_{N,j}}\right)\bigg].\nonumber
$$

And finally the bottom right element is

$$
A_{N+1,N+1}=\frac{1}{2\pi}\sum_{j=1}^{N}\bigg[\sin\theta_{1,j}\ln\left(\frac{r_{1,j+1}}{r_{1,j}}\right)+\sin\theta_{N,j}\ln\left(\frac{r_{N,j+1}}{r_{N,j}}\right)+\cos\theta_{1,j}\beta_{1,j}+\cos\theta_{N,j}\beta_{N,j}\bigg].\nonumber
$$

The final element of $b$ is 
$$b_{N+1}=-V_{\infty}\cos\left(\theta_{1}-\alpha\right)-V_{\infty}\left(\cos\theta_{N}-\alpha\right).$$

We may now solve the linear system $\mathbf{A}\mathbf{x}=\mathbf{b}$ for the vector $\mathbf{x}$. Once this solution is found we may compute the lift coefficient. The method for determining all relevant matrices and vectors are outlined thoroughly in Section 3.16 of the Fundamentals of Aerodynamics by Anderson [[2]](https://books.google.co.uk/books?hl=en&lr=&id=5oVvEAAAQBAJ&oi=fnd&pg=PR2&dq=Fundamentals+of+Aerodynamics&ots=7xS9OVEniC&sig=b0kI42koZynfvb_Z4hSFYXo8pYw&redir_esc=y#v=onepage&q=Fundamentals%20of%20Aerodynamics&f=false).
