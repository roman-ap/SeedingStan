The primary goal of this project is to illustrate the need of seeding HMC algorithms in Stan when dealing with latent variables.  

Let us assume that a population experience a population growth as follows:

$$\begin{aligned}
M_{t} &\sim \mathrm{Bin}\left( N_{t}^{pre}, 1-e^{-m_{t}} \right)\\
A_{t} &\sim \mathrm{Bin}\left( M_{t}, \frac{a_{t}}{m_{t}} \right)\\
N_{t}^{post} &= N_{t}^{pre} - M_{t}\\
N_{t+1}^{pre} &\sim \mathrm{Poi}\left( e^{f_{t}}N_{t}^{post} \right)
\end{aligned}$$

where $m_{t}=a_{t}+b_{t}$.
