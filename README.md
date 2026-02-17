The primary goal of this project is to illustrate the need of seeding HMC algorithms in Stan when dealing with latent variables.  

Let us assume that a population experience a population growth as follows:

$$\begin{aligned}
M_{t} &\sim \mathrm{Bin}\left( N_{t}^{pre}, 1-e^{-m_{t}} \right)\\
A_{t} &\sim \mathrm{Bin}\left( M_{t}, \frac{a_{t}}{m_{t}} \right)\\
N_{t}^{post} &= N_{t}^{pre} - M_{t}\\
N_{t+1}^{pre} &\sim \mathrm{Poi}\left( e^{f_{t}}N_{t}^{post} \right)
\end{aligned}$$

where $m_{t}=a_{t}+b_{t}$,

$$\begin{aligned}
\log(N^{pre}_{1}) &\sim \mathrm{N}\left( m_{LN0}, sd_{LN0} \right)\\
\log(a_{t}) &\sim \mathrm{N}\left( \mu_{la}, \sigma_{la} \right)\\
\log(b_{t}) &\sim \mathrm{N}\left( \mu_{lb}, \sigma_{lb} \right)\\
\log(f_{t}) &\sim \mathrm{N}\left( \mu_{lf}, \sigma_{lf} \right)
\end{aligned}$$


$$\begin{aligned}
\mu_{la} &\sim \mathrm{N}\left( m_{la}, sd_{la} \right)\\
\sigma_{la} &\sim \mathrm{N}\left( 0, sd_{la} \right)\\
\mu_{lb} &\sim \mathrm{N}\left( m_{lb}, sd_{lb} \right)\\
\sigma_{lb} &\sim \mathrm{N}\left( 0, sd_{lb} \right)\\
\mu_{lf} &\sim \mathrm{N}\left( m_{lf}, sd_{lf} \right)\\
\sigma_{lf} &\sim \mathrm{N}\left( 0, sd_{lf} \right)\\
\end{aligned}$$
