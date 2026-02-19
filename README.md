The primary goal of this project is to illustrate the need of seeding HMC algorithms in Stan when dealing with unobserved/latent variables as is the case of state-space models.  

For a discrete-time setting, let us assume that a population experience an exponential growth with time-varying rate $r_{t}$, that is, $N_{t+1}=e^{r_{t}}N_{t}$ for $t=1,2,\ldots,T-1$. We can rewrite the time-varying rate as $r_{t}=f_{t}-m_{t}$, where $f_{t},m_{t}\in\mathbb{R}^{+}$ are the fertility and mortality rates at time $t$, respectively. Furthermore, we can also rewrite the time-vaying mortality rate as $m_{t}=a_{t}+b_{t}$, where $a_{t}\in\mathbb{R}^{+}$ is the rate associated to a specific  mortality cause for which we have measuments whilst $b_{t}\in\mathbb{R}^{+}$ is the corresponding background mortality rate. Now, at each time $t$, we can think of abundance before and after mortality, say $N^{pre}_{t}$ and $N^{post}_{t}$, respectively. Thus, if $M_{t}$ denotes the total mortality we have that $N_{t}^{post}=N_{t}^{pre} - M_{t}$. Finally, let $A_{t}$ be the observed mortality whose underlying rate is $a_{t}$.

Given the aforementioned population demographics, we can propose the following stochastic population mechanism:

$$\begin{aligned}
M_{t} &\sim \mathrm{Bin}\left( N_{t}^{pre}, 1-e^{-m_{t}} \right)\\
A_{t} &\sim \mathrm{Bin}\left( M_{t}, \frac{a_{t}}{m_{t}} \right)\\
N_{t+1}^{pre} &\sim \mathrm{Poi}\left( e^{f_{t}}N_{t}^{post} \right)
\end{aligned}$$

From a Bayesian perspective

$$\begin{aligned}
\log(N^{pre}_{1}) &\sim \mathrm{N}\left( m_{LN0}, sd_{LN0} \right)\\
\log(a_{t}) &\sim \mathrm{N}\left( \mu_{la}, \sigma_{la} \right)\\
\log(b_{t}) &\sim \mathrm{N}\left( \mu_{lb}, \sigma_{lb} \right)\\
\log(f_{t}) &\sim \mathrm{N}\left( \mu_{lf}, \sigma_{lf} \right)
\end{aligned}$$

with

$$\begin{aligned}
\mu_{la} &\sim \mathrm{N}\left( m_{\mu_{la}}, sd_{\mu_{la}} \right)\\
\sigma_{la} &\sim \mathrm{N}\left( 0, sd_{la} \right)\\
\mu_{lb} &\sim \mathrm{N}\left( m_{lb}, sd_{lb} \right)\\
\sigma_{lb} &\sim \mathrm{N}\left( 0, sd_{lb} \right)\\
\mu_{lf} &\sim \mathrm{N}\left( m_{lf}, sd_{lf} \right)\\
\sigma_{lf} &\sim \mathrm{N}\left( 0, sd_{lf} \right)\\
\end{aligned}$$
