# Scaling in local to global condensation of wealth on sparse networks

## Introduction
Generalized Yard-Sale (gYS) model is a wealth exchange model we introduced in our recent study which exhibits some very fascinating dynamics. We carefully delineate our observations in [paper URL], while in here we aim to enhance visualization on what is really happening during the simulations. 

### Model
$N$ individuals start with equal wealth ($w_i=1$). Upon simulation, we update the wealth via following rule: with probability $1-p$, a randomly picked link (say, the end nodes are $i$ and $j$) transfers $\varepsilon \min(w_i,w_j)$ amount of wealth where $\varepsilon \le 1$ is a constant, and with probability $p$, the link transfers $\varepsilon w_i$ amount of wealth if the direction is $i \rightarrow j$. The tranfer direction is randomly chosen with equal probability. We do this update $N$ times during $t \rightarrow t+1$ in simulation time. 

We present simulation results on highly sparse and heterogeneous [scale-free](https://en.wikipedia.org/wiki/Scale-free_network) network with $\gamma=2.5, \langle k \rangle=4$.



## Results

### Emergence of local & global condensate
Here, we plot individual nodes' wealth $w_i$ and the wealth variance $\sigma^2$ of the total individuals ($N=97, p=10^{-6}, \varepsilon=0.05$). 

<img align="center" width="500" alt="Figs_snap" src="https://user-images.githubusercontent.com/73336039/224218576-c6390ebb-98ba-4a54-9f8e-e28576851008.png">

We note that the ***local condensate*** is characterized by almost completely frozen wealth of all nodes (video link below). It is followed by ***relaxation*** in which the frozen wealth dissolves and roam around in a random-walk-like manner across significantly longer time-scales. This transfer of wealth is driven by the aforementioned exchange scheme which occurs with probability $p$, and the average time period for each jump is proportional to $1/\varepsilon p$ as we elucidate in [paper URL].

As a consequence of the dynamic movement, and the coalescence of nearest neighbors' wealth, ***global condensate*** arises in which the total wealth is collected by a single monopolistic node. 
<br />
 
&emsp;&emsp;&emsp;&emsp; **Coalescence** &emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp;&emsp; **Transfer**

<img src="https://user-images.githubusercontent.com/73336039/224312265-6fe85e0a-4bee-4f37-bc9b-9cb0f1d98598.gif" width="500" />


In the following videos, wealth of an individual is respresented by the size of a node, the simulation time is displayed on top left, and the (time-averaged) amount of transferred wealth is visualized with an edge's brightness. 


### Initial growth ($10^1 < t < 10^4$)

[Video link](https://drive.google.com/file/d/1wDLWEYJveX0tX26qWktbd9tYN1gWEwkS/view?usp=share_link) (Google Drive)

<img src="https://user-images.githubusercontent.com/73336039/223454770-ea451211-e225-4506-82f0-6157ff966b3c.png" width="300" height="300" />


### Local condensate ($10^5 < t < 10^6$)

[Video link](https://drive.google.com/file/d/12hHhMHN0Bl-iAV__TnQnWWO1TFTtF0ef/view?usp=share_link) (Google Drive)

<img src="https://user-images.githubusercontent.com/73336039/223454798-75b65b77-4b13-40ae-becd-531accbec6a6.png" width="300" height="300" />


### Relaxation ($10^7 < t < 10^8$)

[Video link](https://drive.google.com/file/d/1-oZu167EcNQZnJlnaDbg2dm_PZQLSlCU/view?usp=share_link) (Google Drive)

<img src="https://user-images.githubusercontent.com/73336039/223454819-97870ff5-b870-4f05-974c-9acb574ba8d3.png" width="300" height="300" />


### Global condensate ($10^8 < t < 10^9$)

[Video link](https://drive.google.com/file/d/1jkMI_qgEmncgAqNV9xwKzFJ3OhkflxvZ/view?usp=share_link) (Google Drive)

<img src="https://user-images.githubusercontent.com/73336039/223454848-c463f3d9-df5d-45e1-89ec-42cc9c689242.png" width="300" height="300" />



