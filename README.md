# Transportation-Network-Case
Code used for a case study regarding (fictional) data on a transportation network


## Problem Description and Data

We are tasked with connecting a network of cities. The network must form **one connected component**, and **Eindhoven** must always be connected.  

We are provided with five data files:

- **cities.txt** – Information about all possible cities.  
- **connections.txt** – All connections between cities with costs for each transport modality.  
- **player.txt** – Budgets available for each transport modality.  
- **targets.txt** – Sets of two cities that yield revenue when connected.  
- **weights.txt** – Weights used to calculate the final score in Case 3.

In total, to evaluate networks the following KPI’s are used: Profit, Emission, and Reach. Profit will be
calculated by subtracting the total investment from the total revenue:

$$Profit = Total\ Revenue - Total\ Investment$$

The total revenue is simply the summation of the revenues of the targets that are contained in a network.
Let $(i, j)$ denote a target, $R_{i,j}$ the revenue of this target and let $C$ denote the set of cities connected in a
network, then:

$$Total\ Revenue = \sum_{i,j \in C} R_{i,j}$$

Investments can be made to increase the available capacity for each modality separately. Investments are
made per 5 units up to a maximum capacity increase of 25. The investment costs differ per each modality.
Let ki be the amount of times we invest into an additional 5 unit of modality i. Then total cost is calculated
as follows:
- **Truck:** $C_t(k_t) = 3 * k_t$
- **Airplane:** $C_a(k_a) = k_a^2$  
- **Ship:** $C_s(k_s) = C_s(k_s-1) + k_s$, with $C_s(1) = 1$  
- **Train:** $C_r(k_r) = 2^(k_r-1)$

Then, total investment is equal to:
$$Total\ Investment = C_t(k_t) + C_a(k_a) + C_s(k_s) + C_r(k_r)$$

When finding a network, the type of modality has to be decided. The transportation types considered
are trucks, ships, planes, and trains. Each unit used of the transport modalities will produce some level of
emission. This level is displayed in the following table:

| Transport modality | Emission per unit |
|------------------|-----------------|
| Truck            | 8               |
| Airplane         | 25              |
| Ship             | 1               |
| Train            | 3               |

For a given network, let $c_t$ denote the total capacity used for trucks, let $c_a$ denote the total capacity used
for airplanes, let $c_s$ denote the total capacity used for ships, and let $c_r$ denote the total capacity for trains.
Then,
$$Total\ Emission = 8*c_t + 25*c_a + 1*c_s + 3*c_r$$

Finally, the reach is the number of regions that are contained in a network. Then, using these three
KPI’s, the performance of a network will be computed as (let $w_P$ ,$w_E$,$w_R$ denote the weights for Profit,
Emission, and Reach respectively): 

$$Score = w_P * P - w_E * E + w_R * R$$

---

## General Use Algorithms and Methods

### Steiner Trees

The shortest route between a set of terminals (including Eindhoven and targets) can be solved as a **minimum cost Steiner tree problem**.  

**Algorithm: Minimum Cost Steiner Tree**

1. Start with set of terminals `T = {Eindhoven}`.  
2. While `T` does not include all required terminals:  
   - Select a terminal not in `T` with minimum distance (Dijkstra) to any city in `T`.  
   - Add the shortest path connecting this terminal to `T`.

### Graph Construction

We construct a graph with one edge per connection. Edge weights are **the average of all modality costs**:

> Example: connection costs = 8, 9, 9, 10 → weight = (8 + 9 + 9 + 10)/4 = 9  

### Modality Choice

A **greedy algorithm** is used to select the transport modality for each connection:

**Algorithm: Modality Choice**

1. Initialize capacities.  
2. For each connection in the Steiner tree:  
   - Choose modality with lowest transport units.  
   - If capacity allows, use it.  
   - Else, invest in capacity expansion.  
   - If investment fails, try next lowest-cost modality.  
3. Return the obtained network.

> We verify networks with a `verifyEdges` function to ensure all connections are included. If not, that Steiner tree is discarded.

---

## Case 1: Proof of Concept

Goal: Minimize emissions while achieving at least **50 million profit**.

**Algorithm: Proof of Concept**

1. Start with terminals `T = {Eindhoven}`.  
2. Sort `targetList` by revenue.  
3. While `Profit < 50`:  
   - Take next feasible target in `targetList`.  
   - Add origin and destination to `T`.  
   - Find corresponding Steiner tree `ST`.  
   - Calculate profit for network `N`.  
4. Output network `N`.

**Special Adjustment:** Edge weights penalize trucks and airplanes by a factor of 2 to reduce emissions.

**Results:**

- Profit: 52 million  
- Emission: 907  
- Reach: 8 regions  

---

## Case 2: Minimum Viable Product

Goal: Maximize profit given that cities are connected from all 10 regions.

**Algorithm: Minimum Viable Product**

1. Start with terminals `T = {Eindhoven}`.  
2. While `Reach < 10`:  
   - Find target with highest profit in an unconnected region.  
   - Add origin city to `T`.  
3. While profitable targets exist:  
   - Add targets to `T` to maximize profit.  
4. Output network using Steiner tree of `T`.

**Results:**

- Profit: 72 million  
- Emission: 2178  
- Reach: 10 regions  

---

## Case 3: Maximum Score

Goal: Maximize **Score = 43*P - 1*E + 292*R**.

**Algorithm: Greedy + Randomized Local Search**

1. Start with terminals `T = {Eindhoven}`.  
2. While score `X` can be increased:  
   - For each target `(i,j)` not in `T`:  
     - Create `T' = T ∪ {i,j}`  
     - Calculate Steiner tree `ST'` and network `N'`  
     - Compute score `X'`  
   - Select target that maximally increases `X`.  
   - Update `T, ST*, N*, X*`.

**Randomized Local Search:**

1. Remove elements from `T` with probability `2/|T|`.  
2. Add remaining targets with probability `1/(total targets - |T|)`.  
3. Compute new network and score.  
4. If score improves, update `T` and `X`.  
5. Repeat 100x, take the best solution.

**Results:**

- Initial greedy score: 3945  
- Final network score: 4044  
- Profit: 56 million  
- Emission: 1284  
- Reach: 10 regions  

---
