C implementations of all the algorithms discussed in my bachelor's degree thesis in physics at the [University of Trento](https://www.unitn.it/en), academic year 2023/2024. My supervisor was professor [Alessandro Roggero](https://webapps.unitn.it/du/it/Persona/PER0016084).

# Efficient quantum circuits for classical decision problems

### Abstract

We go over the basics of quantum computing and error correcting codes; we focus our attention on classical computations and show that these can be implemented on quantum computers using only NOT, CNOT and *Toffoli* gates, with the latter being particularly challenging to implement in a fault tolerant manner. Consequently, we set the goal of minimizing the number of *Toffoli* gates required to implement a classical decision problem. We leverage general properties of boolean functions to establish upper bounds and develop a greedy algorithm that manipulates algebraic expressions representing the targeted decision problem, with the same goal.

### Conclusions

We started our discussion talking about how fault tolerant computation is achieved on quantum computers despite the presence of noise. We presented Shor's code and showed that transversal gates are easier to implement fault tolerantly. As the CNOT gate is transversal in most common error correcting codes while the *Toffoli* gate isn't, we concluded that fault tolerant *Toffoli* gates are more expensive in terms of noise. We set the goal of reducing the number of *Toffoli* gates required to implement a classical decision problem on a quantum computer that can utilize NOT, CNOT, and *Toffoli* gates, and with an infinite supply of *ancilla qubits*. 

Then we showed that finding an expression which contains $P$ products for a function $f\in\mathbb{B}_n$ is a direct way of proving that the corresponding decision problem can be implemented with no more than $P$ *Toffoli* gates. This is why we looked into normals forms: these gave us upper bounds to the number of *Toffoli* gates of the order $\mathcal{O}(n2^n)$, where $n$ is the number of input qubits. Afterwards we managed to improve the theoretical upper bound by leveraging the distributive property of the logical AND over the logical XOR. What we got was the upper bound $\text{Tof}(n) \leq 2^{n-1} - 1$. 

We than looked into general algorithms that can lower the upper bound for a given function. We built on top of the Fast Möbius Transform and obtained a greedy algorithm that recursively collects terms with complexity scaling as $\mathcal{O}(n^2 2^n)$. As showed by Table 3, we demonstrated statistically that the greedy approach performs better on average than using the normal forms directly or collecting variables in a fixed order. We also showed the limitations of working on algebraic expressions instead of circuits by considering structured functions which result in high counts of products most of which repeat. 
