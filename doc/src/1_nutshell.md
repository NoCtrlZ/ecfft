# FFT Nutshell

## Abstract
This documentation explains [`ECFFT`](https://arxiv.org/pdf/2107.08473.pdf). To understand ECFFT, it's necessary to decompose the normal FFT components so starting with FFT and explaining about ECFFT after that.

## Premise
This documentation assumes that normal FFT is Cooley-Tukey FFT algorithm expressed as following.

$$
P(X) = P_0(X^2) + XP_1(X^2)
$$

\\( P_0 \\) and \\( P_1 \\) is respectively even and odd polynomials. It has \\( X^2 \\) as polynomial evaluation argument. Let's think it square map and rewrite formula as following.

$$
P(X) = P_0(ψ(X)) + XP_1(ψ(X))
$$

Define \\( ψ: X \rightarrow X^2 \\) in normal FFT. ECFFT uses different map.

## Components
In above formula, there are two methods regarding as component. One method is that the way to split \\( P \\) into half degree polynomial \\( P_0 \\) and \\( P_1 \\), and the other one is that the way to define the domain \\( X \\). Let former as **Polynomial Decomposition** and latter as **Low Degree Extension**. The terminologies are intended to be aligned with the [original paper](https://www.math.toronto.edu/swastik/ECFFT1.pdf).
