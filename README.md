# A covariant representation for generalized time-frequency analysis of discrete signals

## Detection methodology based on the zeros of the *Kravchuk* spectrogram  

This project contains the Python code associated to the papier

> **Pascal, B.** & Bardenet, R. (2022). ``*A covariant, discrete time-frequency representation tailored for zero-based signal detection*". Submitted.  
>  [hal](https://hal.archives-ouvertes.fr/)

## Project description

Following a recently unorthodox path in time-frequency analysis shedding light on the spectrogram zeros, we introduce a novel generalized time-frequency representation, specifically designed for the analysis of *discrete* signals, particularly amenable to spatial statistics on the zeros thanks to its *compact* phase space.  

> This toolbox provides a stable implementation of this novel *Kravchuk* transform and the code to reproduce `Figures 1, 2 and 6` of the paper ``*A covariant, discrete time-frequency representation tailored for zero-based signal detection*", comparing the standard and the *Kravchuk* spectrograms of noisy chirps, with a peculiar focus on the zeros.  
> A demonstration is given in the notebook [`kravchuk-spectrogram-and-zeros`](/demos/kravchuk-spectrogram-and-zeros.ipynb).

A novel efficient methodology relying on the functional statistics of the point process formed by the zeros of the Kravchuk spectrogram for detecting the presence of some signal is implemented.

> The detection procedure based on the functional statistics of the zeros of the Kravchuk spectrogram is implemented.
> For sake of comparison, we provide also an implementation of the counterpart strategy relying on the zeros of the Short-Time Fourier transform developed in the paper ``*On the zeros of the spectrogram of white noise*" by Bardenet R., Flamant, J. & Chainais, P. (2021) Applied and Computational Harmonic Analysis.
>
> The interested reader can then reproduce `Figures 7, 8 and 9` of the paper ``*A covariant, discrete time-frequency representation tailored for zero-based signal detection*".  
> A demonstration is given in the notebook [`detection-test-Kravchuk-zeros`](/demos/detection-test-Kravchuk-zeros.ipynb).

## Dependencies

The following Python libraries are necessary:
- `matplotlib`
- `numpy`
- `scipy`
- `statsmodels`

Functional statistics of the pattern of zeros of the standard spectrogram are computed using [`SpatStat`](http://spatstat.org/) toolbox developed in [`R`](https://www.r-project.org/).
The incorporation of `R` functions into  `Python`  code relies on the [`spatstat-interface`](https://github.com/For-a-few-DPPs-more/spatstat-interface), developed by [G. Gautier](https://github.com/guilgautier).
