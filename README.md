libllsm2
===

Low Level Speech Model (version 2) for high-quality speech analysis/synthesis

About
---

libllsm2 is a C library providing data structures and routines for analysis (parametrization), modification and synthesis of digital speech signals.

### The model

LLSM is a two-layer model of speech. The first layer (layer 0) is a signal-level parametrization that separately models harmonic and noise (e.g. aspiration and consonants) components; the second layer (layer 1) is an acoustic-level parametrization that performs an approximated glottal inverse filtering on the harmonic component.

* libllsm2 can be viewed as a lossy speech coder, but it differs from conventional speech coders in being designed for modifications (rather than compression).

### Purpose

libllsm2 is designed for the use in speech synthesis and modification applications, where the analysis of speech signals is performed once and for all. Therefore, the library is made to analyze speech *carefully*, and synthesize speech *rapidly*. Things you can do with libllsm2 include,

* Pitch-shifting
* Time-stretching
* Modification of glottal tension (the Rd parameter in a LF model)
* Modification of aspiration noise
* Interpolation/Cross-fading

Using libllsm2 in statistical parameteric speech synthesizers is possible, but you may need to assign default values for parameters typically not modelled by those statistical systems (e.g. phases and the temporal structure of turbulent noise).

### What's new in libllsm2

libllsm2 is a rewrite of [libllsm](https://github.com/Sleepwalking/libllsm) with the following features and changes,

* Clean and well-structured API.
* Better documentation (with Doxygen) and testing.
* New feature: synthesis can be done at a different sampling rate without changing the model.
* New feature: fast and accurate speech analysis/synthesis based on Chirp-Z Transform.
* New feature: real-time synthesis.
* New feature: seamless on-the-fly switching between harmonic synthesis and pulse-by-pulse synthesis.

Compiling
---

Dependencies: [`ciglet`](https://github.com/Sleepwalking/ciglet)

Additional dependencies for the tests: [`libgvps`](https://github.com/Sleepwalking/libgvps), [`libpyin`](https://github.com/Sleepwalking/libpyin)

1. `mkdir external build`
2. Create a symbolic link to `libpyin` under `libllsm2/external`
3. Create a symbolic link to `libgvps` under `libllsm2/external`
4. Run `make single-file` in `ciglet`; create a symbolic link named `ciglet` to `ciglet/single-file` under `libllsm2/external`
5. `make`, `make test`

Note: in all of my C libraries there's a macro named `FP_TYPE`, which is either float or double. It may need to be specified as a compiler flag (i.e., `-DFP_TYPE=float`).

The test speech `test/arctic_a0001.wav` is a sample taken from the CMU Arctic database.

Documentation
---

Run `doxygen doxyfile` to generate the API documentation. Examples of using libllsm2 are given in the tests.

### Additional notes on Pulse-by-Pulse synthesis

Since the Januray 2019 update libllsm2 supports Pulse-by-Pulse (PbP) synthesis (in addition to harmonic bank and CZT based methods) of the periodic component. In PbP mode, the glottal flow of each period is generated directly in frequency domain, filtered by the layer 1 vocal tract response, and transformed to time-domain via overlap-add IFFT. This is the same as [WORLD](https://github.com/mmorise/World) vocoder except for the glottal model part.

It is easy to prove that PbP gives the same result as harmonic synthesis given stationary parameters (F0, spectral envelope and glottal parameters), and even in the case of time-varying parameters the difference is hardly perceptible. However, PbP and harmonic synthesis offer a trade-off between speed and degree of control. PbP runs a few times slower than harmonic synthesis but allows period-level parameter modification. One example usage of PbP is to convert normal speech into irregular speech (e.g. growl, see [`test-pbpeffects`](https://github.com/Sleepwalking/libllsm2/blob/master/test/test-pbpeffects.c)).

libllsm2 supports PbP synthesis to the extent that you may switch between PbP mode and harmonic mode on-the-fly during realtime synthesis by attaching `LLSM_FRAME_PBPSYN` flag to frame containers. Period-level parameter modification is possible via `LLSM_FRAME_PBPEFF` and `llsm_pbpeffect` callback structure.

Licensing
---

libllsm2 is licensed under GPLv3.

I have a pending patent on LLSM-related technologies. Under the terms of GPLv3, the patent license is granted to libllsm2 users, free from royalty.

Commerical version is available upon request (k.hua.kanru [at] ieee [dot] org).

Publications
---

Currently there's no publication directly associated with LLSM. However there is [a poster](https://github.com/Sleepwalking/prometheus-spark/blob/master/writings/pseudo-glottal-inverse-filter-hua-2016.pdf) on the pseudo glottal inverse filtering method in layer 1 LLSM.

K. Hua, "Speech Analysis/Synthesis by Non-parametric Separation of Vocal Source and Tract Responses," presented at Speech Processing Courses in Crete, 2016.

The following are the major publications that LLSM draws inspiration from.

1. G. Degottex, P. Lanchantin, A. Roebel, and X. Rodet, "Mixed source model and its adapted vocal tract filter estimate for voice transformation and synthesis," Speech Communication, vol. 55, no. 2, pp. 278–294, 2013.

2. J. P. Cabral, K. Richmond, J. Yamagishi, and S. Renals, "Glottal Spectral Separation for Speech Synthesis," IEEE Journal of Selected Topics in Signal Processing, vol. 8, no. 2, pp. 195208, 2014.

3. Y. Pantazis and Y. Stylianou, "Improving the modeling of the noise part in the harmonic plus noise model of speech," 2008 IEEE International Conference on Acoustics, Speech and Signal Processing, 2008.

4. I. Saratxaga, Hernáez I., D. Erro, E. Navas, and Sánchez J., "Simple representation of signal phase for harmonic speech models," Electronics Letters, vol. 45, no. 7, p. 381, 2009.

5. J. Bonada. "High quality voice transformations based on modeling radiated voice pulses in frequency domain," in Proceedings of Digital Audio Effects (DAFx), Vol. 3, 2004.

Test Sounds
---

`test/are-you-ready.wav` - https://freesound.org/people/unfa/sounds/258342/

`test/arctic_a0001.wav` - http://festvox.org/cmu_arctic/
