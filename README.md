libllsm2
===

Low Level Speech Model (version 2) for high-quality speech analysis/synthesis

About
---

libllsm2 is a C library for analysis, modification and synthesis of digital speech signals. The original libllsm was designed in the context of concatenative synthesis, where a unified parametrization is desired for pitch shifting, concatenation and cross-fading of short segments of speech. Over the years, the project has evolved to accommodate various technical requirements from Synthesizer V and Dreamtonics' experimental projects. The current libllsm2 caters for both concatenative and statistical parametric synthesis.

### Technical Specs

Upon the realization that speech signals can be understood at different levels of abstraction, I supplied the library with multiple analysis/synthesis routines and a flexible data structure. LLSM is a two-layer model of speech. The first layer (layer 0) is basically a harmonic + noise model (HNM); the second layer (layer 1) reinterprets the harmonic parameters in a source-filter setting. From a compositional point of view, layer 0 decomposes speech into periodic and aperiodic parts, and the periodic part is further decomposed into parts related to glottis and vocal tract in layer 1.

Specifically, the layer 0 parametrization consists of the amplitude and phases of the harmonics, power spectral density of the noise, and another harmonic model describing the temporal shape of the noise. The layer 1 parametrization consists of a temporally and spectrally smooth spectral envelope (for an approximated vocal tract transfer function) and parameters for a glottal model.

The analysis procedure goes as the follows. First the fundamental frequency (F0) is estimated using an external library (e.g. libpyin, Nebula). Given the F0 estimation, libllsm2 extracts layer 0 parameters from the speech. Next and optionally, the user may ask libllsm2 to augment the existing representation with layer 1 parameters.

It goes without saying that synthesis is basically to reverse the analysis steps (layer 1 -> layer 0 -> speech). However, libllsm2 can also directly synthesize from layer 1 without going through the harmonic model, and this is much in the same fashion as inverse FFT based vocoders such as WORLD. The direct synthesis pathway, termed as Pulse-by-Pulse (PbP) synthesis, is carefully implemented to give pretty much the same result as a harmonic model, but has some additional advantages when it comes to parameter modification. Real-time synthesis is also supported for both harmonic and PbP synthesis.

Finally, libllsm2 includes a few helper functions to convert frames into a compact representation (fixed dimensional vectors) for statistical parameteric synthesis. By doing so, you can train HMM or DNN models with libllsm2 features, although this means giving up on some advanced features of libllsm2 (mostly regarding phase preservation and temporal noise shaping).

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

Documentation
---

Run `doxygen doxyfile` to generate the API documentation.

Examples of using libllsm2 are given in the tests.

| File | Description |
| --- | --- |
| `test-layer0-anasynth.c` | Layer 0 analysis and synthesis. |
| `test-layer1-anasynth.c` | Layer 0 -> Layer 1 analysis and synthesis. |
| `test-llsmrt.c` | Real-time synthesis. |
| `test-pbpeffects.c` | Pulse-by-Pulse synthesis with growl effect. |
| `test-coder.c` | Speech coding and decoding. |
| `demo-stretch.c` | Time scale modification example. |

### Additional notes on Pulse-by-Pulse synthesis

Since the Januray 2019 update libllsm2 supports Pulse-by-Pulse (PbP) synthesis (in addition to harmonic bank and CZT based methods) of the periodic component. In PbP mode, the glottal flow of each period is generated directly in frequency domain, filtered by the layer 1 vocal tract response, and transformed to time-domain via overlap-add IFFT. This is the same as [WORLD](https://github.com/mmorise/World) vocoder except for the glottal model part.

It is easy to prove that PbP gives the same result as harmonic synthesis given stationary parameters (F0, spectral envelope and glottal parameters), and even in the case of time-varying parameters the difference is hardly perceptible. However, PbP and harmonic synthesis offer a trade-off between speed and degree of control. PbP runs a few times slower than harmonic synthesis but allows period-level parameter modification. One example usage of PbP is to convert normal speech into irregular speech (e.g. growl, see [`test-pbpeffects`](https://github.com/Sleepwalking/libllsm2/blob/master/test/test-pbpeffects.c)).

libllsm2 supports PbP synthesis to the extent that you may switch between PbP mode and harmonic mode on-the-fly during realtime synthesis by attaching `LLSM_FRAME_PBPSYN` flag to frame containers. Period-level parameter modification is possible via `LLSM_FRAME_PBPEFF` and `llsm_pbpeffect` callback structure.

### How does libllsm2 compare to WORLD?

The original libllsm was quite differently purposed. It was meant more for modification than coding and the reason that I had to invent a new vocoder from scratch was that vocoders such as WORLD and STRAIGHT are spectrally near-lossless but temporally lossy. If you want to be temporally lossless, you have to work with phases, which are very nasty to handle. libllsm is near-lossless in both domains, which means that when you synthesize speech without modification, both spectrogram and waveform will look similar to the input.

But, the project has grown considerably alongside Synthesizer V. It has also absorbed ideas from WORLD and VPM. You can now use libllsm2 just like WORLD. Note that libllsm2 uses a slightly different format for feature vectors:

`| VUV | F0 | Rd | Mel-Spectrum | Band Aperiodicity |`

While WORLD is:

`| VUV | F0 | DCT of Mel-Spectrum | Band Aperiodicity |`

Quality-wise, libllsm2 with coded features is not so different from WORLD because they end up giving approximately the same results, although the internal processing is very different. Before coding, libllsm2 is slightly better because it has phase modeling. Again, the difference is not that huge because *we're really after that little bit of room for improvement*.

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

6. M. Morise, F. Yokomori, and K. Ozawa. "WORLD: a vocoder-based high-quality speech synthesis system for real-time applications", IEICE transactions on information and systems, vol. E99-D, no. 7, pp. 1877-1884, 2016.

Test Sounds
---

`test/are-you-ready.wav` - https://freesound.org/people/unfa/sounds/258342/

`test/arctic_a0001.wav` - http://festvox.org/cmu_arctic/
