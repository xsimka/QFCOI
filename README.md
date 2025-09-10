# QFCOI
##  [Enhancing Objective Quality Assessment for Compressed Omnidirectional Images with Fusion of Measures](https://ieeexplore.ieee.org/document/11121299)

The repository provides an implementation of the objective OIQA method, Quality Fusion of Compressed Omnidirectional Images (QFCOI), optimized for compressed omnidirectional images by emerging codecs such as JPEG XL, AVIF, and HEIC. 
In addition, the repository offers some validation results, such as for the development and testing of performance evaluation methods. Specifically, these are the outputs of the proposed method QFCOI for the application scenario of OMNIQAD images ([OMNIQAD_results](https://github.com/xsimka/QFCOI/blob/main/Data/OMNIQAD_results)) and QFCOI<sub>H</sub> for the CVIQ dataset ([CVIQ_results](https://github.com/xsimka/QFCOI/blob/main/Data/CVIQ_results)). All details are specified in the article. 

The QFCOI [implementation](https://github.com/xsimka/QFCOI/blob/main/QFCOI.m) is provided. An [example usage](https://github.com/xsimka/QFCOI/tree/main/Example_usage) with a basic script and relevant image data (reference and test image) is also available to help demonstrate the usage. 

<!-- Tento text se nezobrazí v README     ![export_diagramu-1](https://github.com/user-attachments/assets/7fc04381-ac06-4968-9fad-0b86f202acdf)                -->

<p align="center">
  <img src="https://github.com/user-attachments/assets/7fc04381-ac06-4968-9fad-0b86f202acdf" alt="export_diagramu-1" width="500"/>
</p>


##
### Datasets for validation
Proposed method performance was validated on [OMNIQAD](https://zenodo.org/doi/10.5281/zenodo.7607070) and [CVIQ](https://github.com/sunwei925/CVIQDatabase) databases. 
- OMNIQAD (75 distorted images by AVIF, HEIC, and JPEG XL based on 5 reference images)
- CVIQ (176 distorted images by HEVC based on 16 reference images)


## 
If you use QFCOI, the data file, or another part of this research, please cite the article:

```
M. Simka, L. Polak, A. Zizien, and K. Fliegel, "QFCOI-Enhancing Objective Quality Assessment for Compressed Omnidirectional Images With Fusion of Measures," in IEEE Access, vol. 13, pp. 140223-140238, 2025, doi: 10.1109/ACCESS.2025.3597214.
```
###### Citation in LaTeX:
```
@ARTICLE{11121299,
  author={Simka, Marek and Polak, Ladislav and Zizien, Adam and Fliegel, Karel},
  journal={IEEE Access}, 
  title={QFCOI-Enhancing Objective Quality Assessment for Compressed Omnidirectional Images With Fusion of Measures}, 
  year={2025},
  volume={13},
  number={},
  pages={140223-140238},
  doi={10.1109/ACCESS.2025.3597214}}
```

##
Copyright (c) 2025 Marek Simka marek.sim2@gmail.com

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.

##
### Authors of the data
>- Marek Simka (xsimka01@vut.cz) (BUT), Ladislav Polak (BUT), Adam Zizien (CTU), Karel Fliegel (CTU)
>- Brno University of Technology (BUT) and Czech Technical University in Prague (CTU), Czech Republic
