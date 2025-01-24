# On the Estimation of Centre of Mass in Periodic Systems

<p align="justify">
Calculation of the centre of mass of a group of particles in a periodically repeating cell is an important aspect of chemical and physical simulation. 
One popular approach, described by Bai and Breen, calculates the centre of mass via the projection of the individual particles' coordinates onto a circle.
This approach, which is mathematically equivalent to finding the first moment of the Fourier series of the mass density, suffers from some numerical error. 
Here we discuss this inaccuracy and propose an extension that overcomes it, enabling improved accuracy across computational simulation. 
</p>

---
<p align="center">
<a href="https://github.com/scams-research/centre-of-mass/actions/workflows/build.yml">
<img src="https://github.com/scams-research/centre-of-mass/actions/workflows/build.yml/badge.svg" alt="Article status"/>
</a>
<a href="https://github.com/scams-research/centre-of-mass/raw/main-pdf/arxiv.tar.gz">
<img src="https://img.shields.io/badge/article-tarball-blue.svg?style=flat" alt="Article tarball"/>
</a>
<a href="https://github.com/scams-research/centre-of-mass/raw/main-pdf/ms.pdf">
<img src="https://img.shields.io/badge/article-pdf-blue.svg?style=flat" alt="Read the article"/>
</a>
<a href="https://arxiv.org/abs/xxxx.xxxxx">
<img src="https://img.shields.io/badge/arXiv-xxxx.xxxxx-orange.svg"/>
</a>
<br><br>
<a href= "https://orcid.org/0009-0002-6808-4573">Harry Richardson</a>,
<a href="https://orcid.org/0000-0003-2659-0806">Josh Dunn</a>,
and 
<a href="https://orcid.org/0000-0003-3381-5911">Andrew R. McCluskey</a>&ast;<br>
&ast;<a href="mailto:andrew.mccluskey@bristol.ac.uk">andrew.mccluskey@bristol.ac.uk</a>
</p>

---

This is the electronic supplementary information (ESI) associated with the publication "On the Estimation of Centre of Mass in Periodic Systems". 
This ESI uses [`showyourwork`](https://show-your.work) to provide a completely reproducible and automated analysis, plotting, and paper generation workflow. 
To run the workflow and generate the paper locally using the cached data run the following: 
```
git clone git@github.com:scams-research/centre-of-mass.git
cd centre-of-mass
pip install showyourwork
showyourwork build 
```
Full details of the workflow can be determined from the [`Snakefile`](https://github.com/scams-research/centre-of-mass/blob/main/Snakefile) and the [`showyourwork.yml`](https://github.com/scams-research/centre-of-mass/blob/main/showyourwork.yml).
