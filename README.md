= Model Order Reduction in C++

**NOTE :** This work has been included in [Feel++](https://github.com/feelpp/feelpp)

:url-url-def: src/SA/README.md

Install Feel++ : https://docs.feelpp.org/user/latest/install/index.html, with mor and toolboxes dependencies.

Configure the compilation : 

```bash
mkdir build
cd build
cmake ..
```

Build the applications :

* Deterministic sensitivity analysis :
```bash
cd src/DSA
make
```

* [Stochastic sensitivity analysis](src/SA/README.md) :
```bash
cd src/SA
make
```
