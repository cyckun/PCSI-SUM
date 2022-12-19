# [PCSI](https://eprint.iacr.org/2021/243.pdf)


### Private Set Operations from Oblivious Switching


### Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [PCSI Source Code](#aby-source-code)
    - [Repository Structure](#repository-structure)
    - [Building](#building-the-aby-framework)
    - [RUN](#testing)

Features
---

This code is provided as a experimental implementation for testing purposes and should not be used in a productive environment. We cannot guarantee security and correctness.

Requirements
---

* Tested  On Ubuntu20.04.
  

PCSI-SUM Source Code
---
#### Repository Structure
* `src/`    - Source code.
@ `test/`   - Test code.
#### Building

1. Clone repository by running:
    ```
   git clone https://github.com/cyckun/PCSI_SUM.git
    ```

2. Enter the Framework directory: `cd PCSI_SUM`

3. Use bazel build with BUILD:
    ```
    make build
    cd build
    cmake ..
    make
    ```

##### Run
```
    ./psi_server
    ./psi_client
```
