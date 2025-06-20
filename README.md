# sypy

A minimal NumPy-like C++ library. 

## Structure

```
.
├── include
│   └── sypy.h
├── src
│   ├── main.cpp
│   └── sypy.cpp
├── Makefile
└── README.md
```

## Build and Run

```sh
make       
./main     
make clean  
```

## Features

- Templated for generic usage
- Matrix operations:
  - Addition, subtraction, multiplication
  - Scalar multiplication
  - Transpose
  - Determinant
  - Trace
  - Hadamard product
- Array operations:
  - `arange`
  - `linspace`
  - `dot` product
  - `reshape`
- Matrix generation:
  - Zeros
  - Ones
  - Random matrix
  - Unit matrix
- All code is free of comments
- And more... 
