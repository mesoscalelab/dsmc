#!/bin/bash

mkdir -p bin
g++ -O3 -Wall -I include src/main.cpp -o bin/dsmc
