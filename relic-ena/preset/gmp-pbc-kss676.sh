#!/bin/bash
cmake -DWSIZE=64 -DTIMER=CYCLE -DFP_PRIME=676 -DARITH=gmp -DFP_QNRES=off $1
