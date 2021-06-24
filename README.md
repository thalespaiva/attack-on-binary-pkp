# attack-on-binary-pkp
This repository contains the implementation of the attack against the Binary variant of the Permuted Kernel Problem proposed by Paiva and Terada at ACNS2021.

# Content
* `Makefile`: Makefile for the C code.
* `attack.py`: SageMath code to estimate attack complexity and to generate challenges for the key search.
* `key_search_optimized.c`: contains the C implementation of the
DFS-based key search algorithm, together with the permutation
extraction after each matching found.
* `sbox32.h`: contains an SBOX that are used to
instantiate a very simple hash function used to compute signatures
of matrices that are invariant under column permutations.


# Initial configuration

## SageMath

To run `attack.py` you need to have [SageMath](https://www.sagemath.org/) installed
in your machine. If you haven't installed it yet, please follow the
[installation guide](https://doc.sagemath.org/html/en/installation/).

The version which I used is given below.
```
$ sage -v
SageMath version 9.3, Release Date: 2021-05-09
```

## C language

To compile and run the C program `key_search_optimized.c`, you will need the
[M4RI library](https://bitbucket.org/malb/m4ri/src/master/).

In Arch Linux, it can be found in th[\*]
[\*]
[\*]
[\*]e `community` repo, and you can install it
using `$ pacman -Ss m4ri`. Notice that if you install M4RI manually from the source,
you will have to change the Makefile to tell `gcc` where it should look for the include
files and library.

# Compiling the C code

To compile `key_search_optimized.c`, just run `$ make` and it will create an executable
called `kso`.
```
$ make
gcc -DUSE_LARGE_NUMBER_OF_LOW_WEIGHT_CODEWORDS=1 -Wall -Werror -pedantic -std=c99 -O3   -fopenmp -L/usr/local/lib -lm4ri -o kso key_search_optimized.c
$ ./kso
Usage: ./kso <path to directory containing challenge files>
```

# Estimating the attack complexity

Let us first run the SageMath script with the `--help` flag.
If you have successfully installed SageMath, this is what you should see.

```
$ sage attack.py --help
usage: attack.py [-h] [--directory DIRECTORY]
                 {generate_challenge,estimate_complexity} m n l attack_param_w
                 attack_param_la

positional arguments:
  {generate_challenge,estimate_complexity}
  m
  n
  l
  attack_param_w
  attack_param_la

optional arguments:
  -h, --help            show this help message and exit
  --directory DIRECTORY
                        Path to the directory where a challenge will be
                        generated.

```

Let us estimate the attack complexity for the following parameters.
* `m = 15`
* `n = 38`
* `l = 10`
* `attack_param_w = 8`
* `attack_param_la = 9`

Then we run
```
$ time sage attack.py estimate_complexity 15 38 10 8 9
{'av_size_of_LWSK': 47766.3164062500,
 'fraction_of_keys': 0.13569077665209361,
 'nchildren_per_level_analytic': [47766.3164062500,
                                  17610.1727987815,
                                  4317.59234071975,
                                  744.465647853847,
                                  55.5385243709078,
                                  5.33747895845528,
                                  0.750384907715224,
                                  0.142191368819397,
                                  0.0364963129288272],
 'nchildren_per_level_simulation': [47766.3164062500,
                                    15907.6847828934,
                                    2040.58370000545,
                                    272.582609277166,
                                    33.8266860325894,
                                    4.15214952427055,
                                    0.738420381168792,
                                    0.140651501175008,
                                    0.0390698614375023],
 'work_factor_analytic': 6.55347028466972e19,
 'work_factor_of_search_analytic': 8.01497511881623e17,
 'work_factor_of_search_simulation': 5.93627104958385e16,
 'work_factor_permutations_analytic': 81.7653228802241,
 'work_factor_with_search_simulation': 4.85381119073750e18}
sage attack.py estimate_complexity 15 38 10 8 9  5.46s user 0.14s system 99% cpu 5.617 total
```

Notice that it returns a detailed response. You can find the following information
* `av_size_of_LWSK`: The expected size of set `LWSK` of vectors of weight `attack_param_w` in the kernel of `Vpub`
* `fraction_of_keys`: The fraction of keys against which the provided attack parameters are effective.
* `nchildren_per_level_analytic`: The number of children nodes in each level (analytical estimation).
* `nchildren_per_level_simulation`: The number of children nodes in each level (estimated with simulations).
* `work_factor_analytic`: The total work factor computed analytically.
* `work_factor_of_search_analytic`: The work factor of the search computed analytically.
* `work_factor_of_search_simulation`: The work factor of the search computed with simulations.
* `work_factor_permutations_analytic`: The work factor of the permutation extraction computed with analytically.
* `work_factor_with_search_simulation`: **Important.** Total work factor computed analytically, except for the search,
for which simulations are used. In general, we use this to generate graphs because it is too slow to simulate the
permutation extraction, and the analytically computed values tend to overestimate the work factor in practice.

# Generating a challenge

To generate a challenge input for the key search in C,
we do the following.
```
$ time sage attack.py generate_challenge 15 38 10 4 14 --directory test_challenge
Running Stern`s algorithm with parameters p=1, l=1, error=1e-06
[100% complete] # Processed codewords: (24 unique, 335 total)
Running Stern`s algorithm with parameters p=1, l=2, error=1e-06
[100% complete] # Processed codewords: (108 unique, 2271 total)
sage attack.py generate_challenge 15 38 10 4 14 --directory test_challenge  10.18s user 0.16s system 99% cpu 10.363 total
```

This will create a random Binary PKP key **that can be attacked with parameters `attack_param_w = 4`
and `attack_param_la = 14`**. The challenge, essentially consisting of matrices `A`, `Vpub`,
and sets `LWSA` and `LWSK`, which are computed using Stern's algorithm, are stored in `test_challenge`
directory.
```
$ ls test_challenge
A_low_weight_head.jcf  la.param  LWSA.sparse_matrix  LWSK.sparse_matrix  secret_subset.indexes  Vpub.jcf
```


# Running the key search

Simply run `./kso test_challenge`. Notice that `kso` will generate a lot of debugging information.
We will reduce the verbosity of this code in the future.
```
$ time ./kso test_challenge
reading 15 x 38 matrix with at most 81 non-zero entries (density at most: 0.14211)
reading 38 x 10 matrix with at most 181 non-zero entries (density at most: 0.47632)
A_low_weight_head
[1   :   1:   1:1   :    :    :    :    :    :  ]
[    :    :  1 :  1 :    :    :    : 1  :    : 1]
[    :    :  1 :1   : 1  :    :    :    :1   :  ]
[ 1  :    : 11 : 1  :    :    :    :    :    :  ]
[    : 11 :  1 :    :    :    :    : 1  :    :  ]
[    :1   :   1:    :    :    :   1:   1:    :  ]
[    :    : 1  : 1 1:    :    :   1:    :    :  ]
[    :   1:    :1   :    :    :   1:    :1   :  ]
[1 1 :    :    :    :    :    : 1  :    :  1 :  ]
[1 1 :1   :    :    :    :    :    :    : 1  :  ]
[1   : 11 :    :    :    :1   :    :    :    :  ]
[    :   1:    :    :1   :  1 : 1  :    :    :  ]
[    :    :    :    :    :    :  1 :  11:  1 :  ]
[    :    :1   :   1:  1 :    :    :    :   1:  ]
[111 :1 1 :1 11:111 :11  :111 :  11:1111:  1 :11]
A_LWSA
[1   :   1:   1:1   :    :    :    :    :    :  ]
[    :    :  1 :  1 :    :    :    : 1  :    : 1]
[    :    :  1 :1   : 1  :    :    :    :1   :  ]
[ 1  :    : 11 : 1  :    :    :    :    :    :  ]
[    : 11 :  1 :    :    :    :    : 1  :    :  ]
[    :1   :   1:    :    :    :   1:   1:    :  ]
[    :    : 1  : 1 1:    :    :   1:    :    :  ]
[    :   1:    :1   :    :    :   1:    :1   :  ]
[1 1 :    :    :    :    :    : 1  :    :  1 :  ]
[1 1 :1   :    :    :    :    :    :    : 1  :  ]
[1   : 11 :    :    :    :1   :    :    :    :  ]
[    :   1:    :    :1   :  1 : 1  :    :    :  ]
[    :    :    :    :    :    :  1 :  11:  1 :  ]
[    :    :1   :   1:  1 :    :    :    :   1:  ]
Vpub
[1111:    :  ]
[ 1  :1   :  ]
[  1 : 1  :  ]
[1  1: 1  :  ]
[  11:  1 :  ]
[1   :1 1 :  ]
[    :1  1:  ]
[1   :1  1:  ]
[11 1: 1 1:  ]
[    :    :1 ]
[ 1 1:    :1 ]
[1 1 :  1 :1 ]
[  1 :1 1 :1 ]
[    : 11 :1 ]
[11  :111 :1 ]
[1 11: 1 1:1 ]
[ 111:11 1:1 ]
[   1:    : 1]
[1111:1   : 1]
[1   :11  : 1]
[111 :  1 : 1]
[   1:1 1 : 1]
[111 :   1: 1]
[   1:11 1: 1]
[11 1:  11: 1]
[1  1:1 11: 1]
[    : 111: 1]
[11 1:1111: 1]
[ 11 : 1  :11]
[ 111:  1 :11]
[1 1 : 11 :11]
[  1 :   1:11]
[1 1 : 1 1:11]
[1 11:11 1:11]
[1 1 :1 11:11]
[  11:1 11:11]
[ 11 : 111:11]
[1 11:1111:11]


Each permutation extraction will take:
(7!/1!)(2!/1!)(2!/1!)(2!/0!)(2!/1!)(2!/1!)(3!/1!)
npossibles[0] = 108
npossibles[1] = 56
npossibles[2] = 12
npossibles[3] = 4
npossibles[4] = 2
npossibles[5] = 9
npossibles[6] = 1
npossibles[7] = 1
npossibles[8] = 1
npossibles[9] = 1
npossibles[10] = 1
npossibles[11] = 1
npossibles[12] = 1
npossibles[13] = 1
[\*] Information on current level
    Worst case total time: 49.912276 seconds

    Current path  21  33  73  43  90  84
[\*] Information on current level
    Worst case total time: 52.021327 seconds

    Current path  39  27  67  90  66  84
[\*] Information on current level
    Worst case total time: 51.968129 seconds

    Current path  63  18  46  21  26 101  13
[\*] Information on current level
    Worst case total time: 52.855516 seconds

    Current path  76  51  87  31  39  70  63
 92  83   5  84  25  23  60 102 100  61 105  45  10   6
LEAF!
[\*] Information on current level
    Worst case total time: 53.393522 seconds

    Current path  94  57  92  72  98  59
    ... :^(
 92  25   5  84  83  23  60 102 100  61   7  45  10   6
LEAF!
---- PARTIALS --------
[    :    :1 ]
[    : 11 :1 ]
[1   :11  : 1]
[ 111:  1 :11]
[11 1:  11: 1]
[   1:    : 1]
[   1:11 1: 1]
[  11:1 11:11]
[1 11: 1 1:1 ]
[1 1 : 11 :11]
[11  :111 :1 ]
[1111:1   : 1]
[111 :  1 : 1]
[1 1 : 1 1:11]
[ 1 1:    :1 ]
[ 1  :1   :  ]
[  11:  1 :  ]
[  1 : 1  :  ]
[111 :   1: 1]
[ 11 : 111:11]
[11 1:1111: 1]
[1 1 :1 11:11]
[    :1  1:  ]
-------------
[1 11:11 1:11]
[    : 111: 1]
[11 1: 1 1:  ]
[1  1: 1  :  ]
[  1 :   1:11]
[1   :1  1:  ]
[   1:1 1 : 1]
[1   :1 1 :  ]
[ 11 : 1  :11]
[1  1:1 11: 1]
[1111:    :  ]
[1 11:1111:11]
[  1 :1 1 :1 ]
[1 1 :  1 :1 ]
[ 111:11 1:1 ]
---- END PARTIALS --------
================== Vsec = KEY ===========
[1 1 : 1 1:11]
[1  1: 1  :  ]
[1111:    :  ]
[    :    :1 ]
[11  :111 :1 ]
[1111:1   : 1]
[111 :  1 : 1]
[ 1  :1   :  ]
[ 111:11 1:1 ]
[   1:1 1 : 1]
[  11:1 11:11]
[1   :1  1:  ]
[ 11 : 1  :11]
[1 11: 1 1:1 ]
[    : 111: 1]
[    :1  1:  ]
[  1 :1 1 :1 ]
[11 1: 1 1:  ]
[11 1:1111: 1]
[    : 11 :1 ]
[1 11:1111:11]
[1 11:11 1:11]
[ 1 1:    :1 ]
[1   :11  : 1]
[ 111:  1 :11]
[  11:  1 :  ]
[1 1 :  1 :1 ]
[1 1 : 11 :11]
[11 1:  11: 1]
[  1 :   1:11]
[  1 : 1  :  ]
[111 :   1: 1]
[1   :1 1 :  ]
[1  1:1 11: 1]
[ 11 : 111:11]
[1 1 :1 11:11]
[   1:    : 1]
[   1:11 1: 1]
============== KEY FOUND ====================
```

# License
    MIT


# Author
    Thales Paiva