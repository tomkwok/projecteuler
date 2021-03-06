# Project Euler problems solution in C

[![Build status][img_build_status]][build_status] [![Language grade: C/C++][img_lgtm]][lgtm] [![Technical debt][img_sonar]][sonar]

[img_build_status]: https://github.com/tomkwok/projecteuler/workflows/projecteuler/badge.svg
[build_status]: https://github.com/tomkwok/projecteuler/actions?query=branch%3Amaster

[img_lgtm]: https://img.shields.io/lgtm/grade/cpp/g/tomkwok/projecteuler.svg?logo=lgtm&logoWidth=18
[lgtm]: https://lgtm.com/projects/g/tomkwok/projecteuler/latest/files/

[img_sonar]: https://img.shields.io/sonar/tech_debt/tomkwok_projecteuler?logo=sonarsource&server=https%3A%2F%2Fsonarcloud.io
[sonar]: https://sonarcloud.io/dashboard?id=tomkwok_projecteuler

This repository contains my solutions to [Project Euler problem](https://projecteuler.net/archives) 1 to 10 written in C. There is some use of C preprocessor macros in code.

The C code was written by November 2014 at the latest. It is mostly intact except for a few fixes (i) for compilation errors and typos in comments; and (ii) to change camel case function names to underscores delimited name phrases.

![Scan of print out of projecteuler.c in 2014](projecteuler.c.jpg)

For each problem, at least one function was written to return the solution to a generalized version of the problem. The parameters given in problem statements were put in function calls in `main()` rather than being hard-coded within the program logic.

For some of the problems solved, a faster or recursive implementation was provided as an alternative to the more straight-forward or naive implementation.

While it is a good exercise to learn C (or any programming language) when trying to solve problems in Project Euler, a cleverer solution is usually derived by doing some amount of mathematics. The C code was annotated with extensive inline comments to explain the methodology of solution with mathematical formulas.

## A list of Project Euler problems solved

The following list contains my one-line re-interpretation of the problems.

1. sum of multiples of a or b inclusive
2. sum of even Fibonacci numbers
3. largest prime factor (prime factorization)
4. largest (even-digit) palindrome product of two n-digit numbers
5. LCM of 1 to n
6. square of sum minus sum of squares
7. the nth prime: given number of primes, find the value of prime
8. n adjacent digits in series having the greatest product
9. Pythagorean triplet given perimeter
10. summation of primes: given (max value of primes + 1), find the sum of values
11. (todo) largest product of adjacent numbers in a 2-dimensional array
