# Near Perfect Numbers Finder

A C++ program designed to identify modified **2-near perfect numbers** within a specified range, satisfying the following condition:

\[
\sigma(n) = 2n \pm d_1 \pm d_2
\]

where \( \sigma(n) \) is the sum of all divisors of \( n \), and \( d_1 \), \( d_2 \) are divisors of \( n \). The program also allows users to focus on numbers of specific forms like \( 2^k*p \). 

## How It Works

The program employs the **Sieve of Eratosthenes** to quickly identify prime numbers up to the specified limit. This sieve helps in verifying whether a number fits the desired forms involving primes.

To calculate the sum of divisors (\( \sigma(n) \)), a sieve-like method is used. For each number up to \( N \), all its divisors are identified by iterating through multiples, allowing for quick computation of their sums.

For narrower searches, the program includes checks to determine if a number fits specific forms:
- **Form 1:** \( n = 2^k \times p \) where \( p \) is an odd prime.
- **Form 2:** \( n = 2^k \times p^2 \) where \( p \) is an odd prime.

Once a candidate number is identified, the program calculates the difference \( \text{diff} = \sigma(n) - 2n \) and checks if this difference can be expressed as \( \pm d_1 \pm d_2 \) using its divisors.

## Usage

1. **Compile the Program:**
   Ensure you have a C++ compiler installed. Compile the program using:
   ```bash
   g++ -o near_perfect_numbers near_perfect_numbers.cpp
   ```
   or similar, depending on your OS and preferred compiler. 

## Results

**Near Perfect Numbers up to 10,000,000 of the form \( n = 2^k \times p^2 \):**

| **n**    | **sigma(n)** | **2n**   | **diff** | **Valid (d1, d2) Combinations**     |
|----------|--------------|----------|----------|--------------------------------------|
| 18       | 39           | 36       | 3        | +1 +2; +2 +1; +6 -3; +9 -6            |
| 36       | 91           | 72       | 19       | +1 +18; +18 +1                        |
| 50       | 93           | 100      | -7       | -2 -5; -5 -2                          |
| 196      | 399          | 392      | 7        | +14 -7                                |
| 200      | 465          | 400      | 65       | +25 +40; +40 +25                      |
| 2312     | 4605         | 4624     | -19      | -17 -2; -2 -17                        |
| 15376    | 30783        | 30752    | 31       | +62 -31                               |
| 1032256  | 2064639      | 2064512  | 127      | +254 -127                             |
| 8454272  | 16908285     | 16908544 | -259     | -2 -257; -257 -2                      |

Near Perfect Numbers up to 1,000,000 of the form \( n = 2^k \times p \):
[Link to CSV](https://raw.githubusercontent.com/0xCUB3/Near-Perfect-Tester/refs/heads/master/one_million_2kp.csv)
