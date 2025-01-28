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
