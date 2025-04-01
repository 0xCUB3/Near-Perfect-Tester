#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // Needed for std::max, std::sort, std::unique, binary_search
#include <limits>
#include <iomanip>
#include <fstream>
#include <cctype>    // For tolower
#include <cstdint>   // For uint64_t
#include <set>       // For unique divisor finding and number generation
#include <stdexcept> // For stoi exceptions
#include <sstream>   // For stringstream in factorization

using namespace std;

// Helper function to safely multiply two uint64_t, returning 0 on overflow
// Using __int128 for intermediate check
uint64_t safe_multiply(uint64_t a, uint64_t b) {
    // Early exit for 0 to avoid unnecessary __int128 calculation
    if (a == 0 || b == 0) return 0;
    // Check if multiplication would obviously overflow based on leading bits
    // Simple check: if a > max/b roughly
    if (b > 0 && a > numeric_limits<uint64_t>::max() / b) return 0; // Check b > 0 to avoid division by zero

    // Use 128-bit integer if available for a precise check
    #ifdef __SIZEOF_INT128__ // Check if the compiler supports 128-bit integers
    unsigned __int128 result = (unsigned __int128)a * b;
    if (result > numeric_limits<uint64_t>::max()) {
        return 0; // Indicate overflow
    }
    return (uint64_t)result;
    #else
    // Fallback if 128-bit integers aren't available (less precise)
    // Already checked a > max/b above, which is a decent heuristic
    return a * b;
    #endif
}

// Helper function to safely raise uint64_t to a power, returning 0 on overflow
uint64_t safe_power(uint64_t base, int exp) {
    uint64_t res = 1;
    if (exp < 0) return 0; // Negative exponent not handled here
    if (base == 0) return (exp == 0) ? 1 : 0; // 0^0=1, 0^positive=0
    if (exp == 0) return 1; // x^0 = 1 for x!=0

    uint64_t current_base = base;
    int current_exp = exp;

    // Exponentiation by squaring
    while (current_exp > 0) {
        if (current_exp % 2 == 1) { // If exponent is odd
            res = safe_multiply(res, current_base);
            if (res == 0 && (current_base != 0 || exp==1) ) return 0; // Check overflow, handle base^1 case
        }
        current_exp /= 2;
        if (current_exp > 0) { // Only need to square if exponent isn't done
             // Check before squaring: if base > sqrt(max), squaring will overflow
             uint64_t sqrt_max = 1ULL << 32; // Approximation of sqrt(ULLONG_MAX)
             if (current_base > sqrt_max) {
                // Base is large, squaring might overflow, use safe_multiply carefully
                 uint64_t next_base_sq = safe_multiply(current_base, current_base);
                 if (next_base_sq == 0 && current_base != 0) return 0; // Overflow detected by safe_multiply
                 current_base = next_base_sq;
             } else {
                 // Base is smaller, direct multiplication is likely safe but use safe_multiply anyway
                 uint64_t next_base_sq = safe_multiply(current_base, current_base);
                  if (next_base_sq == 0 && current_base != 0) return 0; // Overflow detected by safe_multiply
                 current_base = next_base_sq;
             }
        }
    }
    return res;
}


/**
 * Function: simpleSieve
 * ---------------------
 * Finds all prime numbers up to limit using the simple sieve of Eratosthenes.
 */
vector<uint64_t> simpleSieve(uint64_t limit) {
    if (limit < 2) return {};
    vector<bool> isPrime(limit + 1, true);
    isPrime[0] = isPrime[1] = false;

    for (uint64_t i = 2; i * i <= limit; ++i) {
        if (isPrime[i]) {
            for (uint64_t j = i * i; j <= limit; j += i) {
                isPrime[j] = false;
            }
        }
    }

    vector<uint64_t> primes;
    // Reserve approximate space, avoid frequent reallocations
    if (limit > 10) { // log(1) is problematic, avoid for small limit
        try {
             primes.reserve(static_cast<size_t>(1.2 * limit / log(limit))); // Approximation + buffer
        } catch (const std::bad_alloc& e) {
             cerr << "Memory allocation failed for primes vector reservation." << endl;
             // Continue without reservation
        }
    } else {
        primes.reserve(5); // Small number
    }

    for (uint64_t i = 2; i <= limit; ++i) {
        if (isPrime[i]) {
            primes.push_back(i);
        }
    }
    return primes;
}

/**
 * Function: isPrime
 * -----------------
 * Checks if a number is prime using trial division with the provided list of primes.
 */
bool isPrime(uint64_t n, const vector<uint64_t>& primes) {
    if (n < 2) return false;
    if (n == 2 || n == 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;

    uint64_t sqrt_n = static_cast<uint64_t>(sqrt(static_cast<long double>(n)));

    // Check against provided primes
    for (const auto& p : primes) {
        if (p > sqrt_n) break;
        // Skip 2 and 3 as they were already checked
        if (p > 3 && n % p == 0) return false;
    }

    // If sqrt_n is larger than the largest prime checked, continue checking
    // Optimized trial division starting from 5, stepping by 6 (5, 7, 11, 13, ...)
    // Start from the largest prime checked + step or 5 if no primes were checked
    uint64_t start_check = 5;
    if (!primes.empty() && primes.back() > 3) {
        // Start checking from the next potential prime after the largest one we used
        start_check = (primes.back() / 6) * 6 + 5;
        if (start_check <= primes.back()) {
             start_check += 6; // Ensure we start after the last checked prime
        }
    }

    for (uint64_t i = start_check; i <= sqrt_n; i += 6) {
         if (n % i == 0 || n % (i + 2) == 0) {
             return false;
         }
     }

    return true;
}


/**
 * Function: computeSigma
 * ----------------------
 * Computes the sum of divisors (sigma function) of a number using its prime factorization.
 * sigma(p^a) = (p^(a+1) - 1) / (p - 1)
 */
uint64_t computeSigma(uint64_t n, const vector<uint64_t>& primes) {
    if (n == 0) return 0;
    if (n == 1) return 1;

    uint64_t sigma = 1;
    uint64_t temp_n = n; // Use a temporary variable for modification

    // ADDED '-> uint64_t' here
    auto process_prime_factor = [&](uint64_t p) -> uint64_t {
        if (temp_n % p == 0) {
            uint64_t p_power_a_plus_1 = p; // Start with p^1
            int a = 0;
            while (temp_n % p == 0) {
                a++;
                // Prevent infinite loop if temp_n is not reduced
                uint64_t next_temp_n = temp_n / p;
                if (next_temp_n >= temp_n && p > 1) { // Division didn't reduce value (e.g., temp_n was 0 or small)
                    cerr << "\nError: Division by p=" << p << " did not reduce temp_n=" << temp_n << " in computeSigma. Aborting factor." << endl;
                    return (uint64_t)0; // Signal error/overflow state
                }
                temp_n = next_temp_n;

                // Need to calculate p^(a+1) for the formula safely
                uint64_t next_power = safe_multiply(p_power_a_plus_1, p);
                 if (next_power == 0 && p!=0 && a > 0) { // Check a > 0 because p^1 might be 0 if p is 0 (though p should be prime)
                    // Overflow calculating p^(a+1)
                    return (uint64_t)0; // Signal overflow
                }
                 // Handle case where p=1 somehow (shouldn't happen for primes) or base was 0 initially
                 if (p_power_a_plus_1 == 0 && a == 0 && p > 0) {
                    p_power_a_plus_1 = p; // Re-initialize if first power calc was issue
                 } else {
                    p_power_a_plus_1 = next_power;
                 }

            }

            // Calculate term = (p^(a+1) - 1) / (p - 1)
             if (p_power_a_plus_1 == 0 && a > 0) { // Check if power calculation overflowed on last step
                 return (uint64_t)0; // Signal overflow
             }
             uint64_t term_numerator = p_power_a_plus_1 - 1;
             uint64_t term_denominator = p - 1;

             // Check for division by zero (p=1 case, shouldn't happen for primes)
             if (term_denominator == 0) {
                 if (p == 1) return sigma; // Sigma doesn't change for factor 1
                 return (uint64_t)0; // Error for p=1 if it wasn't expected
             }

             // Check if division is exact (it must be for sigma formula)
             // Note: If p_power_a_plus_1 overflowed, numerator might wrap around, making division non-exact
             if (term_numerator == numeric_limits<uint64_t>::max() && p_power_a_plus_1==0) {
                 // This means p^(a+1) wrapped around to 0. Indicates definite overflow.
                 return (uint64_t)0; // Signal overflow
             }
             if (term_numerator % term_denominator != 0) {
                 // This likely indicates an issue, possibly overflow.
                 // Fallback: Calculate sum 1 + p + p^2 + ... + p^a iteratively
                 // This fallback might *also* overflow but is worth trying if formula failed
                 cerr << "\nWarning: Sigma formula division failed for p=" << p << ", a=" << a << ", n=" << n <<". Trying iterative sum." << endl;
                 uint64_t term = 1;
                 uint64_t current_p_power = 1;
                 for (int i=0; i<a; ++i) {
                     current_p_power = safe_multiply(current_p_power, p);
                     if (current_p_power == 0) return (uint64_t)0; // Overflow in power

                     // Safe addition check
                     #ifdef __SIZEOF_INT128__
                     unsigned __int128 next_term = (unsigned __int128)term + current_p_power;
                     if (next_term > numeric_limits<uint64_t>::max()) return (uint64_t)0; // <<< CORRECTED: return uint64_t
                     term = (uint64_t)next_term;
                     #else
                     // Fallback check (less precise but necessary without 128-bit)
                     if (term > numeric_limits<uint64_t>::max() - current_p_power) return (uint64_t)0; // <<< CORRECTED: return uint64_t
                     term += current_p_power;
                     #endif
                 }
                 // Multiply the main sigma by the calculated term safely
                 uint64_t next_sigma = safe_multiply(sigma, term);
                  if (next_sigma == 0 && sigma != 0 && term != 0) return (uint64_t)0; // Overflow
                  sigma = next_sigma;

             } else {
                 // Division was exact, proceed with formula result
                 uint64_t term = term_numerator / term_denominator;
                 uint64_t next_sigma = safe_multiply(sigma, term);
                 if (next_sigma == 0 && sigma != 0 && term != 0) return (uint64_t)0; // Overflow
                 sigma = next_sigma;
             }
             // Check if the overall sigma calculation resulted in 0 unexpectedly
             if (sigma == 0 && n!=0) return (uint64_t)0; // Signal overflow in sigma product
        }
        return (uint64_t)1; // Signal success (or prime not a factor), using 1 to indicate no overflow occurred here
    };

    // Handle factor 2
    if (process_prime_factor(2) == 0) return 0; // Check if the lambda signaled overflow

    // Handle odd prime factors
    uint64_t sqrt_temp_n_orig = (temp_n > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n))) : 0;
    for (const auto& p : primes) {
        if (p == 2) continue; // Already handled
        if (temp_n == 1) break; // Fully factored
        // Optimization: Check primes up to sqrt of *original* remaining n after factoring 2
        if (p > sqrt_temp_n_orig && p*p > temp_n) break; // Combine checks: if p > original sqrt AND p*p > current n, we can stop.

        if (temp_n % p == 0) {
             if (process_prime_factor(p) == 0) return 0; // Process factor and check overflow
        }
    }

    // If temp_n is still greater than 1, it must be a prime factor itself
    if (temp_n > 1) {
        // Calculate term = 1 + temp_n safely (sigma(p^1))
        uint64_t term;
        if (temp_n == numeric_limits<uint64_t>::max()) return 0; // 1 + max overflows
        term = temp_n + 1; // Safe since temp_n < max

        uint64_t next_sigma = safe_multiply(sigma, term);
        if (next_sigma == 0 && sigma != 0 && term != 0) return 0; // Overflow in sigma product
        sigma = next_sigma;
    }

    return sigma;
}


/**
 * Function: computeDivisors
 * -------------------------
 * Computes all distinct divisors of a number n.
 */
vector<uint64_t> computeDivisors(uint64_t n, const vector<uint64_t>& primes) {
    set<uint64_t> divisors_set;
    if (n == 0) return {};
    divisors_set.insert(1);
    if (n == 1) return {1};

    uint64_t temp_n = n;
    bool overflow_detected = false; // Flag to track overflow

    // Function to add multiples of a prime power to the divisor set
    auto add_prime_power_divisors = [&](uint64_t p) {
        vector<uint64_t> current_divisors(divisors_set.begin(), divisors_set.end());
        uint64_t power = 1;
        while (temp_n % p == 0) {
            power = safe_multiply(power, p);
            if (power == 0) { // Check for overflow
                cerr << "\nWarning: Overflow computing power for divisor generation at n=" << n << ", p=" << p << ". Results might be incomplete." << endl;
                overflow_detected = true;
                return; // Stop processing this prime
            }
            for (uint64_t d : current_divisors) {
                 uint64_t new_divisor = safe_multiply(d, power);
                 if (new_divisor == 0 && d != 0 && power != 0) { // Check overflow, avoid 0*x or x*0 case false positive
                     cerr << "\nWarning: Overflow computing new divisor (" << d << "*" << power << ") for n=" << n << ". Results might be incomplete." << endl;
                     overflow_detected = true;
                     // Decide whether to continue adding other divisors for this power, or stop entirely
                     // Let's stop processing this prime factor entirely on overflow for safety.
                     return;
                 }
                 divisors_set.insert(new_divisor);
            }
            temp_n /= p;
            if (temp_n == 0 && n!=0) { // Should not happen if p is a factor, indicates problem
                cerr << "\nError: temp_n became 0 unexpectedly during divisor generation for n=" << n << endl;
                overflow_detected = true; // Mark as overflow/error state
                return;
            }
        }
    };

    // Handle factor 2
    add_prime_power_divisors(2);
    if (overflow_detected) return {}; // Return empty on overflow

    // Handle odd prime factors
    uint64_t sqrt_temp_n_orig = (temp_n > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n))) : 0;
    for (const auto& p : primes) {
        if (p == 2) continue;
        if (temp_n == 1) break; // Optimization: if temp_n reduced to 1, no more factors
        if (p > sqrt_temp_n_orig) break; // Optimization

        uint64_t p_squared = safe_multiply(p, p);
         if (p_squared == 0 && p > 1) { /* Overflow, rely on p > sqrt_temp_n */ }
         else if (p_squared > temp_n && temp_n > 1) {
             // If p*p > current temp_n, p cannot be a factor unless p = temp_n.
             // The temp_n > 1 check after loop handles the p = temp_n case.
             // We can break the loop based on the original sqrt check.
         }

        if (temp_n % p == 0) {
            add_prime_power_divisors(p);
             if (overflow_detected) return {}; // Return empty on overflow
             // Recompute sqrt? Not strictly necessary due to sqrt_temp_n_orig check
             // sqrt_temp_n = (temp_n > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n))) : 0;
        }
    }

    // If temp_n is still greater than 1, it's a prime factor
    if (temp_n > 1 && !overflow_detected) {
        add_prime_power_divisors(temp_n);
        if (overflow_detected) return {}; // Return empty on overflow
    }

    // Final check before returning
    if (overflow_detected) {
         return {};
    }

    return vector<uint64_t>(divisors_set.begin(), divisors_set.end());
}


/**
 * Enum: SignOption
 */
enum SignOption {
    ALL = 1,
    MIXED_SIGNS,
    ONLY_POSITIVE,
    ONLY_NEGATIVE
};

/**
 * Function: findValidCombinations
 * -------------------------------
 * Finds all valid combinations of DISTINCT divisors d1 and d2 with their respective signs
 * that satisfy the equation: sigma(n) = 2n ± d1 ± d2
 *
 * Parameters:
 *  - D: Reference to the vector containing all distinct divisors of n (must include 1).
 *  - diff: The difference sigma(n) - 2n.
 *  - signOption: The user's choice for sign combinations.
 *
 * Returns:
 *  - A vector of strings, each representing a valid combination of distinct d1 and d2 with signs.
 */
vector<string> findValidCombinations(const vector<uint64_t> &D, long long diff, SignOption signOption) {
    vector<string> validPairs;
    set<string> uniquePairs; // Use set to automatically handle duplicates and sorting

    // Iterate through all possible pairs (d1, d2)
    for (uint64_t d1 : D) {
        // Skip d1=0 if it somehow exists (shouldn't for n>0)
        if (d1 == 0) continue;

        long long ld1 = (d1 <= (uint64_t)numeric_limits<long long>::max()) ? static_cast<long long>(d1) : -1;
        if (ld1 == -1) continue; // Skip if d1 divisor exceeds long long max

        for (uint64_t d2 : D) {
            // Skip d2=0 if it somehow exists
            if (d2 == 0) continue;

            // <<<--- ADD THIS CHECK ---<<<
            if (d1 == d2) {
                continue; // Ensure divisors are distinct
            }
            // >>>---------------------->>>

            long long ld2 = (d2 <= (uint64_t)numeric_limits<long long>::max()) ? static_cast<long long>(d2) : -1;
            if (ld2 == -1) continue; // Skip if d2 divisor exceeds long long max


            // Check combinations based on signOption
            if (signOption == ALL || signOption == ONLY_POSITIVE) {
                // ++
                 if (ld1 <= numeric_limits<long long>::max() - ld2) { // Check for potential overflow before addition
                     if ((ld1 + ld2) == diff) {
                         // Ensure consistent ordering in the set (e.g., smaller number first) to avoid "+a +b" and "+b +a" duplicates if needed?
                         // The set handles string duplicates "as is". If "+1 +2" and "+2 +1" should be treated as one, more logic needed.
                         // For now, keep them distinct as the order might matter based on d1/d2 definition.
                         uniquePairs.insert("+" + to_string(d1) + " +" + to_string(d2));
                     }
                 }
            }
            if (signOption == ALL || signOption == MIXED_SIGNS) {
                // +-
                if ((ld1 - ld2) == diff) {
                    uniquePairs.insert("+" + to_string(d1) + " -" + to_string(d2));
                }
                // -+
                if ((ld2 - ld1) == diff) { // Equivalent to (-ld1 + ld2)
                    uniquePairs.insert("-" + to_string(d1) + " +" + to_string(d2));
                }
            }
            if (signOption == ALL || signOption == ONLY_NEGATIVE) {
                // --
                 if (ld1 <= numeric_limits<long long>::max() - ld2) {
                     // Check underflow: (-ld1 - ld2) == diff
                     // Use the non-overflowing positive sum ld1+ld2 for comparison
                     long long neg_sum = - (ld1 + ld2); // Calculate the negative sum safely
                      if (neg_sum == diff) {
                           // Consistent ordering for "--" pairs? Similar consideration as "++".
                           uniquePairs.insert("-" + to_string(d1) + " -" + to_string(d2));
                     }
                 } else {
                     // ld1 + ld2 overflowed positive range, so -ld1 - ld2 definitely underflows
                      // No need to check against diff unless diff could be LLONG_MIN etc.
                      // Assume no match on overflow/underflow for simplicity.
                 }
            }
        }
    }

    validPairs.assign(uniquePairs.begin(), uniquePairs.end());
    return validPairs;
}

/**
 * Function: generateNumbers
 * -------------------------
 * Generates a list of numbers to process based on the selected form.
 */
vector<uint64_t> generateNumbers(uint64_t N, int choice, const vector<uint64_t> &primes, int exponent_m_choice3 = 1) {
    set<uint64_t> numbers_set; // Use set to avoid duplicates automatically

    if (choice == 1) {
        cout << "\nWarning: Generating all numbers up to N=" << N << " can be very slow and memory intensive." << endl;
        cout << "Consider using a smaller N or choosing a specific form (2, 3, or 4)." << endl;
        // All numbers from 1 to N
        for (uint64_t n = 1; n <= N; ++n) {
            numbers_set.insert(n);
            if (n == numeric_limits<uint64_t>::max()) break; // Prevent overflow in loop
        }
    } else if (choice == 2) {
        // Numbers of the form n = 2^k * p
        for (size_t i = 0; i < primes.size(); ++i) {
            uint64_t p = primes[i];
            if (p < 3) continue; // Skip p=2

            if (p > N) break; // Optimization: if p > N, no 2^k * p can be <= N (for k>=0)

            uint64_t current_n = p; // Start with k=0 (n=p)
            if (current_n <= N) {
                 numbers_set.insert(current_n);
            } else {
                 continue; // p itself > N
            }


            uint64_t k_power_of_2 = 1;
            while(true) {
                 // Prepare for next power of 2 (k >= 1)
                 uint64_t next_k_power = safe_multiply(k_power_of_2, 2);
                 if (next_k_power == 0 || next_k_power < k_power_of_2 ) { // Overflow check 2^k
                     break;
                 }
                 k_power_of_2 = next_k_power; // Now 2^k for k>=1

                // Calculate n = p * 2^k
                 current_n = safe_multiply(p, k_power_of_2);

                 if (current_n == 0 || current_n > N) {
                     break; // Overflow or exceeded limit N. Stop for this p.
                 }
                 numbers_set.insert(current_n);
            }
        }
    } else if (choice == 3) {
        // Numbers of the form n = 2^k * p^m (using exponent_m_choice3)
         for (size_t i = 0; i < primes.size(); ++i) {
            uint64_t p = primes[i];
            if (p < 3) continue; // Skip p=2

            uint64_t p_power_m = safe_power(p, exponent_m_choice3);
            // If p^m overflows or is > N, stop for this p and subsequent larger primes
            if (p_power_m == 0 || p_power_m > N) {
                 // Optimization check: If p > N^(1/m), then p^m will exceed N.
                 // Can we break the outer loop early?
                 // If p > N, we can break.
                 if (p > N) break;
                 // If p^m > N, maybe a smaller p with a larger power of 2 works?
                 // But if p^m > N, then p^m * 2^k > N for all k>=0. So continue is correct.
                 continue;
            }

            // Process n = p^m * 2^k
            uint64_t base = p_power_m;

            // Add case k=0
            if (base <= N) {
                numbers_set.insert(base);
            } else {
                // base (p^m) > N, no need to continue with powers of 2
                continue; // Already handled by the check before this block, but safe
            }

            uint64_t k_power_of_2 = 1; // Represents 2^k
            while(true) {
                 // Prepare for next power of 2 (k >= 1)
                  uint64_t next_k_power = safe_multiply(k_power_of_2, 2);
                  if (next_k_power == 0 || next_k_power < k_power_of_2 ) { // Overflow check 2^k
                      break;
                  }
                  k_power_of_2 = next_k_power; // 2^k for k>=1

                 // Calculate n = base * 2^k
                 uint64_t current_n = safe_multiply(base, k_power_of_2);

                 if (current_n == 0 || current_n > N) {
                     break; // Overflow or exceeded limit N
                 }
                 numbers_set.insert(current_n);
            }
             if (p > N) break; // Optimization (redundant check?)
        }
    } else if (choice == 4) {
         // Numbers of the form n = 2^k * p^m * q^n, where q > p (both odd primes), k>=0, m>=1, n>=1
         // Loop through all m >= 1 and n >= 1 such that the result is <= N
        for (size_t i = 0; i < primes.size(); ++i) {
            uint64_t p = primes[i];
            if (p < 3) continue; // p must be an odd prime
            if (p > N) break; // Optimization: if p > N, no power m>=1 will work

            uint64_t current_p_power_m = 1; // Start with p^0 to enter the loop easily
            for (int m = 1; ; ++m) { // Loop for exponent m >= 1
                // Calculate p^m for this iteration
                uint64_t next_p_power = safe_multiply(current_p_power_m, p);
                 if (next_p_power == 0 || next_p_power < current_p_power_m ) { // Overflow check for p^m
                    break; // Stop m loop due to p^m overflow
                 }
                 current_p_power_m = next_p_power; // This is p^m

                // Check if current_p_power_m (p^m) exceeds N
                if (current_p_power_m > N) {
                    break; // Stop increasing m for this p
                }

                // Now loop through q > p
                for (size_t j = i + 1; j < primes.size(); ++j) { // Iterate q > p
                    uint64_t q = primes[j];
                    // q is already > p and odd

                     // Optimization: Check if p^m * q^1 > N
                     uint64_t base_check = safe_multiply(current_p_power_m, q);
                     if (base_check == 0 || base_check > N) {
                         // If p^m * q > N, then p^m * q^n (n>=1) will also be > N.
                         // Also, any larger q will also result in p^m * q'^n > N.
                         // So, we can break the q loop for this p^m.
                         break; // Stop trying larger q for this p^m
                     }

                     uint64_t current_q_power_n = 1; // Start with q^0
                     for(int n_exp = 1; ; ++n_exp) { // Loop for exponent n >= 1
                         // Calculate q^n for this iteration
                         uint64_t next_q_power = safe_multiply(current_q_power_n, q);
                         if (next_q_power == 0 || next_q_power < current_q_power_n ) { // Overflow check for q^n
                             break; // Stop n loop due to q^n overflow
                         }
                         current_q_power_n = next_q_power; // This is q^n

                         // Calculate base = p^m * q^n
                         uint64_t base = safe_multiply(current_p_power_m, current_q_power_n);
                         if (base == 0 || base > N) {
                             // Overflow or exceeded limit N. Stop increasing n for this q.
                             break; // Stop inner loop (powers of q)
                         }

                        // Now multiply by powers of 2 (2^k)
                        // Add n = base * 2^0 first (k=0)
                        if (base <= N) { // Double check base <= N
                           numbers_set.insert(base);
                        } else {
                             // If base itself > N, no need for k loop. Also, break n loop.
                             // This case should theoretically be caught by the check before the k loop, but added for safety.
                              break;
                        }


                        uint64_t k_power_of_2 = 1; // Represents 2^k starting power
                        while (true) { // Loop for exponent k >= 1
                            // Prepare for next power of 2
                            uint64_t next_k_power_2 = safe_multiply(k_power_of_2, 2);
                            if (next_k_power_2 == 0 || next_k_power_2 < k_power_of_2) { // Overflow check for 2^k
                                break; // Stop k loop due to 2^k overflow
                            }
                            k_power_of_2 = next_k_power_2; // This is now 2^k for k>=1

                            // Calculate n = base * 2^k
                            uint64_t current_n = safe_multiply(base, k_power_of_2);

                            if (current_n == 0 || current_n > N) {
                                break; // Overflow or exceeded limit N, stop k loop
                            }
                            numbers_set.insert(current_n);
                        } // End k loop
                     } // End n loop
                } // End q loop
            } // End m loop
        } // End p loop
    } // End choice 4

    // Convert set to vector for return
    vector<uint64_t> numbers(numbers_set.begin(), numbers_set.end());
    // The set ensures elements are unique and sorted.
    return numbers;
}

/**
 * Function: getPrimeFactorizationString
 * -------------------------------------
 * Computes the prime factorization of a number n and returns it as a string.
 * Format: p1^a1 * p2^a2 * ...
 */
string getPrimeFactorizationString(uint64_t n, const vector<uint64_t>& primes) {
    if (n == 0) return "0";
    if (n == 1) return "1";

    stringstream ss;
    uint64_t temp_n = n;
    bool first_factor = true;

    // Factor out 2
    if (temp_n % 2 == 0) {
        int count = 0;
        while (temp_n > 0 && temp_n % 2 == 0) { // Ensure temp_n doesn't become 0 unexpectedly
            count++;
            temp_n /= 2;
        }
        if (count > 0) {
            if (!first_factor) ss << " * ";
            ss << "2^" << count;
            first_factor = false;
        }
    }

    // Factor out odd primes using the provided list
    uint64_t sqrt_temp_n_orig = (temp_n > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n))) : 0;
    for (const auto& p : primes) {
        if (p == 2) continue; // Already handled
        if (temp_n == 1) break; // Fully factored
        if (p > sqrt_temp_n_orig) break; // Optimization

        // Check if p*p > current temp_n as another optimization point?
        uint64_t p_squared = safe_multiply(p, p);
         if (p_squared == 0 && p > 1) { /* Overflow, rely on p > sqrt_temp_n */ }
         else if (p_squared > temp_n && temp_n > 1) {
             // If p*p > current temp_n, p can only be a factor if p == temp_n.
             // This will be handled by the check after the loop. We can't break the loop
             // itself based on this, only on p > sqrt(original temp_n).
         }


        if (temp_n % p == 0) {
            int count = 0;
            while (temp_n > 0 && temp_n % p == 0) {
                count++;
                temp_n /= p;
            }
             if (count > 0) {
                 if (!first_factor) ss << " * ";
                 ss << p << "^" << count;
                 first_factor = false;
                 // Update original sqrt check? Not needed if we break loop early enough
                 // sqrt_temp_n_orig = (temp_n > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n))) : 0;
             }
        }
    }

    // If temp_n is still greater than 1, it must be a prime factor itself
    if (temp_n > 1) {
        if (!first_factor) ss << " * ";
        ss << temp_n << "^1"; // Remaining factor is prime
    }

    // If factorization resulted in empty string (e.g., error or n=1 case missed)
    if (ss.str().empty() && n==1) return "1"; // Should be caught earlier
    if (ss.str().empty() && n > 1) return "ErrorFactorizing"; // Indicate failure

    return ss.str();
}


/**
 * Function: printTableHeader
 * --------------------------
 * Prints the header of the formatted table with vertical splitters.
 */
void printTableHeader() {
    // Define column widths
    const int width_n = 15;
    const int width_sigma = 22;
    const int width_factorization = 25; // New column for prime factorization
    const int width_diff = 12;
    const int width_combinations = 50;

    // Print the top border
    cout << "+" << string(width_n, '-')
         << "+" << string(width_sigma, '-')
         << "+" << string(width_factorization, '-')
         << "+" << string(width_diff, '-')
         << "+" << string(width_combinations, '-') << "+\n";

    // Print the header row
    cout << "| " << left << setw(width_n - 1) << "n"
         << "| " << left << setw(width_sigma - 1) << "sigma(n)"
         << "| " << left << setw(width_factorization - 1) << "Prime Factorization" // Updated Header
         << "| " << left << setw(width_diff - 1) << "diff"
         << "| " << left << setw(width_combinations - 1) << "Valid (d1, d2) Combinations"
         << "|\n";

    // Print the separator
    cout << "+" << string(width_n, '=')
         << "+" << string(width_sigma, '=')
         << "+" << string(width_factorization, '=') // Updated separator
         << "+" << string(width_diff, '=')
         << "+" << string(width_combinations, '=') << "+\n";
}

/**
 * Function: printTableRow
 * -----------------------
 * Prints a single row of the formatted table with vertical splitters.
 */
void printTableRow(uint64_t n, uint64_t sigma, const string& factorization, long long diff, const string& combinations) {
    // Define column widths (must match printTableHeader)
    const int width_n = 15;
    const int width_sigma = 22;
    const int width_factorization = 25;
    const int width_diff = 12;
    const int width_combinations = 50;

    // Limit factorization string length for display
    string display_factorization = factorization;
     if (display_factorization.length() > width_factorization - 3) {
        display_factorization = display_factorization.substr(0, width_factorization - 6) + "...";
    }


    // Limit combination string length for display
    string display_combinations = combinations;
    if (display_combinations.length() > width_combinations - 3) { // Account for "| " and " |"
        display_combinations = display_combinations.substr(0, width_combinations - 6) + "...";
    }


    cout << "| " << left << setw(width_n - 1) << n
         << "| " << left << setw(width_sigma - 1) << sigma
         << "| " << left << setw(width_factorization - 1) << display_factorization // Use factorization string
         << "| " << left << setw(width_diff - 1) << diff
         << "| " << left << setw(width_combinations - 1) << display_combinations
         << "|\n";
}

// Function to check if a string contains only digits
bool is_digits(const std::string &str) {
    return !str.empty() && str.find_first_not_of("0123456789") == std::string::npos;
}


int main() {
    uint64_t N;
    cout << "Enter the upper limit N (e.g., 1000000): ";
    while(!(cin >> N) || N < 1) {
        cout << "Please enter a positive integer greater than 0: ";
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    // Clear the input buffer
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    cout << "\nChoose the form to search for near perfect numbers:\n";
    cout << "1. All numbers (can be slow for large N)\n";
    cout << "2. Numbers of the form n = 2^k * p (where p is an odd prime, k>=0)\n";
    cout << "3. Numbers of the form n = 2^k * p^m (where p is an odd prime, k>=0, m>=1)\n";
    cout << "4. Numbers of the form n = 2^k * p^m * q^n (where p, q are odd primes, q>p, k>=0, m>=1, n>=1)\n"; // Updated description
    cout << "Enter 1, 2, 3, or 4 (Press Enter for default - Option 2): ";

    string choiceStr;
    getline(cin, choiceStr);
    int choice;
    int exponent_m = 1; // Default exponent m (used only for choice 3)
    // Exponent n is no longer needed here as choice 4 iterates all n>=1

    if(choiceStr.empty()) {
        choice = 2; // Defaulting to option 2
        cout << "No choice entered. Defaulting to option 2 (n = 2^k * p).\n";
    } else {
        try {
            choice = stoi(choiceStr);
            if(choice < 1 || choice > 4) {
                cout << "Invalid choice. Defaulting to option 2 (n = 2^k * p).\n";
                choice = 2;
            }
        } catch(const std::invalid_argument& ia) {
             cout << "Invalid input. Defaulting to option 2 (n = 2^k * p).\n";
             choice = 2;
        } catch(const std::out_of_range& oor) {
             cout << "Input out of range. Defaulting to option 2 (n = 2^k * p).\n";
             choice = 2;
        }
    }

    // If choice 3 is selected, prompt for exponent m
    if(choice == 3) {
        cout << "Enter the exponent m for p (form n = 2^k * p^m, m >= 1, default 1): ";
        string exponentMStr;
        getline(cin, exponentMStr);
        if(exponentMStr.empty()) {
            exponent_m = 1; // Default if empty
        } else if (!is_digits(exponentMStr)) {
            cout << "Invalid input for m (not a number). Defaulting to 1.\n";
            exponent_m = 1;
        } else {
            try {
                unsigned long long temp_m = stoull(exponentMStr);
                // Limit exponent to prevent easy overflow in safe_power
                if (temp_m < 1 || temp_m > 63) {
                     cout << "Exponent m must be between 1 and 63. Defaulting to 1.\n";
                     exponent_m = 1;
                 } else {
                    exponent_m = static_cast<int>(temp_m);
                 }
            } catch(const std::invalid_argument& ia) {
                cout << "Invalid input format for m. Defaulting to 1.\n";
                exponent_m = 1;
            } catch(const std::out_of_range& oor) {
                cout << "Exponent m out of range. Defaulting to 1.\n";
                exponent_m = 1;
            }
        }
        if (exponent_m < 1) exponent_m = 1; // Ensure m>=1
    }
    // No need to ask for exponents for choice 4 anymore


    // Prompt to exclude primes
    cout << "Do you want to exclude prime numbers from the output? (y/n, default n): ";
    string excludePrimesStr;
    getline(cin, excludePrimesStr);
    bool excludePrimes = false;
    if(!excludePrimesStr.empty()) {
        char response = tolower(excludePrimesStr[0]);
        if(response == 'y') {
            excludePrimes = true;
        }
    }

    // Prompt to select sign combination condition
    cout << "\nSelect the sign combinations to include for (d1, d2):\n";
    cout << "1. All combinations (++ , +-, -+, --)\n";
    cout << "2. Only mixed signs (+-, -+)\n";
    cout << "3. Only positive signs (++)\n";
    cout << "4. Only negative signs (--)\n";
    cout << "Enter 1, 2, 3, or 4 (Press Enter for default - All combinations): ";

    string signOptionStr;
    getline(cin, signOptionStr);
    SignOption signOption = ALL; // Default
    if(!signOptionStr.empty()) {
        try {
            int signChoice = stoi(signOptionStr);
            if(signChoice >=1 && signChoice <=4) {
                signOption = static_cast<SignOption>(signChoice);
            } else {
                cout << "Invalid choice. Defaulting to All combinations.\n";
                signOption = ALL;
            }
        } catch(const std::invalid_argument& ia) {
            cout << "Invalid input. Defaulting to All combinations.\n";
            signOption = ALL;
        } catch(const std::out_of_range& oor) {
             cout << "Input out of range. Defaulting to All combinations.\n";
             signOption = ALL;
        }
    }

    // Prompt to filter pairs (d1, d2) where |d1| == |d2|
    cout << "Do you want to exclude results where ALL valid (d1, d2) pairs have |d1| == |d2|? (y/n, default y): ";
    string filterSameDivisorStr;
    getline(cin, filterSameDivisorStr);
    bool filterSameAbsDivisorPairs = true; // Default to yes
    if (!filterSameDivisorStr.empty()) {
        char response = tolower(filterSameDivisorStr[0]);
        if (response == 'n') {
            filterSameAbsDivisorPairs = false;
        }
    }


    // Generate base primes up to sqrt(N) using simple sieve
    uint64_t sqrtN = (N > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(N))) + 1 : 2;
    // Clamp sqrtN to prevent excessive memory usage if N is huge
    uint64_t max_sieve_limit = 200000000; // e.g., Limit sieve to 200 million
    if (sqrtN > max_sieve_limit) {
        cout << "Warning: sqrt(N) is very large (" << sqrtN << "). Limiting prime sieve to " << max_sieve_limit << " for efficiency." << endl;
        sqrtN = max_sieve_limit;
    }
    cout << "\nGenerating base primes up to " << sqrtN << " using simple sieve...\n";
    vector<uint64_t> basePrimes = simpleSieve(sqrtN);
    cout << "Number of base primes found: " << basePrimes.size() << "\n";

    // Generate numbers to process based on user's choice
    cout << "Generating numbers up to N=" << N << " based on the selected form...\n";
    // Pass exponent_m only for choice 3 generation. generateNumbers ignores it for choice 4.
    vector<uint64_t> numbersToProcess = generateNumbers(N, choice, basePrimes, exponent_m);
    cout << "Number of candidate numbers generated: " << numbersToProcess.size() << "\n";
    if (numbersToProcess.empty() && N > 1) {
         cout << "Warning: No candidate numbers were generated. Check N and choice." << endl;
         // Potentially exit early? Or let the loop run (it will do nothing).
    }


    // Display the type of numbers being processed
    cout << "\nSearching for Near Perfect Numbers up to " << N;
    switch(choice) {
        case 1: cout << ":"; break;
        case 2: cout << " of the form n = 2^k * p (k>=0):"; break;
        case 3: cout << " of the form n = 2^k * p^" << exponent_m << " (k>=0, m=" << exponent_m << "):"; break;
        case 4: cout << " of the form n = 2^k * p^m * q^n (q>p, k>=0, m>=1, n>=1):"; break; // Updated message
    }
    if(excludePrimes) {
        cout << " (Primes excluded)";
    }
    cout << "\nApplying condition: sigma(n) = 2n + (sign1*d1) + (sign2*d2)";

    // Display selected sign combinations
    cout << "\nSign Combinations Included: ";
    switch(signOption) {
        case ALL: cout << "All (+ +, + -, - +, - -)"; break;
        case MIXED_SIGNS: cout << "Mixed (+ -, - +)"; break;
        case ONLY_POSITIVE: cout << "Positive (+ +)"; break;
        case ONLY_NEGATIVE: cout << "Negative (- -)"; break;
    }
    // Display filter status
    if (filterSameAbsDivisorPairs) {
         cout << "\nFiltering out numbers where all valid pairs have |d1| == |d2|.";
    } else {
         cout << "\nIncluding numbers even if all valid pairs have |d1| == |d2|.";
    }
    cout << "\n\n";


    // Print the table header
    printTableHeader();
    bool anyNearPerfect = false;

    // Open CSV file for writing (overwrites existing file)
    ofstream csvFile("near_perfect_numbers.csv");
    if(!csvFile.is_open()) {
        cerr << "Error: Could not create/write to 'near_perfect_numbers.csv'.\n";
        return 1;
    }

    // Write CSV headers
    csvFile << "n,sigma(n),Prime Factorization,diff,Valid (d1, d2) Combinations\n"; // Updated header

    // Processing numbers
    uint64_t processed = 0;
    size_t totalNumbers = numbersToProcess.size();
    // Calculate report interval, ensure it's at least 1
    const uint64_t report_interval = (totalNumbers > 100) ? std::max((uint64_t)1, (uint64_t)(totalNumbers / 100)) : 10; // Report less often if few numbers
     uint64_t next_report = report_interval;


    for(uint64_t n : numbersToProcess) {
        processed++;
        // Display progress
        if (totalNumbers > 0 && (processed >= next_report || processed == totalNumbers)) {
             double percent = (double)processed / totalNumbers * 100.0;
             cout << "\rProcessing: " << processed << " / " << totalNumbers << " (" << fixed << setprecision(1) << percent << "%)" << flush;
             next_report += report_interval;
        }

        // Skip 0 or 1 (sigma(0)=0, sigma(1)=1)
        if (n < 2) continue;

        // Skip prime numbers if exclusion is enabled
        if(excludePrimes) {
            bool is_n_prime = false;
            // Optimization: Check if n is in the basePrimes list (if n <= sqrtN)
            if (n <= sqrtN && !basePrimes.empty() && n >= basePrimes[0] ) {
                 if (binary_search(basePrimes.begin(), basePrimes.end(), n)) {
                     is_n_prime = true;
                 }
            }
            // If not found in basePrimes or n > sqrtN, run the full primality test
            // isPrime is only efficient if basePrimes covers up to sqrt(n)
            // If n > sqrtN*sqrtN, isPrime might need extension or a better algorithm (like Miller-Rabin)
            // For now, rely on the provided basePrimes. It might be slow if n is large and sqrtN was capped.
            if (!is_n_prime && n > (basePrimes.empty() ? 1 : basePrimes.back()) ) {
                // Only run isPrime if n wasn't found AND n is larger than the largest base prime
                // We pass basePrimes to isPrime for trial division optimization up to sqrt(n) or basePrimes.back()
                 is_n_prime = isPrime(n, basePrimes);
            }

            if(is_n_prime) {
                 continue; // Skip this prime number
            }
        }


        // Compute sigma(n)
        uint64_t sigma = computeSigma(n, basePrimes);
        if (sigma == 0 && n > 0) { // Check if computeSigma overflowed (sigma(0)=0 is valid)
             // cerr << "\nWarning: Sigma calculation overflowed for n=" << n << ". Skipping." << endl;
             continue;
        }

        // Calculate difference: sigma(n) - 2n (as long long)
        long long diff;
        uint64_t two_n_val = safe_multiply(n, 2);
        if (two_n_val == 0 && n != 0) { // Check 2*n overflow
            // cerr << "\nWarning: 2*n calculation overflowed for n=" << n << ". Skipping." << endl;
             continue;
        }

        // Check for potential overflow/underflow before subtraction
        if (sigma >= two_n_val) {
            uint64_t diff_abs = sigma - two_n_val;
             if (diff_abs > (uint64_t)numeric_limits<long long>::max()) {
                 //cerr << "\nWarning: Positive difference calculation overflowed long long for n=" << n << ". Skipping." << endl;
                 continue; // Difference too large positive to fit in long long
             }
             diff = static_cast<long long>(diff_abs);
        } else {
            // sigma < two_n_val, result is negative
            uint64_t diff_abs = two_n_val - sigma;
            // Check if absolute difference exceeds magnitude of LLONG_MIN
             if (diff_abs > static_cast<uint64_t>(numeric_limits<long long>::max()) + 1ULL ) { // Compare against max+1 magnitude
                 //cerr << "\nWarning: Negative difference calculation overflowed long long for n=" << n << ". Skipping." << endl;
                 continue; // Difference too large negative
             }
             // Need to cast carefully to avoid implementation-defined behavior for -LLONG_MIN
             if (diff_abs == static_cast<uint64_t>(numeric_limits<long long>::max()) + 1ULL) {
                  diff = numeric_limits<long long>::min();
             } else {
                 diff = -static_cast<long long>(diff_abs);
             }
        }


        // Compute divisors for finding combinations
        vector<uint64_t> divisors = computeDivisors(n, basePrimes);
        if (divisors.empty() && n > 1) { // Check if computeDivisors indicated overflow/error
             // cerr << "\nWarning: Divisor computation failed (likely overflow) for n=" << n << ". Skipping." << endl;
             continue;
        }


        // Find valid (d1, d2) combinations based on signOption
        vector<string> validPairs = findValidCombinations(divisors, diff, signOption);

        if(!validPairs.empty()) {
            // Apply the filter: exclude if ALL pairs have |d1| == |d2|
            if (filterSameAbsDivisorPairs) {
                bool hasDistinctPair = false;
                for(const auto &pairStr : validPairs) {
                    // Minimal parsing: find the space separating the two signed numbers
                    size_t first_space = pairStr.find(' ');
                    if (first_space != string::npos && first_space > 0 && first_space + 1 < pairStr.length()) {
                         size_t second_sign_pos = pairStr.find_first_of("+-", first_space + 1);
                         if (second_sign_pos != string::npos && second_sign_pos > first_space + 1 && second_sign_pos + 1 < pairStr.length()) {
                             string d1_abs_str = pairStr.substr(1, first_space - 1); // number after first sign
                             string d2_abs_str = pairStr.substr(second_sign_pos + 1); // number after second sign
                             if (d1_abs_str != d2_abs_str) {
                                  hasDistinctPair = true;
                                  break; // Found one distinct pair, no need to check further
                             }
                         } else { hasDistinctPair = true; break; } // Parsing failed - assume distinct
                    } else { hasDistinctPair = true; break; } // Parsing failed - assume distinct
                }
                // If after checking all pairs, none had distinct absolute values, skip this n.
                if (!hasDistinctPair) {
                    continue;
                }
            }

            // If we reached here, the number is a valid near-perfect number based on criteria
            anyNearPerfect = true;

            // Combine all valid pairs into a single string separated by "; "
            string combinationsCombinedStr;
            for(size_t i = 0; i < validPairs.size(); ++i) {
                combinationsCombinedStr += validPairs[i];
                if (i < validPairs.size() - 1) combinationsCombinedStr += "; ";
            }

            // Get prime factorization string
            string factorizationStr = getPrimeFactorizationString(n, basePrimes);

            // Print the row in the table
            printTableRow(n, sigma, factorizationStr, diff, combinationsCombinedStr);

            // Write the row to the CSV file
            // Enclose factorization and combinations in quotes to handle potential special characters
            csvFile << n << "," << sigma << ",\"" << factorizationStr << "\"," << diff << ",\"" << combinationsCombinedStr << "\"\n";
        }
    } // End number processing loop

    // Ensure the final progress update is cleared or finalized
    cout << "\rProcessing: " << totalNumbers << " / " << totalNumbers << " (100.0%)" << endl; // Final progress line


    // Print the bottom border of the table
    cout << "+" << string(15, '-') // width_n
         << "+" << string(22, '-') // width_sigma
         << "+" << string(25, '-') // width_factorization
         << "+" << string(12, '-') // width_diff
         << "+" << string(50, '-') // width_combinations
         << "+\n";

    // Close the CSV file
    csvFile.close();

    if(anyNearPerfect) {
        cout << "\nResults have been saved to 'near_perfect_numbers.csv'.\n";
    } else {
        cout << "\nNo near perfect numbers found matching the criteria";
        switch(choice) {
            case 1: cout << "."; break;
            case 2: cout << " of the form n = 2^k * p."; break;
            case 3: cout << " of the form n = 2^k * p^" << exponent_m << "."; break;
            case 4: cout << " of the form n = 2^k * p^m * q^n (q > p, m>=1, n>=1)."; break; // Updated message
        }
        if(excludePrimes) {
            cout << " (Primes were excluded)";
        }
        if (filterSameAbsDivisorPairs) {
             cout << " (Numbers where all pairs had |d1|=|d2| were excluded)";
        }
        cout << "\n";
    }

    return 0;
}