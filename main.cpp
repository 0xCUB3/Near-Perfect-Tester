#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm> // Needed for std::max, std::sort, std::unique
#include <limits>
#include <iomanip>
#include <fstream>
#include <cctype>    // For tolower
#include <cstdint>   // For uint64_t
#include <set>       // For unique divisor finding in computeDivisors
#include <stdexcept> // For stoi exceptions

using namespace std;

// Helper function to safely multiply two uint64_t, returning 0 on overflow
// Using __int128 for intermediate check
uint64_t safe_multiply(uint64_t a, uint64_t b) {
    // Early exit for 0 to avoid unnecessary __int128 calculation
    if (a == 0 || b == 0) return 0;
    // Check if multiplication would obviously overflow based on leading bits
    // (More advanced check, requires bit manipulation or could use compiler builtins)
    // Simple check: if a > max/b roughly
    if (a > numeric_limits<uint64_t>::max() / b) return 0;

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
            if (res == 0 && current_base != 0) return 0; // Overflow detected by safe_multiply
        }
        current_exp /= 2;
        if (current_exp > 0) { // Only need to square if exponent isn't done
             // Check before squaring: if base > sqrt(max), squaring will overflow
             if (current_base > (uint64_t)sqrt(numeric_limits<uint64_t>::max())) return 0;
             current_base = safe_multiply(current_base, current_base);
             if (current_base == 0) return 0; // Overflow detected by safe_multiply
        }
    }
    return res;
}


/**
 * Function: simpleSieve
 * ---------------------
 * Finds all prime numbers up to sqrtN using the simple sieve of Eratosthenes.
 *
 * Parameters:
 *  - limit: The upper limit up to which primes are identified.
 *
 * Returns:
 *  - A vector of primes up to limit.
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
 * Optimizes by checking against precomputed primes up to sqrt(n).
 * If n is larger than the square of the largest prime in `primes`, it finishes the check.
 *
 * Parameters:
 *  - n: The number to check for primality.
 *  - primes: A vector of known primes, ideally up to at least sqrt(n).
 *
 * Returns:
 *  - True if n is prime, False otherwise.
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
    for (uint64_t i = 5; i <= sqrt_n; i += 6) {
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
 *
 * Parameters:
 *  - n: The number for which sigma(n) is computed.
 *  - primes: A vector of primes up to sqrt(n) (can be smaller, but less efficient).
 *
 * Returns:
 *  - The sum of divisors of n. Returns 0 if intermediate calculation overflows.
 */
uint64_t computeSigma(uint64_t n, const vector<uint64_t>& primes) {
    if (n == 0) return 0;
    if (n == 1) return 1;

    uint64_t sigma = 1;
    uint64_t temp_n = n; // Use a temporary variable for modification

    // Handle factor 2 separately
    if (temp_n % 2 == 0) {
        uint64_t term = 1;
        uint64_t power = 1;
        while (temp_n % 2 == 0) {
            power = safe_multiply(power, 2);
            if (power == 0) return 0; // Overflow in power calculation
            #ifdef __SIZEOF_INT128__
            unsigned __int128 next_term = (unsigned __int128)term + power;
            if (next_term > numeric_limits<uint64_t>::max()) return 0; // Overflow in term sum
            term = (uint64_t)next_term;
            #else
            // Fallback check (less precise)
            if (term > numeric_limits<uint64_t>::max() - power) return 0;
            term += power;
            #endif
            temp_n /= 2;
        }
        sigma = safe_multiply(sigma, term);
        if (sigma == 0 && term != 0) return 0; // Overflow in sigma product
    }

    // Handle odd prime factors
    uint64_t sqrt_temp_n = (temp_n > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n))) : 0;
    for (const auto& p : primes) {
        if (p == 2) continue; // Already handled
        if (p > sqrt_temp_n) break; // Optimization: check primes up to sqrt(remaining n)
        if (p*p > temp_n) break; // Alternative optimization based on current prime (requires temp_n>0 check)

        if (temp_n % p == 0) {
            uint64_t term = 1;
            uint64_t power = 1;
            while (temp_n % p == 0) {
                power = safe_multiply(power, p);
                if (power == 0) return 0; // Overflow in power calculation
                 #ifdef __SIZEOF_INT128__
                 unsigned __int128 next_term = (unsigned __int128)term + power;
                 if (next_term > numeric_limits<uint64_t>::max()) return 0; // Overflow in term sum
                 term = (uint64_t)next_term;
                 #else
                 if (term > numeric_limits<uint64_t>::max() - power) return 0;
                 term += power;
                 #endif
                temp_n /= p;
            }
             sigma = safe_multiply(sigma, term);
             if (sigma == 0 && term != 0) return 0; // Overflow in sigma product
            // Update sqrt_temp_n as temp_n decreases
             sqrt_temp_n = (temp_n > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n))) : 0;
        }
    }

    // If temp_n is still greater than 1, it must be a prime factor itself
    if (temp_n > 1) {
        // Calculate term = 1 + temp_n safely
        uint64_t term;
        if (temp_n == numeric_limits<uint64_t>::max()) return 0; // 1 + max overflows
        term = temp_n + 1;
        sigma = safe_multiply(sigma, term);
        if (sigma == 0 && term != 0) return 0; // Overflow in sigma product
    }

    return sigma;
}


/**
 * Function: computeDivisors
 * -------------------------
 * Computes all distinct divisors of a number n.
 *
 * Parameters:
 *  - n: The number for which divisors are computed.
 *  - primes: A vector of primes up to sqrt(n) (can be smaller, but less efficient).
 *
 * Returns:
 *  - A sorted vector of distinct divisors of n. Returns empty vector if overflow occurs during calculation.
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
        }
    };

    // Handle factor 2
    add_prime_power_divisors(2);
    if (overflow_detected) return {}; // Return empty on overflow

    // Handle odd prime factors
    uint64_t sqrt_temp_n = (temp_n > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n))) : 0;
    for (const auto& p : primes) {
        if (p == 2) continue;
        if (p > sqrt_temp_n) break;
        if (temp_n == 1) break; // Optimization: if temp_n reduced to 1, no more factors
        // Use p*p check to avoid repeated sqrt calculation inside loop
        // Need to be careful with large p near sqrt(ULLONG_MAX)
        uint64_t p_squared = safe_multiply(p, p);
        if (p_squared == 0 && p > 1) { // p*p overflowed
            // Cannot use p*p check, rely on p > sqrt_temp_n check
        } else if (p_squared > temp_n) {
            break; // If p*p > remaining temp_n, p cannot be a factor (unless p = temp_n)
        }


        if (temp_n % p == 0) {
            add_prime_power_divisors(p);
             if (overflow_detected) return {}; // Return empty on overflow
            // Update sqrt_temp_n as temp_n decreases
             sqrt_temp_n = (temp_n > 1) ? static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n))) : 0;
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
 * ----------------
 * Enumerates the possible sign combination options.
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
 * Finds all valid combinations of divisors d1 and d2 with their respective signs
 * that satisfy the equation: sigma(n) = 2n ± d1 ± d2
 *
 * Parameters:
 *  - D: Reference to the vector containing all distinct divisors of n.
 *  - diff: The difference sigma(n) - 2n.
 *  - signOption: The user's choice for sign combinations.
 *
 * Returns:
 *  - A vector of strings, each representing a valid combination of d1 and d2 with signs.
 *
 * Sign Combinations:
 *  - ALL: Includes ++, +-, -+, --
 *  - MIXED_SIGNS: Includes +-, -+
 *  - ONLY_POSITIVE: Includes ++
 *  - ONLY_NEGATIVE: Includes --
 */
vector<string> findValidCombinations(const vector<uint64_t> &D, long long diff, SignOption signOption) {
    vector<string> validPairs;
    set<string> uniquePairs; // Use set to automatically handle duplicates and sorting

    // Iterate through all possible pairs (d1, d2)
    for (uint64_t d1 : D) {
        for (uint64_t d2 : D) {
            // Skip if any of the divisors is zero (shouldn't happen with computeDivisors)
            if (d1 == 0 || d2 == 0) {
                continue;
            }

            // Cast to long long for safe arithmetic with diff
            // Note: This cast itself can be problematic if d1/d2 > LLONG_MAX
            // However, diff is long long, so comparisons should be okay if d1/d2 are large positive.
            // The checks below handle potential overflow in the *combination* arithmetic.
            long long ld1 = (d1 <= (uint64_t)numeric_limits<long long>::max()) ? static_cast<long long>(d1) : -1; // Mark as invalid if too large
            long long ld2 = (d2 <= (uint64_t)numeric_limits<long long>::max()) ? static_cast<long long>(d2) : -1;

            if (ld1 == -1 || ld2 == -1) continue; // Skip if divisor exceeds long long max


            // Check combinations based on signOption
            if (signOption == ALL || signOption == ONLY_POSITIVE) {
                // ++
                 // Check for potential overflow before addition: ld1 + ld2 > LLONG_MAX ?
                 if (ld1 <= numeric_limits<long long>::max() - ld2) {
                     if ((ld1 + ld2) == diff) {
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
                if ((-ld1 + ld2) == diff) {
                    uniquePairs.insert("-" + to_string(d1) + " +" + to_string(d2));
                }
            }
            if (signOption == ALL || signOption == ONLY_NEGATIVE) {
                // --
                 // Check for potential underflow before subtraction: -ld1 - ld2 < LLONG_MIN ?
                 // Equivalent to: -ld1 < LLONG_MIN + ld2
                 if (ld2 > numeric_limits<long long>::min() + ld1) { // Check avoids direct LLONG_MIN + ld2 if ld2 is large neg
                     if ((-ld1 - ld2) == diff) {
                          uniquePairs.insert("-" + to_string(d1) + " -" + to_string(d2));
                     }
                 } else if (ld1 == 0 && ld2 == 0 && diff == 0) { // Handle 0 case specifically if needed
                     // This case shouldn't occur if d1, d2 > 0
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
 *
 * Parameters:
 *  - N: The upper limit up to which numbers are generated.
 *  - choice: The form selected by the user.
 *  - primes: Vector containing relevant primes (e.g., up to sqrt(N)).
 *  - exponent_m: The exponent m for the prime p (only relevant for choice 3 and 4).
 *  - exponent_n: The exponent n for the prime q (only relevant for choice 4).
 *
 * Returns:
 *  - A vector of unique, sorted integers representing the numbers to process.
 */
vector<uint64_t> generateNumbers(uint64_t N, int choice, const vector<uint64_t> &primes, int exponent_m = 1, int exponent_n = 1) {
    set<uint64_t> numbers_set; // Use set to avoid duplicates automatically

    if (choice == 1) {
        cout << "\nWarning: Generating all numbers up to N=" << N << " can be very slow and memory intensive." << endl;
        cout << "Consider using a smaller N or choosing a specific form (2, 3, or 4)." << endl;
        // All numbers from 1 to N
        for (uint64_t n = 1; n <= N; ++n) {
            // Check if n fits the structure if needed, although choice 1 implies all numbers
            numbers_set.insert(n);
            if (n == numeric_limits<uint64_t>::max()) break; // Prevent overflow in loop
        }
    } else if (choice == 2) {
        // Numbers of the form n = 2^k * p
        for (size_t i = 0; i < primes.size(); ++i) {
            uint64_t p = primes[i];
            if (p < 3) continue; // Skip p=2

            uint64_t k_power_of_2 = 1; // Represents 2^k
            while(true) {
                // Calculate n = p * 2^k
                 uint64_t current_n = safe_multiply(p, k_power_of_2);

                 if (current_n == 0 || current_n > N) {
                     // Overflow or exceeded limit N. Stop for this p.
                     // If k_power_of_2 was already 1, and p*1 > N, then break the outer loop too.
                     if (k_power_of_2 == 1 && (current_n > N || p > N)) break; // Optimization
                     break; // Stop inner loop (powers of 2)
                 }

                 numbers_set.insert(current_n);

                 // Prepare for next power of 2
                 uint64_t next_k_power = safe_multiply(k_power_of_2, 2);
                 if (next_k_power == 0 || next_k_power < k_power_of_2 ) { // Overflow check
                     break;
                 }
                 k_power_of_2 = next_k_power;
            }
             // Optimization: if p itself is > N, subsequent primes will be too.
             // Need to be careful if N is very large, p could exceed sqrt(N) but still be < N.
             if (p > N) break;
        }
    } else if (choice == 3) {
        // Numbers of the form n = 2^k * p^m
         for (size_t i = 0; i < primes.size(); ++i) {
            uint64_t p = primes[i];
            if (p < 3) continue; // Skip p=2

            uint64_t p_power_m = safe_power(p, exponent_m);
            // If p^m overflows or is > N, stop for this p and subsequent larger primes
            if (p_power_m == 0 || p_power_m > N) {
                // Optimization check: If p > N^(1/m), then p^m will exceed N.
                // Approximate check: If p*p > N and m >= 2, we can often break early.
                 uint64_t p_sq = safe_multiply(p,p);
                 if (p_sq == 0 && p > 1) {} // p*p overflowed, can't use this check
                 else if (p_sq > N && exponent_m >=2) break;

                 if (p > N && exponent_m == 1) break; // Break if p > N for m=1

                 // If p_power_m overflowed or > N, continue to next prime p, unless we can break outer loop
                 continue;
            }


            uint64_t k_power_of_2 = 1; // Represents 2^k
            while(true) {
                 // Calculate n = p^m * 2^k
                 uint64_t current_n = safe_multiply(p_power_m, k_power_of_2);

                 if (current_n == 0 || current_n > N) {
                     break; // Overflow or exceeded limit N
                 }
                 numbers_set.insert(current_n);

                 // Prepare for next power of 2
                  uint64_t next_k_power = safe_multiply(k_power_of_2, 2);
                  if (next_k_power == 0 || next_k_power < k_power_of_2 ) { // Overflow check
                      break;
                  }
                  k_power_of_2 = next_k_power;
            }
             if (p > N) break; // Optimization
        }
    } else if (choice == 4) {
        // Check if exponent_n is 0. If so, this case is redundant with Choice 3.
        if (exponent_n <= 0) { // Changed to <= 0 for robustness
            cerr << "\nWarning: Exponent n must be >= 1 for choice 4. Skipping generation for this choice." << endl;
        } else {
            // Numbers of the form n = 2^k * p^m * q^n, where q > p (both odd primes)
            for (size_t i = 0; i < primes.size(); ++i) {
                uint64_t p = primes[i];
                if (p < 3) continue; // p must be an odd prime

                uint64_t p_power_m = safe_power(p, exponent_m);
                if (p_power_m == 0) continue; // Overflow computing p^m, try next p
                if (p_power_m > N) break;     // Optimization: if p^m > N, no larger p or any q will work

                for (size_t j = i + 1; j < primes.size(); ++j) { // Iterate q > p
                    uint64_t q = primes[j];
                    // q is already > p and odd

                    uint64_t q_power_n = safe_power(q, exponent_n);
                    if (q_power_n == 0) continue; // Overflow computing q^n, try next q

                    // Calculate base = p^m * q^n
                    uint64_t base = safe_multiply(p_power_m, q_power_n);
                    if (base == 0 || base > N) {
                         // If base > N, increasing q further won't help for this p
                         // If base calculation overflows, increasing q further will likely also overflow.
                         break; // Stop trying larger q for this p
                    }

                    // Now multiply by powers of 2 (2^k)
                    uint64_t k_power_of_2 = 1; // Represents 2^k
                    while (true) {
                        // Calculate n = base * 2^k
                        uint64_t current_n = safe_multiply(base, k_power_of_2);

                        if (current_n == 0 || current_n > N) {
                            break; // Overflow or exceeded limit N
                        }
                        numbers_set.insert(current_n);

                        // Prepare for next power of 2
                        uint64_t next_k_power = safe_multiply(k_power_of_2, 2);
                        if (next_k_power == 0 || next_k_power < k_power_of_2) { // Overflow check
                            break;
                        }
                        k_power_of_2 = next_k_power;
                    } // End k loop
                } // End q loop

                // Optimization: If p^m * 3^n (smallest q^n) is already large, break outer loop early
                 uint64_t smallest_q_power_n = safe_power(3, exponent_n); // Smallest q is 3
                 if (smallest_q_power_n > 0) { // Check if 3^n didn't overflow
                     uint64_t min_base = safe_multiply(p_power_m, smallest_q_power_n);
                     // If p^m * 3^n > N, no q will work for this p, so break p loop
                     if (min_base > N && min_base != 0) {
                         break;
                     }
                 }

            } // End p loop
        } // End else block for exponent_n > 0 check
    } // End choice 4

    // Convert set to vector for return
    vector<uint64_t> numbers(numbers_set.begin(), numbers_set.end());
    // The set ensures elements are unique and sorted.
    return numbers;
}


/**
 * Function: printTableHeader
 * --------------------------
 * Prints the header of the formatted table with vertical splitters.
 */
void printTableHeader() {
    // Define column widths
    const int width_n = 15; // Increased width for larger N
    const int width_sigma = 22; // Increased width
    const int width_2n = 22;    // Increased width
    const int width_diff = 12;  // Increased width
    const int width_combinations = 50;

    // Print the top border
    cout << "+" << string(width_n, '-')
         << "+" << string(width_sigma, '-')
         << "+" << string(width_2n, '-')
         << "+" << string(width_diff, '-')
         << "+" << string(width_combinations, '-') << "+\n";

    // Print the header row
    cout << "| " << left << setw(width_n - 1) << "n"
         << "| " << left << setw(width_sigma - 1) << "sigma(n)"
         << "| " << left << setw(width_2n - 1) << "2n"
         << "| " << left << setw(width_diff - 1) << "diff"
         << "| " << left << setw(width_combinations - 1) << "Valid (d1, d2) Combinations"
         << "|\n";

    // Print the separator
    cout << "+" << string(width_n, '=')
         << "+" << string(width_sigma, '=')
         << "+" << string(width_2n, '=')
         << "+" << string(width_diff, '=')
         << "+" << string(width_combinations, '=') << "+\n";
}

/**
 * Function: printTableRow
 * -----------------------
 * Prints a single row of the formatted table with vertical splitters.
 *
 * Parameters:
 *  - n: The number.
 *  - sigma: The sum of divisors of n.
 *  - two_n: Twice the number.
 *  - diff: The difference sigma(n) - 2n.
 *  - combinations: The valid (d1, d2) combinations as a string.
 */
void printTableRow(uint64_t n, uint64_t sigma, uint64_t two_n, long long diff, const string& combinations) {
    // Define column widths (must match printTableHeader)
    const int width_n = 15;
    const int width_sigma = 22;
    const int width_2n = 22;
    const int width_diff = 12;
    const int width_combinations = 50;

    // Limit combination string length for display
    string display_combinations = combinations;
    if (display_combinations.length() > width_combinations - 3) { // Account for "| " and " |"
        display_combinations = display_combinations.substr(0, width_combinations - 6) + "...";
    }


    cout << "| " << left << setw(width_n - 1) << n
         << "| " << left << setw(width_sigma - 1) << sigma
         << "| " << left << setw(width_2n - 1) << two_n
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
    cout << "4. Numbers of the form n = 2^k * p^m * q^n (where p, q are odd primes, q>p, k>=0, m>=1, n>=1)\n";
    cout << "Enter 1, 2, 3, or 4 (Press Enter for default - Option 2): ";

    string choiceStr;
    getline(cin, choiceStr);
    int choice;
    int exponent_m = 1; // Default exponent m
    int exponent_n = 1; // Default exponent n

    if(choiceStr.empty()) {
        choice = 2; // Defaulting to option 2 as 'All numbers' is often too slow
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

    // If choice 3 or 4 is selected, prompt for exponents
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
                if (temp_m < 1 || temp_m > 60) { // Exponent > ~60 likely overflows uint64_t even for p=3
                     cout << "Exponent m must be between 1 and 60. Defaulting to 1.\n";
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
    } else if (choice == 4) {
         cout << "Enter the exponent m for p (form n = 2^k * p^m * q^n, m >= 1, default 1): ";
         string exponentMStr;
         getline(cin, exponentMStr);
         if(exponentMStr.empty()) {
             exponent_m = 1;
         } else if (!is_digits(exponentMStr)) {
             cout << "Invalid input for m (not a number). Defaulting to 1.\n";
             exponent_m = 1;
         } else {
             try {
                 unsigned long long temp_m = stoull(exponentMStr);
                  if (temp_m < 1 || temp_m > 60) {
                     cout << "Exponent m must be between 1 and 60. Defaulting to 1.\n";
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

         cout << "Enter the exponent n for q (form n = 2^k * p^m * q^n, n >= 1, default 1): ";
         string exponentNStr;
         getline(cin, exponentNStr);
         if(exponentNStr.empty()) {
             exponent_n = 1;
         } else if (!is_digits(exponentNStr)) {
             cout << "Invalid input for n (not a number). Defaulting to 1.\n";
             exponent_n = 1;
         } else {
             try {
                 unsigned long long temp_n = stoull(exponentNStr);
                 if (temp_n < 1 || temp_n > 60) {
                     cout << "Exponent n must be between 1 and 60. Defaulting to 1.\n";
                     exponent_n = 1;
                 } else {
                     exponent_n = static_cast<int>(temp_n);
                 }
             } catch(const std::invalid_argument& ia) {
                 cout << "Invalid input format for n. Defaulting to 1.\n";
                 exponent_n = 1;
             } catch(const std::out_of_range& oor) {
                 cout << "Exponent n out of range. Defaulting to 1.\n";
                 exponent_n = 1;
             }
         }
          // Double check in case user entered 0 and bypassed digit check somehow
          if (exponent_n < 1) exponent_n = 1;
          if (exponent_m < 1) exponent_m = 1;

    }


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
    cout << "3. Only positive signs (++ )\n";
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
    cout << "\nGenerating base primes up to " << sqrtN << " using simple sieve...\n";
    vector<uint64_t> basePrimes = simpleSieve(sqrtN);
    cout << "Number of base primes found: " << basePrimes.size() << "\n";

    // Generate numbers to process based on user's choice
    cout << "Generating numbers up to N=" << N << " based on the selected form...\n";
    vector<uint64_t> numbersToProcess = generateNumbers(N, choice, basePrimes, exponent_m, exponent_n);
    cout << "Number of candidate numbers generated: " << numbersToProcess.size() << "\n";

    // Display the type of numbers being processed
    cout << "\nSearching for Near Perfect Numbers up to " << N;
    switch(choice) {
        case 1: cout << ":"; break;
        case 2: cout << " of the form n = 2^k * p (k>=0):"; break;
        case 3: cout << " of the form n = 2^k * p^" << exponent_m << " (k>=0, m=" << exponent_m << "):"; break;
        case 4: cout << " of the form n = 2^k * p^" << exponent_m << " * q^" << exponent_n << " (q>p, k>=0, m=" << exponent_m << ", n=" << exponent_n << "):"; break;
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
    csvFile << "n,sigma(n),2n,diff,Valid (d1, d2) Combinations\n";

    // Processing numbers
    uint64_t processed = 0;
    size_t totalNumbers = numbersToProcess.size();
    // Calculate report interval, ensure it's at least 1
    const uint64_t report_interval = (totalNumbers > 100) ? std::max((uint64_t)1, (uint64_t)(totalNumbers / 100)) : 1;

    for(uint64_t n : numbersToProcess) {
        processed++;
        // Display progress
        if (totalNumbers > 0 && (processed % report_interval == 0 || processed == totalNumbers)) {
             double percent = (double)processed / totalNumbers * 100.0;
             cout << "\rProcessing: " << processed << " / " << totalNumbers << " (" << fixed << setprecision(1) << percent << "%)" << flush;
        }

        // Skip 0 or 1 (sigma(0)=0, sigma(1)=1)
        if (n < 2) continue;

        // Skip prime numbers if exclusion is enabled
        // Note: isPrime can be slow for very large n if basePrimes is relatively small
        if(excludePrimes) {
            // Optimization: if n is in the basePrimes list (up to sqrtN), it's prime
            // This is faster than running the full isPrime check for smaller primes
            bool maybe_prime = false;
            if (n <= sqrtN) {
                 if (binary_search(basePrimes.begin(), basePrimes.end(), n)) {
                     maybe_prime = true;
                 }
            }
            // If not found in basePrimes or n > sqrtN, run the full check
            if (maybe_prime || (n > sqrtN && isPrime(n, basePrimes))) {
                 continue;
            }
        }


        // Compute sigma(n)
        uint64_t sigma = computeSigma(n, basePrimes);
        if (sigma == 0 && n > 0) { // Check if computeSigma overflowed (sigma(0)=0 is valid)
             // cerr << "\nWarning: Sigma calculation overflowed for n=" << n << ". Skipping." << endl;
             continue;
        }

        // Calculate 2n safely
        uint64_t two_n = safe_multiply(n, 2);
         if (two_n == 0 && n != 0) { // Check overflow from safe_multiply
             // cerr << "\nWarning: 2n calculation overflowed for n=" << n << ". Skipping." << endl;
             continue;
         }


        // Calculate difference: sigma(n) - 2n (as long long)
        long long diff;
        // Check for potential overflow/underflow before subtraction
        if (sigma >= two_n) {
            uint64_t diff_abs = sigma - two_n;
             if (diff_abs > (uint64_t)numeric_limits<long long>::max()) {
                 //cerr << "\nWarning: Positive difference calculation overflowed long long for n=" << n << ". Skipping." << endl;
                 continue; // Difference too large positive to fit in long long
             }
             diff = static_cast<long long>(diff_abs);
        } else {
            // sigma < two_n, result is negative
            uint64_t diff_abs = two_n - sigma;
            // Check if absolute difference exceeds magnitude of LLONG_MIN
             if (diff_abs > (uint64_t)numeric_limits<long long>::max() + 1ULL ) { // Compare against max+1 magnitude
                 //cerr << "\nWarning: Negative difference calculation overflowed long long for n=" << n << ". Skipping." << endl;
                 continue; // Difference too large negative
             }
            diff = -static_cast<long long>(diff_abs);
        }


        // Compute divisors for finding combinations
        // This can be the bottleneck for large n
        vector<uint64_t> divisors = computeDivisors(n, basePrimes);
        if (divisors.empty() && n > 1) { // Check if computeDivisors indicated overflow (n=1 has divisor {1})
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
                    // Ensure space exists and there's something after it
                    if (first_space != string::npos && first_space + 1 < pairStr.length()) {
                         size_t second_sign_pos = pairStr.find_first_of("+-", first_space + 1);
                         // Ensure second sign exists and there's something after it
                         if (second_sign_pos != string::npos && second_sign_pos + 1 < pairStr.length()) {
                             string d1_abs_str = pairStr.substr(1, first_space - 1); // number after first sign
                             string d2_abs_str = pairStr.substr(second_sign_pos + 1); // number after second sign
                             if (d1_abs_str != d2_abs_str) {
                                  hasDistinctPair = true;
                                  break; // Found one distinct pair, no need to check further
                             }
                         } else {
                              // Parsing failed (e.g., no second sign or number) - treat as potentially distinct
                              hasDistinctPair = true; break;
                         }
                    } else {
                         // Parsing failed (e.g., no space) - treat as potentially distinct
                         hasDistinctPair = true; break;
                    }
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

            // Print the row in the table
            printTableRow(n, sigma, two_n, diff, combinationsCombinedStr);

            // Write the row to the CSV file
            // Enclose combinations in quotes to handle commas or potential special characters
            csvFile << n << "," << sigma << "," << two_n << "," << diff << ",\"" << combinationsCombinedStr << "\"\n";
        }
    } // End number processing loop

    // Ensure the final progress update is cleared or finalized
    cout << "\rProcessing: " << totalNumbers << " / " << totalNumbers << " (100.0%)" << endl; // Final progress line


    // Print the bottom border of the table
    cout << "+" << string(15, '-') // Match updated width_n
         << "+" << string(22, '-') // Match updated width_sigma
         << "+" << string(22, '-') // Match updated width_2n
         << "+" << string(12, '-') // Match updated width_diff
         << "+" << string(50, '-') << "+\n";

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
            case 4:
                 if (exponent_n > 0) { // Check if choice 4 was actually run
                    cout << " of the form n = 2^k * p^" << exponent_m << " * q^" << exponent_n << " (q > p).";
                 } else {
                    cout << " (Choice 4 skipped as exponent n was not >= 1).";
                 }
                 break;
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