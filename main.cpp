#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cctype> // For tolower
#include <cstdint> // For uint64_t
#include <set> // For unique divisor finding in computeDivisors

using namespace std;

// Helper function to safely multiply two uint64_t, returning 0 on overflow
// Using __int128 for intermediate check
uint64_t safe_multiply(uint64_t a, uint64_t b) {
    unsigned __int128 result = (unsigned __int128)a * b;
    if (result > numeric_limits<uint64_t>::max()) {
        return 0; // Indicate overflow
    }
    return (uint64_t)result;
}

// Helper function to safely raise uint64_t to a power, returning 0 on overflow
uint64_t safe_power(uint64_t base, int exp) {
    uint64_t res = 1;
    if (exp < 0) return 0; // Negative exponent not handled here
    if (exp == 0) return (base == 0) ? 0 : 1; // 0^0 is undefined, let's return 1, but 0^positive is 0. x^0 is 1 for x!=0

    for (int i = 0; i < exp; ++i) {
         unsigned __int128 temp_res = (unsigned __int128)res * base;
         if (temp_res > numeric_limits<uint64_t>::max()) {
             return 0; // Overflow
         }
         res = (uint64_t)temp_res;
    }
    return res;
}


/**
 * Function: simpleSieve
 * ---------------------
 * Finds all prime numbers up to sqrtN using the simple sieve of Eratosthenes.
 *
 * Parameters:
 *  - sqrtN: The upper limit up to which primes are identified.
 *
 * Returns:
 *  - A vector of primes up to sqrtN.
 */
vector<uint64_t> simpleSieve(uint64_t sqrtN) {
    uint64_t limit = sqrtN + 1;
    vector<bool> isPrime(limit, true);
    isPrime[0] = isPrime[1] = false;

    for (uint64_t i = 2; i * i <= sqrtN; ++i) {
        if (isPrime[i]) {
            for (uint64_t j = i * i; j <= sqrtN; j += i) {
                isPrime[j] = false;
            }
        }
    }

    vector<uint64_t> primes;
    // Reserve approximate space, avoid frequent reallocations
    if (sqrtN > 10) { // log(1) is problematic, avoid for small N
        primes.reserve(static_cast<size_t>(1.2 * sqrtN / log(sqrtN))); // Approximation + buffer
    } else {
        primes.reserve(5); // Small number
    }
    for (uint64_t i = 2; i <= sqrtN; ++i) {
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
        if (n % p == 0) return false;
    }

    // If sqrt_n is larger than the largest prime checked, continue checking
    // This part is less efficient if primes vector is too small, but ensures correctness.
    // Optimized trial division starting from the next prime after the largest in `primes` or 5
    uint64_t start_check = primes.empty() ? 5 : primes.back() + (primes.back() % 6 == 1 ? 4 : 2) ;
    if (start_check < 5) start_check = 5; // Ensure we start at 5 or later

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
            unsigned __int128 next_term = (unsigned __int128)term + power;
            if (next_term > numeric_limits<uint64_t>::max()) return 0; // Overflow in term sum
            term = (uint64_t)next_term;
            temp_n /= 2;
        }
        sigma = safe_multiply(sigma, term);
        if (sigma == 0) return 0; // Overflow in sigma product
    }

    // Handle odd prime factors
    uint64_t sqrt_temp_n = static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n)));
    for (const auto& p : primes) {
        if (p == 2) continue; // Already handled
        if (p > sqrt_temp_n) break; // Optimization: check primes up to sqrt(remaining n)
        if (p*p > temp_n) break; // Alternative optimization based on current prime

        if (temp_n % p == 0) {
            uint64_t term = 1;
            uint64_t power = 1;
            while (temp_n % p == 0) {
                power = safe_multiply(power, p);
                if (power == 0) return 0; // Overflow in power calculation
                unsigned __int128 next_term = (unsigned __int128)term + power;
                 if (next_term > numeric_limits<uint64_t>::max()) return 0; // Overflow in term sum
                term = (uint64_t)next_term;
                temp_n /= p;
            }
             sigma = safe_multiply(sigma, term);
             if (sigma == 0) return 0; // Overflow in sigma product
            // Update sqrt_temp_n as temp_n decreases
            sqrt_temp_n = static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n)));
        }
    }

    // If temp_n is still greater than 1, it must be a prime factor itself
    if (temp_n > 1) {
        uint64_t term = 1 + temp_n; // sigma(p) = 1 + p
        if (term < temp_n) return 0; // Check for overflow in 1+p
        sigma = safe_multiply(sigma, term);
        if (sigma == 0) return 0; // Overflow in sigma product
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
 *  - A sorted vector of distinct divisors of n.
 */
vector<uint64_t> computeDivisors(uint64_t n, const vector<uint64_t>& primes) {
    set<uint64_t> divisors_set;
    divisors_set.insert(1);
    if (n == 0) return {};
    if (n == 1) return {1};

    uint64_t temp_n = n;

    // Function to add multiples of a prime power to the divisor set
    auto add_prime_power_divisors = [&](uint64_t p) {
        vector<uint64_t> current_divisors(divisors_set.begin(), divisors_set.end());
        uint64_t power = 1;
        while (temp_n % p == 0) {
            power = safe_multiply(power, p);
            if (power == 0) { // Check for overflow
                cerr << "\nWarning: Overflow computing power for divisor generation at n=" << n << ", p=" << p << ". Results might be incomplete." << endl;
                return false; // Indicate failure due to overflow
            }
            for (uint64_t d : current_divisors) {
                 uint64_t new_divisor = safe_multiply(d, power);
                 if (new_divisor == 0) {
                     cerr << "\nWarning: Overflow computing new divisor for n=" << n << ". Results might be incomplete." << endl;
                     continue; // Skip this divisor but continue others
                 }
                 divisors_set.insert(new_divisor);
            }
            temp_n /= p;
        }
        return true; // Success
    };

    // Handle factor 2
    if (!add_prime_power_divisors(2)) return {}; // Return empty on overflow

    // Handle odd prime factors
    uint64_t sqrt_temp_n = static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n)));
    for (const auto& p : primes) {
        if (p == 2) continue;
        if (p > sqrt_temp_n) break;
        if (p*p > temp_n) break;

        if (temp_n % p == 0) {
            if (!add_prime_power_divisors(p)) return {}; // Return empty on overflow
            // Update sqrt_temp_n as temp_n decreases
            sqrt_temp_n = static_cast<uint64_t>(sqrt(static_cast<long double>(temp_n)));
        }
    }

    // If temp_n is still greater than 1, it's a prime factor
    if (temp_n > 1) {
        if (!add_prime_power_divisors(temp_n)) return {}; // Return empty on overflow
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
            long long ld1 = static_cast<long long>(d1);
            long long ld2 = static_cast<long long>(d2);

            // Check combinations based on signOption
            if (signOption == ALL || signOption == ONLY_POSITIVE) {
                // ++
                if ((ld1 + ld2) == diff) {
                    uniquePairs.insert("+" + to_string(d1) + " +" + to_string(d2));
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
                if ((-ld1 - ld2) == diff) {
                     uniquePairs.insert("-" + to_string(d1) + " -" + to_string(d2));
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
            numbers_set.insert(n);
            if (n == numeric_limits<uint64_t>::max()) break; // Prevent overflow in loop
        }
    } else if (choice == 2) {
        // Numbers of the form n = 2^k * p
        for (size_t i = 0; i < primes.size(); ++i) {
            uint64_t p = primes[i];
            if (p < 3) continue; // Skip p=2

            uint64_t current_n = p; // Start with k=0 (2^0 * p)
            if (current_n > N) break; // If p itself > N, no need to continue with this p

            uint64_t k_power_of_2 = 1; // Represents 2^k
            while(true) {
                // Calculate n = p * 2^k
                 unsigned __int128 n_128 = static_cast<unsigned __int128>(p) * k_power_of_2;

                 if (n_128 == 0 || n_128 > N) {
                     break; // Overflow or exceeded limit N
                 }

                 numbers_set.insert(static_cast<uint64_t>(n_128));

                 // Prepare for next power of 2
                 unsigned __int128 next_k_power = (unsigned __int128)k_power_of_2 * 2;
                 if (next_k_power > numeric_limits<uint64_t>::max()) break; // Next power of 2 overflows uint64_t
                 k_power_of_2 = (uint64_t)next_k_power;
                 if (k_power_of_2 == 0) break; // Should not happen if check above works, but safety first
            }
        }
    } else if (choice == 3) {
        // Numbers of the form n = 2^k * p^m
         for (size_t i = 0; i < primes.size(); ++i) {
            uint64_t p = primes[i];
            if (p < 3) continue; // Skip p=2

            uint64_t p_power_m = safe_power(p, exponent_m);
            if (p_power_m == 0) continue; // Overflow computing p^m or p^m > N initially
            if (p_power_m > N) continue;  // If p^m itself > N, no need for powers of 2

            uint64_t k_power_of_2 = 1; // Represents 2^k
            while(true) {
                 // Calculate n = p^m * 2^k
                 unsigned __int128 n_128 = static_cast<unsigned __int128>(p_power_m) * k_power_of_2;

                 if (n_128 == 0 || n_128 > N) {
                     break; // Overflow or exceeded limit N
                 }
                 numbers_set.insert(static_cast<uint64_t>(n_128));

                 // Prepare for next power of 2
                 unsigned __int128 next_k_power = (unsigned __int128)k_power_of_2 * 2;
                 if (next_k_power > numeric_limits<uint64_t>::max()) break; // Next power of 2 overflows uint64_t
                 k_power_of_2 = (uint64_t)next_k_power;
                 if (k_power_of_2 == 0) break;
            }
        }
    } else if (choice == 4) {
        // Numbers of the form n = 2^k * p^m * q^n, where q > p (both odd primes)
        for (size_t i = 0; i < primes.size(); ++i) {
            uint64_t p = primes[i];
            if (p < 3) continue; // p must be an odd prime

            uint64_t p_power_m = safe_power(p, exponent_m);
            if (p_power_m == 0) continue; // Overflow computing p^m
            if (p_power_m > N) break;     // Optimization: if p^m > N, no larger p will work either

            for (size_t j = i + 1; j < primes.size(); ++j) { // Iterate q > p
                uint64_t q = primes[j];
                // q is already > p and odd (since primes list > 2 only contains odd primes)

                uint64_t q_power_n = safe_power(q, exponent_n);
                if (q_power_n == 0) continue; // Overflow computing q^n

                // Calculate base = p^m * q^n
                unsigned __int128 base_128 = static_cast<unsigned __int128>(p_power_m) * q_power_n;
                if (base_128 == 0 || base_128 > N) {
                     // If base > N, increasing q further won't help for this p
                     // If base calculation overflows, increasing q further will also overflow.
                     break;
                }
                uint64_t base = static_cast<uint64_t>(base_128);

                // Now multiply by powers of 2 (2^k)
                uint64_t k_power_of_2 = 1; // Represents 2^k
                while (true) {
                    // Calculate n = base * 2^k
                    unsigned __int128 n_128 = static_cast<unsigned __int128>(base) * k_power_of_2;

                    if (n_128 == 0 || n_128 > N) {
                        break; // Overflow or exceeded limit N
                    }
                    numbers_set.insert(static_cast<uint64_t>(n_128));

                    // Prepare for next power of 2
                    unsigned __int128 next_k_power = (unsigned __int128)k_power_of_2 * 2;
                    if (next_k_power > numeric_limits<uint64_t>::max()) break; // Next power of 2 overflows uint64_t
                    k_power_of_2 = (uint64_t)next_k_power;
                    if (k_power_of_2 == 0) break;
                }
            }
        }
    }

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

    cout << "| " << left << setw(width_n - 1) << n
         << "| " << left << setw(width_sigma - 1) << sigma
         << "| " << left << setw(width_2n - 1) << two_n
         << "| " << left << setw(width_diff - 1) << diff
         << "| " << left << setw(width_combinations - 1) << combinations
         << "|\n";
}

// Function to check if a string contains only digits
bool is_digits(const std::string &str) {
    return str.find_first_not_of("0123456789") == std::string::npos;
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
    cout << "2. Numbers of the form n = 2^k * p (where p is an odd prime)\n";
    cout << "3. Numbers of the form n = 2^k * p^m (where p is an odd prime, m >= 1)\n";
    cout << "4. Numbers of the form n = 2^k * p^m * q^n (where p, q are odd primes, q > p, m,n >= 1)\n";
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
        } catch(...) {
            cout << "Invalid input. Defaulting to option 2 (n = 2^k * p).\n";
            choice = 2;
        }
    }

    // If choice 3 or 4 is selected, prompt for exponents
    if(choice == 3) {
        cout << "Enter the exponent m for p (form n = 2^k * p^m, m >= 1, default 1): ";
        string exponentMStr;
        getline(cin, exponentMStr);
        if(exponentMStr.empty()) {
            exponent_m = 1; // Default
        } else {
            try {
                exponent_m = stoi(exponentMStr);
                if(exponent_m < 1) {
                    cout << "Invalid exponent m. Must be >= 1. Defaulting to 1.\n";
                    exponent_m = 1;
                }
            } catch(...) {
                cout << "Invalid input for m. Defaulting to 1.\n";
                exponent_m = 1;
            }
        }
    } else if (choice == 4) {
         cout << "Enter the exponent m for p (form n = 2^k * p^m * q^n, m >= 1, default 1): ";
         string exponentMStr;
         getline(cin, exponentMStr);
         if(exponentMStr.empty()) {
             exponent_m = 1; // Default
         } else {
             try {
                 exponent_m = stoi(exponentMStr);
                 if(exponent_m < 1) {
                     cout << "Invalid exponent m. Must be >= 1. Defaulting to 1.\n";
                     exponent_m = 1;
                 }
             } catch(...) {
                 cout << "Invalid input for m. Defaulting to 1.\n";
                 exponent_m = 1;
             }
         }

         cout << "Enter the exponent n for q (form n = 2^k * p^m * q^n, n >= 1, default 1): ";
         string exponentNStr;
         getline(cin, exponentNStr);
         if(exponentNStr.empty()) {
             exponent_n = 1; // Default
         } else {
             try {
                 exponent_n = stoi(exponentNStr);
                 if(exponent_n < 1) {
                     cout << "Invalid exponent n. Must be >= 1. Defaulting to 1.\n";
                     exponent_n = 1;
                 }
             } catch(...) {
                 cout << "Invalid input for n. Defaulting to 1.\n";
                 exponent_n = 1;
             }
         }
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
        } catch(...) {
            cout << "Invalid input. Defaulting to All combinations.\n";
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
    uint64_t sqrtN = static_cast<uint64_t>(sqrt(static_cast<long double>(N))) + 1;
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
        case 2: cout << " of the form n = 2^k * p:"; break;
        case 3: cout << " of the form n = 2^k * p^" << exponent_m << ":"; break;
        case 4: cout << " of the form n = 2^k * p^" << exponent_m << " * q^" << exponent_n << " (q > p):"; break;
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
    const uint64_t report_interval = std::max((uint64_t)1, (uint64_t)(totalNumbers / 100)); // Report roughly 100 times

    for(uint64_t n : numbersToProcess) {
        processed++;
        // Display progress
        if (processed % report_interval == 0 || processed == totalNumbers) {
             double percent = (double)processed / totalNumbers * 100.0;
             cout << "\rProcessing: " << processed << " / " << totalNumbers << " (" << fixed << setprecision(1) << percent << "%)" << flush;
        }

        // Skip 0 or 1 (sigma(0)=0, sigma(1)=1)
        if (n < 2) continue;

        // Skip prime numbers if exclusion is enabled
        // Note: isPrime can be slow for very large n if basePrimes is relatively small
        if(excludePrimes && isPrime(n, basePrimes)) {
            continue;
        }

        // Compute sigma(n)
        uint64_t sigma = computeSigma(n, basePrimes);
        if (sigma == 0) { // Check if computeSigma overflowed
             // cerr << "\nWarning: Sigma calculation overflowed for n=" << n << ". Skipping." << endl;
             continue;
        }

        // Calculate 2n safely
        unsigned __int128 two_n_128 = (unsigned __int128)n * 2;
        if (two_n_128 > numeric_limits<uint64_t>::max()) {
             // cerr << "\nWarning: 2n calculation overflowed for n=" << n << ". Skipping." << endl;
             continue;
        }
        uint64_t two_n = (uint64_t)two_n_128;

        // Calculate difference: sigma(n) - 2n
        long long diff;
        // Check for potential overflow/underflow before subtraction
        if (sigma >= two_n) {
            diff = static_cast<long long>(sigma - two_n);
        } else {
            // sigma < two_n, result is negative
            // Ensure the difference fits within long long limits
            if (two_n - sigma > static_cast<uint64_t>(numeric_limits<long long>::max()) + 1ULL ) {
                 //cerr << "\nWarning: Negative difference calculation overflowed for n=" << n << ". Skipping." << endl;
                 continue; // Difference too large negative
            }
            diff = -static_cast<long long>(two_n - sigma);
        }


        // Compute divisors for finding combinations
        // This can be the bottleneck for large n
        vector<uint64_t> divisors = computeDivisors(n, basePrimes);
        if (divisors.empty() && n > 0) { // Check if computeDivisors indicated overflow
             cerr << "\nWarning: Divisor computation failed (likely overflow) for n=" << n << ". Skipping." << endl;
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
                    size_t second_sign = pairStr.find_first_of("+-", first_space + 1);
                    if (first_space != string::npos && second_sign != string::npos) {
                        string d1_abs_str = pairStr.substr(1, first_space - 1); // number after first sign
                        string d2_abs_str = pairStr.substr(second_sign + 1); // number after second sign
                        if (d1_abs_str != d2_abs_str) {
                             hasDistinctPair = true;
                             break; // Found one distinct pair, no need to check further
                        }
                    } else {
                        // Parsing failed? Log error or assume distinct to be safe?
                        // Let's assume it might be distinct if parsing is weird.
                        hasDistinctPair = true;
                        break;
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
    }
    cout << "\nFinished processing." << endl; // End progress line

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
            case 4: cout << " of the form n = 2^k * p^" << exponent_m << " * q^" << exponent_n << " (q > p)."; break;
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