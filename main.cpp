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

using namespace std;

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
    primes.reserve(static_cast<size_t>(sqrtN / log(sqrtN)) + 1); // Approximation
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
 *
 * Parameters:
 *  - n: The number to check for primality.
 *  - primes: A vector of known primes up to sqrt(n).
 *
 * Returns:
 *  - True if n is prime, False otherwise.
 */
bool isPrime(uint64_t n, const vector<uint64_t>& primes) {
    if (n < 2) return false;
    uint64_t sqrt_n = static_cast<uint64_t>(sqrt(static_cast<long double>(n)));
    for (const auto& p : primes) {
        if (p > sqrt_n) break;
        if (n % p == 0) return false;
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
 *  - primes: A vector of primes up to sqrt(n).
 *
 * Returns:
 *  - The sum of divisors of n.
 */
uint64_t computeSigma(uint64_t n, const vector<uint64_t>& primes) {
    if (n == 0) return 0;

    uint64_t sigma = 1;
    uint64_t original_n = n;

    for (const auto& p : primes) {
        if (p * p > n) break;
        if (n % p == 0) {
            uint64_t term = 1;
            uint64_t power = 1;
            while (n % p == 0) {
                power *= p;
                term += power;
                n /= p;
            }
            sigma *= term;
        }
    }
    if (n > 1) {
        sigma *= (1 + n);
    }
    return sigma;
}

/**
 * Function: isTwoPowerTimesPrime
 * ------------------------------
 * Checks if a given number n can be expressed in the form n = 2^k * p,
 * where p is an odd prime and k >= 0.
 *
 * Parameters:
 *  - n: The number to be checked.
 *  - primes: Vector containing all primes up to sqrt(n).
 *
 * Returns:
 *  - The prime p if n = 2^k * p; otherwise, 0.
 */
uint64_t isTwoPowerTimesPrime(uint64_t n, const vector<uint64_t>& primes) {
    if (n < 1) return 0;
    while (n % 2 == 0) {
        n /= 2;
    }
    if (n >= 3 && isPrime(n, primes)) {
        return n;
    }
    return 0;
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
 *  - D: Reference to the vector containing all prime factors of n.
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

    // Iterate through all possible pairs (d1, d2)
    for (auto it1 = D.begin(); it1 != D.end(); ++it1) {
        for (auto it2 = D.begin(); it2 != D.end(); ++it2) {
            uint64_t d1 = *it1;
            uint64_t d2 = *it2;

            // Skip if any of the divisors is zero
            if (d1 == 0 || d2 == 0) {
                continue;
            }

            // Depending on the signOption, include or exclude certain combinations
            switch(signOption) {
                case ALL:
                    // ++
                    if ((static_cast<long long>(d1) + static_cast<long long>(d2)) == diff) {
                        validPairs.emplace_back("+" + to_string(d1) + " +" + to_string(d2));
                    }
                    // +-
                    if ((static_cast<long long>(d1) - static_cast<long long>(d2)) == diff) {
                        validPairs.emplace_back("+" + to_string(d1) + " -" + to_string(d2));
                    }
                    // -+
                    if ((-static_cast<long long>(d1) + static_cast<long long>(d2)) == diff) {
                        validPairs.emplace_back("-" + to_string(d1) + " +" + to_string(d2));
                    }
                    // --
                    if ((-static_cast<long long>(d1) - static_cast<long long>(d2)) == diff) {
                        validPairs.emplace_back("-" + to_string(d1) + " -" + to_string(d2));
                    }
                    break;

                case MIXED_SIGNS:
                    // +-
                    if ((static_cast<long long>(d1) - static_cast<long long>(d2)) == diff) {
                        validPairs.emplace_back("+" + to_string(d1) + " -" + to_string(d2));
                    }
                    // -+
                    if ((-static_cast<long long>(d1) + static_cast<long long>(d2)) == diff) {
                        validPairs.emplace_back("-" + to_string(d1) + " +" + to_string(d2));
                    }
                    break;

                case ONLY_POSITIVE:
                    // ++
                    if ((static_cast<long long>(d1) + static_cast<long long>(d2)) == diff) {
                        validPairs.emplace_back("+" + to_string(d1) + " +" + to_string(d2));
                    }
                    break;

                case ONLY_NEGATIVE:
                    // --
                    if ((-static_cast<long long>(d1) - static_cast<long long>(d2)) == diff) {
                        validPairs.emplace_back("-" + to_string(d1) + " -" + to_string(d2));
                    }
                    break;

                default:
                    // Do nothing for undefined options
                    break;
            }
        }
    }

    // Remove duplicate combinations
    sort(validPairs.begin(), validPairs.end());
    validPairs.erase(unique(validPairs.begin(), validPairs.end()), validPairs.end());

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
 *  - primes: Vector containing all primes up to sqrt(N).
 *  - exponent: The exponent m for the prime p (only relevant for choice 3).
 *
 * Returns:
 *  - A vector of integers representing the numbers to process.
 */
vector<uint64_t> generateNumbers(uint64_t N, int choice, const vector<uint64_t> &primes, int exponent = 2) {
    vector<uint64_t> numbers;

    if (choice == 1) {
        // All numbers from 1 to N - Note: Processing all numbers up to 1e12 is impractical.
        // Consider limiting N or optimizing further if choice 1 is necessary.
        for (uint64_t n = 1; n <= N; ++n) {
            numbers.push_back(n);
            if (n == numeric_limits<uint64_t>::max()) break; // Prevent overflow
        }
    }
    else if (choice == 2) {
        // Numbers of the form n = 2^k * p
        for (const auto& p : primes) {
            if (p < 3) continue; // Skip even primes
            uint64_t power = 1;
            while (true) {
                unsigned __int128 n = static_cast<unsigned __int128>(power) * p;
                if (n > N) break;
                numbers.emplace_back(static_cast<uint64_t>(n));
                power *= 2;
                if (power == 0) break; // Prevent overflow
            }
        }
    }
    else if (choice == 3) {
        // Numbers of the form n = 2^k * p^m
        for (const auto& p : primes) {
            if (p < 3) continue; // Skip even primes
            unsigned __int128 p_power = 1;
            for(int i = 0; i < exponent; ++i) {
                p_power *= p;
                if(p_power > N) break;
            }
            if(p_power > N) continue;

            uint64_t power = 1;
            while (true) {
                unsigned __int128 n = p_power * power;
                if (n > N) break;
                numbers.emplace_back(static_cast<uint64_t>(n));
                power *= 2;
                if (power == 0) break; // Prevent overflow
            }
        }
    }

    // Remove duplicates and sort
    sort(numbers.begin(), numbers.end());
    numbers.erase(unique(numbers.begin(), numbers.end()), numbers.end());

    return numbers;
}

/**
 * Function: printTableHeader
 * --------------------------
 * Prints the header of the formatted table with vertical splitters.
 */
void printTableHeader() {
    // Define column widths
    const int width_n = 12;
    const int width_sigma = 20;
    const int width_2n = 20;
    const int width_diff = 20;
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
    // Define column widths
    const int width_n = 12;
    const int width_sigma = 20;
    const int width_2n = 20;
    const int width_diff = 20;
    const int width_combinations = 50;

    cout << "| " << left << setw(width_n - 1) << n
         << "| " << left << setw(width_sigma - 1) << sigma
         << "| " << left << setw(width_2n - 1) << two_n
         << "| " << left << setw(width_diff - 1) << diff
         << "| " << left << setw(width_combinations - 1) << combinations
         << "|\n";
}

int main() {
    uint64_t N;
    cout << "Enter the upper limit N (up to 1,000,000,000,000): ";
    while(!(cin >> N) || N < 1) {
        cout << "Please enter a positive integer greater than 0: ";
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    // Clear the input buffer
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    cout << "\nChoose the form to search for near perfect numbers:\n";
    cout << "1. All numbers (default)\n";
    cout << "2. Numbers of the form n = 2^k * p (where p is an odd prime)\n";
    cout << "3. Numbers of the form n = 2^k * p^m (where p is an odd prime and m is user-specified)\n";
    cout << "Enter 1, 2, or 3 (Press Enter for default - All numbers): ";

    string choiceStr;
    getline(cin, choiceStr);
    int choice;
    int exponent = 2; // Default exponent for choice 3
    if(choiceStr.empty()) {
        choice = 1;
    }
    else {
        try {
            choice = stoi(choiceStr);
            if(choice < 1 || choice > 3) {
                cout << "Invalid choice. Defaulting to all numbers.\n";
                choice = 1;
            }
        }
        catch(...) {
            cout << "Invalid input. Defaulting to all numbers.\n";
            choice = 1;
        }
    }

    // If choice 3 is selected, prompt for exponent
    if(choice == 3) {
        cout << "Enter the exponent m for the form n = 2^k * p^m (m must be a positive integer >=1): ";
        string exponentStr;
        getline(cin, exponentStr);
        if(exponentStr.empty()) {
            cout << "No exponent entered. Defaulting to m = 2.\n";
            exponent = 2;
        }
        else {
            try {
                exponent = stoi(exponentStr);
                if(exponent < 1) {
                    cout << "Invalid exponent. Defaulting to m = 2.\n";
                    exponent = 2;
                }
            }
            catch(...) {
                cout << "Invalid input. Defaulting to m = 2.\n";
                exponent = 2;
            }
        }
    }

    // Prompt to exclude primes
    cout << "Do you want to exclude prime numbers from the output? (y/n): ";
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
    cout << "\nSelect the sign combinations to include in the search:\n";
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
            }
            else {
                cout << "Invalid choice. Defaulting to All combinations.\n";
                signOption = ALL;
            }
        }
        catch(...) {
            cout << "Invalid input. Defaulting to All combinations.\n";
            signOption = ALL;
        }
    }

    // Generate base primes up to sqrt(N) using simple sieve
    uint64_t sqrtN = static_cast<uint64_t>(sqrt(static_cast<long double>(N))) + 1;
    cout << "\nGenerating base primes up to " << sqrtN << " using simple sieve...\n";
    vector<uint64_t> basePrimes = simpleSieve(sqrtN);
    cout << "Number of base primes found: " << basePrimes.size() << "\n";

    // Generate numbers to process based on user's choice
    cout << "Generating numbers based on the selected form...\n";
    vector<uint64_t> numbersToProcess = generateNumbers(N, choice, basePrimes, (choice == 3) ? exponent : 2);
    cout << "Number of numbers to process: " << numbersToProcess.size() << "\n";

    // Display the type of numbers being processed
    cout << "\nNear Perfect Numbers up to " << N;
    switch(choice) {
        case 1: cout << ":"; break;
        case 2: cout << " of the form n = 2^k * p:"; break;
        case 3: cout << " of the form n = 2^k * p^" << exponent << ":"; break;
    }
    if(excludePrimes) {
        cout << " (Primes excluded from output)";
    }

    // Display selected sign combinations
    cout << "\nSign Combinations Included: ";
    switch(signOption) {
        case ALL:
            cout << "All combinations (++ , +-, -+, --)";
            break;
        case MIXED_SIGNS:
            cout << "Only mixed signs (+-, -+)";
            break;
        case ONLY_POSITIVE:
            cout << "Only positive signs (++ )";
            break;
        case ONLY_NEGATIVE:
            cout << "Only negative signs (--)";
            break;
        default:
            cout << "All combinations (++ , +-, -+, --)";
            break;
    }
    cout << "\n";

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
    for(auto n : numbersToProcess) {
        processed++;
        // Display progress every 1 million
        if (processed % 1000000 == 0) {
            cout << "\rProcessed " << processed << " / " << totalNumbers << " numbers..." << flush;
        }

        // Skip prime numbers if exclusion is enabled
        if(excludePrimes && isPrime(n, basePrimes)) {
            continue;
        }

        // Compute sigma(n)
        uint64_t sigma = computeSigma(n, basePrimes);
        uint64_t two_n = 2ULL * n;
        long long diff = static_cast<long long>(sigma) - static_cast<long long>(two_n);

        // Find valid (d1, d2) combinations based on signOption
        // To find d1 and d2, you need the list of divisors
        // Instead of computing all divisors, you can derive them from prime factors
        // However, for simplicity, we'll compute the divisors here

        // Compute divisors
        // Since sigma(n) is already computed, but you need the actual divisors for combinations
        // To optimize, consider memoizing or using other strategies if this becomes a bottleneck
        // For now, we'll compute them on-the-fly
        // Note: For large n, this can be optimized further

        // Function to compute all divisors from prime factors
        vector<uint64_t> divisors = {1};
        uint64_t temp = n;
        for(const auto& p : basePrimes) {
            if(p * p > temp) break;
            if(temp % p == 0) {
                vector<uint64_t> tempDivisors;
                uint64_t power = 1;
                while(temp % p == 0) {
                    power *= p;
                    temp /= p;
                    // Multiply existing divisors with the new power
                    for(auto d : divisors) {
                        tempDivisors.push_back(d * power);
                    }
                }
                // Merge and remove duplicates
                divisors.insert(divisors.end(), tempDivisors.begin(), tempDivisors.end());
            }
        }
        if(temp > 1) {
            vector<uint64_t> tempDivisors;
            for(auto d : divisors) {
                tempDivisors.push_back(d * temp);
            }
            divisors.insert(divisors.end(), tempDivisors.begin(), tempDivisors.end());
        }
        sort(divisors.begin(), divisors.end());

        // Now find valid combinations
        vector<string> validPairs = findValidCombinations(divisors, diff, signOption);

        // **New Logic: Exclude numbers where all validPairs involve the same divisor**
        bool hasDistinctPair = false;
        string combinationsCombinedStr;
        for(const auto &pairStr : validPairs) {
            // Parse the pair string to extract d1 and d2
            size_t spacePos = pairStr.find(' ');
            if(spacePos != string::npos) {
                // Extract d1
                string d1Str = pairStr.substr(1, spacePos - 1); // Skip the '+' or '-'

                // Extract d2 with sign intact
                string d2Str = pairStr.substr(spacePos + 1); // Includes the sign

                try {
                    int d1 = stoi(d1Str);
                    int d2 = stoi(d2Str); // Corrected Line
                    // Check if d1 and d2 have different absolute values
                    if(abs(d1) != abs(d2)) {
                        hasDistinctPair = true;
                        break;
                    }
                }
                catch(const std::invalid_argument& e) {
                    cerr << "Invalid pair format: " << pairStr << ". Skipping this pair.\n";
                    continue; // Skip this pair and continue with the next
                }
                catch(const std::out_of_range& e) {
                    cerr << "Pair values out of range: " << pairStr << ". Skipping this pair.\n";
                    continue;
                }
            }
        }

        // Proceed only if there's at least one distinct pair
        if(!hasDistinctPair) {
            continue; // Skip this number as it doesn't meet the criteria
        }

        if(!validPairs.empty()) {
            anyNearPerfect = true;
            // Combine all valid pairs into a single string separated by "; "
            string combinationsCombinedStr;
            for(const auto &combination : validPairs) {
                if(!combinationsCombinedStr.empty()) combinationsCombinedStr += "; ";
                combinationsCombinedStr += combination;
            }

            // Print the row in the table
            printTableRow(n, sigma, two_n, diff, combinationsCombinedStr);

            // Write the row to the CSV file
            // Enclose combinations in quotes to handle commas or special characters
            csvFile << n << "," << sigma << "," << two_n << "," << diff << ",\"" << combinationsCombinedStr << "\"\n";
        }
    }
    cout << "\n"; // For clean newline after progress updates
    // Print the bottom border of the table
    cout << "+" << string(12, '-')
         << "+" << string(20, '-')
         << "+" << string(20, '-')
         << "+" << string(20, '-')
         << "+" << string(50, '-') << "+\n";

    // Close the CSV file
    csvFile.close();

    if(anyNearPerfect) {
        cout << "\nResults have been saved to 'near_perfect_numbers.csv'.\n";
    }
    else {
        cout << "\nNo near perfect numbers found";
        switch(choice) {
            case 1: cout << "."; break;
            case 2: cout << " of the form n = 2^k * p."; break;
            case 3: cout << " of the form n = 2^k * p^" << exponent << "."; break;
        }
        if(excludePrimes) {
            cout << " Primes were excluded from the output.";
        }
        cout << "\n";
    }

    return 0;
}