#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <fstream>
#include <cctype> // For tolower
using namespace std;

/**
 * Function: sieveOfEratosthenes
 * -----------------------------
 * Implements the Sieve of Eratosthenes to identify prime numbers up to N.
 *
 * Parameters:
 *  - N: The upper limit up to which primes are identified.
 *
 * Returns:
 *  - A vector<bool> where isPrime[i] is true if i is a prime, false otherwise.
 */
vector<bool> sieveOfEratosthenes(int N) {
    vector<bool> isPrime(N + 1, true);
    if (N >= 0) isPrime[0] = isPrime[1] = false;
    for (int i = 2; i <= N; ++i) {
        if (isPrime[i]) {
            for (int j = 2 * i; j <= N; j += i) {
                isPrime[j] = false;
            }
        }
    }
    return isPrime;
}

/**
 * Function: isTwoPowerTimesPrime
 * ------------------------------
 * Checks if a given number n can be expressed in the form n = 2^k * p,
 * where p is an odd prime and k >= 0.
 *
 * Parameters:
 *  - n: The number to be checked.
 *  - isPrime: A vector indicating primality of numbers up to N.
 *
 * Returns:
 *  - The prime p if n = 2^k * p; otherwise, 0.
 */
int isTwoPowerTimesPrime(int n, const vector<bool> &isPrime) {
    if (n < 1) return 0;
    while (n % 2 == 0) {
        n /= 2;
    }
    if (n >= 3 && isPrime[n]) {
        return n;
    }
    return 0;
}

/**
 * Function: isTwoPowerTimesPrimePower
 * -----------------------------------
 * Checks if a given number n can be expressed in the form n = 2^k * p^m,
 * where p is an odd prime, m >= 1, and k >= 0.
 *
 * Parameters:
 *  - n: The number to be checked.
 *  - isPrime: A vector indicating primality of numbers up to N.
 *  - exponent: The exponent m for the prime p.
 *
 * Returns:
 *  - The prime p if n = 2^k * p^m; otherwise, 0.
 */
int isTwoPowerTimesPrimePower(int n, const vector<bool> &isPrime, int exponent) {
    if (n < 1 || exponent < 1) return 0;
    while (n % 2 == 0) {
        n /= 2;
    }
    // Now, n should be p^m
    double root = pow(static_cast<double>(n), 1.0 / exponent);
    int p = round(root);
    // Verify p^exponent == n and p is prime
    long long p_power = 1;
    for(int i = 0; i < exponent; ++i) {
        p_power *= p;
        if(p_power > n) break;
    }
    if (p_power == n && p >= 3 && isPrime[p]) {
        return p;
    }
    return 0;
}

/**
 * Function: computeDivisorsAndSigmaSingle
 * ---------------------------------------
 * Computes the divisors of a single number and their sum (sigma).
 *
 * Parameters:
 *  - n: The number for which divisors and sigma are computed.
 *
 * Returns:
 *  - A pair containing:
 *      - A vector of divisors of n.
 *      - The sum of divisors (sigma(n)).
 */
pair<vector<int>, long long> computeDivisorsAndSigmaSingle(int n) {
    vector<int> divisors;
    long long sigma = 0;
    for (int i = 1; i <= static_cast<int>(sqrt(n)); ++i) {
        if (n % i == 0) {
            divisors.push_back(i);
            sigma += i;
            if (i != n / i) {
                divisors.push_back(n / i);
                sigma += n / i;
            }
        }
    }
    sort(divisors.begin(), divisors.end());
    return {divisors, sigma};
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
 *  - D: Reference to the vector containing all divisors of n.
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
vector<string> findValidCombinations(const vector<int> &D, long long diff, SignOption signOption) {
    vector<string> validPairs;
    // Iterate through all possible pairs (d1, d2)
    for (auto it1 = D.begin(); it1 != D.end(); ++it1) {
        for (auto it2 = D.begin(); it2 != D.end(); ++it2) {
            int d1 = *it1;
            int d2 = *it2;

            // Skip if any of the divisors is zero
            if (d1 == 0 || d2 == 0) {
                continue;
            }

            // Depending on the signOption, include or exclude certain combinations
            switch(signOption) {
                case ALL:
                    // Check all four possible sign combinations
                    // ++
                    if ((d1 + d2) == diff) {
                        validPairs.emplace_back("+" + to_string(d1) + " +" + to_string(d2));
                    }
                    // +-
                    if ((d1 - d2) == diff) {
                        validPairs.emplace_back("+" + to_string(d1) + " -" + to_string(d2));
                    }
                    // -+
                    if ((-d1 + d2) == diff) {
                        validPairs.emplace_back("-" + to_string(d1) + " +" + to_string(d2));
                    }
                    // --
                    if ((-d1 - d2) == diff) {
                        validPairs.emplace_back("-" + to_string(d1) + " -" + to_string(d2));
                    }
                    break;

                case MIXED_SIGNS:
                    // Include only +-, -+
                    // +-
                    if ((d1 - d2) == diff) {
                        validPairs.emplace_back("+" + to_string(d1) + " -" + to_string(d2));
                    }
                    // -+
                    if ((-d1 + d2) == diff) {
                        validPairs.emplace_back("-" + to_string(d1) + " +" + to_string(d2));
                    }
                    break;

                case ONLY_POSITIVE:
                    // Include only ++
                    if ((d1 + d2) == diff) {
                        validPairs.emplace_back("+" + to_string(d1) + " +" + to_string(d2));
                    }
                    break;

                case ONLY_NEGATIVE:
                    // Include only --
                    if ((-d1 - d2) == diff) {
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
 *  - isPrime: A vector indicating primality of numbers up to N.
 *  - exponent: The exponent m for the prime p (only relevant for choice 3).
 *
 * Returns:
 *  - A vector of integers representing the numbers to process.
 */
vector<int> generateNumbers(int N, int choice, const vector<bool> &isPrime, int exponent = 2) {
    vector<int> numbers;

    if (choice == 1) {
        // All numbers from 1 to N
        for (int n = 1; n <= N; ++n) {
            numbers.push_back(n);
        }
    }
    else if (choice == 2) {
        // Numbers of the form n = 2^k * p
        for (int p = 3; p <= N; ++p) {
            if (isPrime[p]) {
                int power = 1;
                while (true) {
                    long long n = static_cast<long long>(power) * p;
                    if (n > N) break;
                    numbers.push_back(static_cast<int>(n));
                    power *= 2;
                }
            }
        }
    }
    else if (choice == 3) {
        // Numbers of the form n = 2^k * p^m
        for (int p = 3; p <= N; ++p) {
            if (isPrime[p]) {
                // Compute p^m
                long long p_power = 1;
                for(int i = 0; i < exponent; ++i) {
                    p_power *= p;
                    if(p_power > N) break;
                }
                if(p_power > N) continue;

                int power = 1;
                while (true) {
                    long long n = p_power * power;
                    if (n > N) break;
                    numbers.push_back(static_cast<int>(n));
                    power *= 2;
                }
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
    const int width_n = 10;
    const int width_sigma = 15;
    const int width_2n = 15;
    const int width_diff = 15;
    const int width_combinations = 35;

    // Print the top border
    cout << "+" << string(width_n, '-') << "+"
         << string(width_sigma, '-') << "+"
         << string(width_2n, '-') << "+"
         << string(width_diff, '-') << "+"
         << string(width_combinations, '-') << "+\n";

    // Print the header row
    cout << "| " << left << setw(width_n - 1) << "n"
         << "| " << left << setw(width_sigma - 1) << "sigma(n)"
         << "| " << left << setw(width_2n - 1) << "2n"
         << "| " << left << setw(width_diff - 1) << "diff"
         << "| " << left << setw(width_combinations - 1) << "Valid (d1, d2) Combinations"
         << "|\n";

    // Print the separator
    cout << "+" << string(width_n, '=') << "+"
         << string(width_sigma, '=') << "+"
         << string(width_2n, '=') << "+"
         << string(width_diff, '=') << "+"
         << string(width_combinations, '=') << "+\n";
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
void printTableRow(int n, long long sigma, long long two_n, long long diff, const string& combinations) {
    // Define column widths
    const int width_n = 10;
    const int width_sigma = 15;
    const int width_2n = 15;
    const int width_diff = 15;
    const int width_combinations = 35;

    cout << "| " << left << setw(width_n - 1) << n
         << "| " << left << setw(width_sigma - 1) << sigma
         << "| " << left << setw(width_2n - 1) << two_n
         << "| " << left << setw(width_diff - 1) << diff
         << "| " << left << setw(width_combinations - 1) << combinations
         << "|\n";
}

int main() {
    int N;
    cout << "Enter the upper limit N: ";
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

    // Generate primes up to N
    vector<bool> isPrime = sieveOfEratosthenes(N);

    // Generate numbers to process based on user's choice
    vector<int> numbersToProcess = generateNumbers(N, choice, isPrime, (choice == 3) ? exponent : 2);

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

    for(auto n : numbersToProcess) {
        // Skip prime numbers if exclusion is enabled
        if(excludePrimes && isPrime[n]) {
            continue;
        }

        // Compute sigma(n)
        auto [divisors, sigma] = computeDivisorsAndSigmaSingle(n);
        long long two_n = 2LL * n;
        long long diff = sigma - two_n;

        // Find valid (d1, d2) combinations based on signOption
        vector<string> validPairs = findValidCombinations(divisors, diff, signOption);

        // **New Logic: Exclude numbers where all validPairs involve the same divisor**
        bool hasDistinctPair = false;
        string combinations; // Initialize here to avoid scope issues
        for(const auto &pairStr : validPairs) {
            // Parse the pair string to extract d1 and d2
            // Possible formats:
            // "+d1 +d2", "+d1 -d2", "-d1 +d2", "-d1 -d2"
            size_t spacePos = pairStr.find(' ');
            if(spacePos != string::npos) {
                // Extract d1
                string d1Str = pairStr.substr(1, spacePos - 1); // Skip the '+' or '-'

                // Extract d2 with sign intact
                string d2Str = pairStr.substr(spacePos + 1); // Includes the sign

                // No need to remove leading '+', as stoi handles it
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
            string combinationsCombined;
            for(const auto &combination : validPairs) {
                if(!combinationsCombined.empty()) combinationsCombined += "; ";
                combinationsCombined += combination;
            }

            // Print the row in the table
            printTableRow(n, sigma, two_n, diff, combinationsCombined);

            // Write the row to the CSV file
            // Enclose combinations in quotes to handle commas or special characters
            csvFile << n << "," << sigma << "," << two_n << "," << diff << ",\"" << combinationsCombined << "\"\n";
        }
    }

    // Print the bottom border of the table
    cout << "+" << string(10, '-') << "+"
         << string(15, '-') << "+"
         << string(15, '-') << "+"
         << string(15, '-') << "+"
         << string(35, '-') << "+\n";

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