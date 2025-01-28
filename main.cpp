#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

using namespace std;

// we love AI for making descriptiong :)
/**
 * Function: computeDivisorsAndSigma
 * --------------------
 * Precomputes the list of divisors and the sum of divisors (sigma) for each number up to N.
 *
 * Parameters:
 * - N: The upper limit up to which divisors and sigma values are computed.
 * - divisors: A reference to a vector of vectors where divisors[i] will contain all divisors of i.
 * - sigma: A reference to a vector where sigma[i] will store the sum of divisors of i.
 *
 * Returns:
 * - void (the results are stored in the passed vectors).
 */
void computeDivisorsAndSigma(int N, vector<vector<int>> &divisors, vector<long long> &sigma) {
    // Initialize the divisors and sigma vectors
    divisors.assign(N + 1, vector<int>());
    sigma.assign(N + 1, 0);

    // Sieve-like approach to populate divisors and compute sigma
    for(int i = 1; i <= N; ++i){
        for(int j = i; j <= N; j += i){
            divisors[j].push_back(i); // i is a divisor of j
            sigma[j] += i;            // Add i to the sum of divisors of j
        }
    }
}

/**
 * Function: sieveOfEratosthenes
 * --------------------
 * Implements the Sieve of Eratosthenes to identify prime numbers up to N.
 *
 * Parameters:
 * - N: The upper limit up to which primes are identified.
 *
 * Returns:
 * - A vector<bool> where isPrime[i] is true if i is a prime, false otherwise.
 */
vector<bool> sieveOfEratosthenes(int N) {
    vector<bool> isPrime(N + 1, true);
    if(N >= 0) isPrime[0] = isPrime[1] = false;
    for(int i = 2; i <= N; ++i){
        if(isPrime[i]){
            for(int j = 2 * i; j <= N; j += i){
                isPrime[j] = false;
            }
        }
    }
    return isPrime;
}

/**
 * Function: isTwoPowerTimesPrime
 * --------------------
 * Checks if a given number n can be expressed in the form n = 2^k * p,
 * where p is an odd prime and k >= 0.
 *
 * Parameters:
 * - n: The number to be checked.
 * - isPrime: A vector indicating primality of numbers up to N.
 *
 * Returns:
 * - The prime p if n = 2^k * p; otherwise, 0.
 */
int isTwoPowerTimesPrime(int n, const vector<bool> &isPrime){
    if(n < 1) return 0;
    while(n % 2 == 0){
        n /= 2;
    }
    if(n >= 3 && isPrime[n]){
        return n;
    }
    return 0;
}

/**
 * Function: isTwoPowerTimesPrimeSquared
 * --------------------
 * Checks if a given number n can be expressed in the form n = 2^k * p^2,
 * where p is an odd prime and k >= 0.
 *
 * Parameters:
 * - n: The number to be checked.
 * - isPrime: A vector indicating primality of numbers up to N.
 *
 * Returns:
 * - The prime p if n = 2^k * p^2; otherwise, 0.
 */
int isTwoPowerTimesPrimeSquared(int n, const vector<bool> &isPrime){
    if(n < 1) return 0;
    while(n % 2 == 0){
        n /= 2;
    }
    // Now, n should be p^2
    double sqrt_n = sqrt((double)n);
    int p = round(sqrt_n);
    if(p * p == n && p >= 3 && isPrime[p]){
        return p;
    }
    return 0;
}

/**
 * Function: findValidCombinations
 * --------------------
 * Finds all valid combinations of divisors d1 and d2 with their respective signs
 * that satisfy the equation: sigma(n) = 2n ± d1 ± d2
 *
 * Parameters:
 * - D: Reference to the vector containing all divisors of n.
 * - diff: The difference sigma(n) - 2n.
 *
 * Returns:
 * - A vector of strings, each representing a valid combination of d1 and d2 with signs.
 */
vector<string> findValidCombinations(const vector<int> &D, long long diff){
    vector<string> validPairs;
    // Iterate through all possible pairs (d1, d2)
    for(auto it1 = D.begin(); it1 != D.end(); ++it1){
        for(auto it2 = D.begin(); it2 != D.end(); ++it2){
            int d1 = *it1;
            int d2 = *it2;

            // Skip if any of the divisors is zero
            if(d1 == 0 || d2 == 0){
                continue;
            }

            // Check all four possible sign combinations
            if( (d1 + d2) == diff ){
                validPairs.emplace_back("+" + to_string(d1) + " +" + to_string(d2));
            }
            if( (d1 - d2) == diff ){
                validPairs.emplace_back("+" + to_string(d1) + " -" + to_string(d2));
            }
            if( (d2 - d1) == diff ){
                validPairs.emplace_back("+" + to_string(d2) + " -" + to_string(d1));
            }
            if( (-d1 - d2) == diff ){
                validPairs.emplace_back("-" + to_string(d1) + " -" + to_string(d2));
            }
        }
    }

    // Remove duplicate combinations
    sort(validPairs.begin(), validPairs.end());
    validPairs.erase(unique(validPairs.begin(), validPairs.end()), validPairs.end());

    return validPairs;
}

int main(){
    int N;
    cout << "Enter the upper limit N: ";
    cin >> N;

    if(N < 1){
        cout << "Please enter a positive integer greater than 0.\n";
        return 1;
    }

    // Precompute primes up to N using Sieve of Eratosthenes
    vector<bool> isPrime = sieveOfEratosthenes(N);

    // Precompute divisors and sigma(n) for all numbers up to N
    vector<vector<int>> divisors;
    vector<long long> sigma;
    computeDivisorsAndSigma(N, divisors, sigma);

    // Prompt user for form selection
    cout << "\nChoose the form to search for near perfect numbers:\n";
    cout << "1. All numbers (default)\n";
    cout << "2. Numbers of the form n = 2^k * p (where p is an odd prime)\n";
    cout << "3. Numbers of the form n = 2^k * p^2 (where p is an odd prime)\n";
    cout << "Enter 1, 2, or 3 (Press Enter for default - All numbers): ";

    // Flush the input buffer
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    string choiceStr;
    getline(cin, choiceStr);
    int choice;
    if(choiceStr.empty()){
        choice = 1; // Default to all numbers
    }
    else{
        choice = stoi(choiceStr);
        if(choice <1 || choice >3){
            cout << "Invalid choice. Defaulting to all numbers.\n";
            choice =1;
        }
    }

    // Header for output
    cout << "\nNear Perfect Numbers up to " << N << ":\n";
    cout << "----------------------------------------\n";

    bool anyNearPerfect = false;

    // Iterate through each number to check the near perfect condition based on the selected form
    for(int n = 1; n <= N; ++n){
        bool process = false;
        switch(choice){
            case 1:
                process = true; // Process all numbers
                break;
            case 2:
                if(isTwoPowerTimesPrime(n, isPrime)){
                    process = true;
                }
                break;
            case 3:
                if(isTwoPowerTimesPrimeSquared(n, isPrime)){
                    process = true;
                }
                break;
            default:
                process = true;
                break;
        }

        if(process){
            long long s = sigma[n];
            long long two_n = 2LL * n;
            long long diff = s - two_n;

            const vector<int> &D = divisors[n];

            // Find all valid (d1, d2) combinations that satisfy sigma(n) = 2n ± d1 ± d2
            vector<string> validPairs = findValidCombinations(D, diff);

            if(!validPairs.empty()){
                anyNearPerfect = true;
                cout << "n = " << n << "\n";
                cout << "sigma(n) = " << s << "\n";
                cout << "2n      = " << two_n << "\n";
                cout << "diff    = " << diff << " (sigma(n) - 2n)\n";
                cout << "Valid combinations of d1 and d2:\n";

                for(const auto &combination : validPairs){
                    cout << "  " << combination << "\n";
                }
                cout << "----------------------------------------\n";
            }
        }
    }

    if(!anyNearPerfect){
        cout << "No near perfect numbers found up to " << N;
        switch(choice){
            case 1:
                cout << ".\n";
                break;
            case 2:
                cout << " in the form n = 2^k * p.\n";
                break;
            case 3:
                cout << " in the form n = 2^k * p^2.\n";
                break;
            default:
                cout << ".\n";
                break;
        }
    }

    return 0;
}