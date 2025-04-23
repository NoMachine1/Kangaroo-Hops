#include <gmp.h>
#include <gmpxx.h>
#include <chrono>
#include <ctime>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>
#include <array>
#include <set>
#include <sstream>
#include <omp.h>
#include <immintrin.h>

using namespace std;

typedef array<mpz_class, 2> Point;

const mpz_class modulo("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 16);
const mpz_class Gx("79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 16);
const mpz_class Gy("483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8", 16);
const Point PG = {Gx, Gy};
const Point Z = {0, 0};

auto starttime = chrono::high_resolution_clock::now();

// Point addition on the elliptic curve
Point add(const Point& P, const Point& Q, const mpz_class& p = modulo) {
    if (P == Z) return Q;
    if (Q == Z) return P;
    const mpz_class& P0 = P[0];
    const mpz_class& P1 = P[1];
    const mpz_class& Q0 = Q[0];
    const mpz_class& Q1 = Q[1];
    mpz_class lmbda, num, denom, inv;
    if (P != Q) {
        num = Q1 - P1;
        denom = Q0 - P0;
    } else {
        if (P1 == 0) return Z;
        num = 3 * P0 * P0;
        denom = 2 * P1;
    }
    mpz_invert(inv.get_mpz_t(), denom.get_mpz_t(), p.get_mpz_t());
    lmbda = (num * inv) % p;
    mpz_class x = (lmbda * lmbda - P0 - Q0) % p;
    if (x < 0) x += p;
    mpz_class y = (lmbda * (P0 - x) - P1) % p;
    if (y < 0) y += p;
    return {x, y};
}

// Scalar multiplication on the elliptic curve
Point mul(const mpz_class& k, const Point& P = PG, const mpz_class& p = modulo) {
    Point R0 = Z;
    Point R1 = P;

    // Determine the number of bits in k
    size_t bit_length = mpz_sizeinbase(k.get_mpz_t(), 2);

    #pragma omp parallel
    {
        // Each thread will work on its own copy of R0 and R1
        Point local_R0 = Z;
        Point local_R1 = P;

        #pragma omp for
        for (int i = bit_length - 1; i >= 0; --i) {
            // Check if the i-th bit is set
            if (mpz_tstbit(k.get_mpz_t(), i)) {
                local_R0 = add(local_R0, local_R1, p);
                local_R1 = add(local_R1, local_R1, p);
            } else {
                local_R1 = add(local_R0, local_R1, p);
                local_R0 = add(local_R0, local_R0, p);
            }
        }

        // Combine results from all threads
        #pragma omp critical
        {
            R0 = add(R0, local_R0, p);
        }
    }

    return R0;
}

// Point Subtraction
Point point_subtraction(const Point& P, const Point& Q, const mpz_class& p = modulo) {
    Point Q_neg = {Q[0], (-Q[1]) % p};
    return add(P, Q_neg, p);
}

// Compute Y from X and parity
mpz_class X2Y(const mpz_class& X, int y_parity, const mpz_class& p = modulo) {
    mpz_class X_cubed = (X * X * X) % p;
    mpz_class tmp = (X_cubed + mpz_class(7)) % p;
    mpz_class Y;
    mpz_class exp = (p + mpz_class(1)) / mpz_class(4);
    mpz_powm(Y.get_mpz_t(), tmp.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
    if ((Y % 2) != y_parity) {
        Y = p - Y;
    }
    return Y;
}

// Handle Solution
bool handle_solution(const mpz_class& solution) {
    string HEX = solution.get_str(16);
    mpz_class dec = solution;
    auto end = chrono::system_clock::now();
    time_t end_time = chrono::system_clock::to_time_t(end);
    cout << "\n\033[01;33m[+]\033[32m PUZZLE SOLVED: \033[32m" << ctime(&end_time) << "\r";
    cout << "\r\033[01;33m[+]\033[32m Private key (dec): \033[32m" << dec << "\033[0m" << endl;
    static std::ofstream file("KEYFOUNDKEYFOUND.txt", std::ios::app);
    file << "\n" << string(140, '-') << endl;
    file << "SOLVED " << ctime(&end_time);
    file << "Private Key (decimal): " << dec << endl;
    file << "Private Key (hex): " << dec.get_str(16) << endl;
    file << string(140, '-') << endl;
    return true;
}

// Generate powers of two for hopping
vector<mpz_class> generate_powers_of_two(int hop_modulo) {
    vector<mpz_class> powers(hop_modulo);
    #pragma omp parallel for
    for (int pw = 0; pw < hop_modulo; ++pw) {
        powers[pw] = mpz_class(1) << pw;
    }
    return powers;
}

// Main kangaroo function
// Main kangaroo function with fixed speed calculation
string kangaroo(const vector<Point>& P, const Point& W0, const Point& W1, const mpz_class& DP_rarity, int Nw, int Nt, int hop_modulo,
              const mpz_class& upper_range_limit, const mpz_class& lower_range_limit, const vector<mpz_class>& powers_of_two) {
    alignas(32) vector<Point> T(Nt, Z), W(Nw, Z), Wi(Nw, Z);
    alignas(32) vector<mpz_class> t(Nt), w(Nw), wi(Nw);
    gmp_randclass rand(gmp_randinit_default);

    // Initialize T, W, and Wi with random points
    #pragma omp parallel for
    for (int k = 0; k < Nt; ++k) {
        t[k] = lower_range_limit + rand.get_z_range(upper_range_limit - lower_range_limit);
        T[k] = mul(t[k], P[0], modulo);
    }

    #pragma omp parallel for
    for (int k = 0; k < Nw; ++k) {
        w[k] = rand.get_z_range(upper_range_limit - lower_range_limit);
        W[k] = add(W0, mul(w[k], P[0], modulo), modulo);
    }

    #pragma omp parallel for
    for (int k = 0; k < Nw; ++k) {
        wi[k] = rand.get_z_range(upper_range_limit - lower_range_limit);
        Wi[k] = add(W1, mul(wi[k], P[0], modulo), modulo);
    }

    long long Hops = 0, Hops_old = 0;
    auto t0 = chrono::high_resolution_clock::now();
    bool solved = false;

    while (!solved) {
        // Thread-local counters
        vector<long long> local_hops(omp_get_max_threads(), 0);

        #pragma omp parallel for
        for (int k = 0; k < (Nt + 2 * Nw); ++k) {
            int thread_id = omp_get_thread_num();
            ++local_hops[thread_id];

            if (k < Nt) {
                // Process T[k]
                mpz_class pw = T[k][0] % hop_modulo;
                mpz_class result;
                mpz_fdiv_r(result.get_mpz_t(), T[k][0].get_mpz_t(), DP_rarity.get_mpz_t());
                if (result == 0) {
                    // Check if T[k] exists in W or Wi
                    bool found = false;
                    #pragma omp critical
                    {
                        T.push_back(T[k]);
                        t.push_back(t[k]);
                    }
                    std::set<Point> T_set(T.begin(), T.end());

                    for (size_t i = 0; i < W.size(); ++i) {
                        if (T_set.find(W[i]) != T_set.end()) {
                            #pragma omp critical
                            {
                                auto it = find(T.begin(), T.end(), W[i]);
                                if (it != T.end()) {
                                    int index = distance(T.begin(), it);
                                    mpz_class tP = t[index];
                                    mpz_class wP = w[i];
                                    mpz_class dec = abs(tP - wP);
                                    handle_solution(dec);
                                    solved = true;
                                    found = true;
                                }
                            }
                        }
                        if (found) break;
                    }

                    if (!found) {
                        for (size_t i = 0; i < Wi.size(); ++i) {
                            if (T_set.find(Wi[i]) != T_set.end()) {
                                #pragma omp critical
                                {
                                    auto it = find(T.begin(), T.end(), Wi[i]);
                                    if (it != T.end()) {
                                        int index = distance(T.begin(), it);
                                        mpz_class tP = t[index];
                                        mpz_class wiP = wi[i];
                                        mpz_class dec = abs(tP - wiP);
                                        handle_solution(dec);
                                        solved = true;
                                        found = true;
                                    }
                                }
                            }
                            if (found) break;
                        }
                    }
                }
                if (solved) continue;
                t[k] += powers_of_two[pw.get_ui()];
                T[k] = add(P[pw.get_ui()], T[k], modulo);
            } else if (k < Nt + Nw) {
                // Process W[k - Nt]
                int n = k - Nt;
                mpz_class pw = W[n][0] % hop_modulo;
                mpz_class result;
                mpz_fdiv_r(result.get_mpz_t(), W[n][0].get_mpz_t(), DP_rarity.get_mpz_t());
                if (result == 0) {
                    bool found = false;
                    #pragma omp critical
                    {
                        W.push_back(W[n]);
                        w.push_back(w[n]);
                    }
                    std::set<Point> W_set(W.begin(), W.end());

                    for (size_t i = 0; i < T.size(); ++i) {
                        if (W_set.find(T[i]) != W_set.end()) {
                            #pragma omp critical
                            {
                                auto it = find(W.begin(), W.end(), T[i]);
                                if (it != W.end()) {
                                    int index = distance(W.begin(), it);
                                    mpz_class tP = t[i];
                                    mpz_class wP = w[index];
                                    mpz_class dec = abs(tP - wP);
                                    handle_solution(dec);
                                    solved = true;
                                    found = true;
                                }
                            }
                        }
                        if (found) break;
                    }
                }
                if (solved) continue;
                w[n] += powers_of_two[pw.get_ui()];
                W[n] = add(P[pw.get_ui()], W[n], modulo);
            } else {
                // Process Wi[k - Nt - Nw]
                int n = k - Nt - Nw;
                mpz_class pw = Wi[n][0] % hop_modulo;
                mpz_class result;
                mpz_fdiv_r(result.get_mpz_t(), Wi[n][0].get_mpz_t(), DP_rarity.get_mpz_t());
                if (result == 0) {
                    bool found = false;
                    #pragma omp critical
                    {
                        Wi.push_back(Wi[n]);
                        wi.push_back(wi[n]);
                    }
                    std::set<Point> Wi_set(Wi.begin(), Wi.end());

                    for (size_t i = 0; i < T.size(); ++i) {
                        if (Wi_set.find(T[i]) != Wi_set.end()) {
                            #pragma omp critical
                            {
                                auto it = find(Wi.begin(), Wi.end(), T[i]);
                                if (it != Wi.end()) {
                                    int index = distance(Wi.begin(), it);
                                    mpz_class tP = t[i];
                                    mpz_class wiP = wi[index];
                                    mpz_class dec = abs(tP - wiP);
                                    handle_solution(dec);
                                    solved = true;
                                    found = true;
                                }
                            }
                        }
                        if (found) break;
                    }
                }
                if (solved) continue;
                wi[n] += powers_of_two[pw.get_ui()];
                Wi[n] = add(P[pw.get_ui()], Wi[n], modulo);
            }
        }

        // Aggregate thread-local counters into the global counter
        for (int i = 0; i < local_hops.size(); ++i) {
            Hops += local_hops[i];
        }

        // Update progress
        auto t1 = chrono::high_resolution_clock::now();
        double elapsed_seconds = chrono::duration_cast<chrono::duration<double>>(t1 - t0).count();
        if (elapsed_seconds > 1.0) {
            double hops_per_second = static_cast<double>(Hops - Hops_old) / elapsed_seconds;
            auto elapsed_time = chrono::duration_cast<chrono::seconds>(t1 - starttime);
            int hours = chrono::duration_cast<chrono::hours>(elapsed_time).count();
            int minutes = chrono::duration_cast<chrono::minutes>(elapsed_time % chrono::hours(1)).count();
            int seconds = (elapsed_time % chrono::minutes(1)).count();
            stringstream elapsed_time_str;
            elapsed_time_str << setfill('0') << setw(2) << hours << ":"
                             << setfill('0') << setw(2) << minutes << ":"
                             << setfill('0') << setw(2) << seconds;
            double p_2 = log2(Hops);
            cout << "\r[+] [Hops: 2^" << fixed << setprecision(2) << p_2 << " <-> " << fixed << setprecision(0)
                 << hops_per_second << " h/s] ["
                 << elapsed_time_str.str() << "]" << flush;
            t0 = t1;
            Hops_old = Hops;
        }
    }

    cout << "\r[+] Hops: " << Hops << endl;
    auto end = chrono::high_resolution_clock::now();
    double elapsed_seconds = chrono::duration_cast<chrono::duration<double>>(end - starttime).count();
    return "\r[+] Solution time: " + to_string(elapsed_seconds) + " sec";
}

int main() {
    int puzzle = 50;
    string compressed_public_key = "03f46f41027bbf44fafd6b059091b900dad41e6845b2241dc3254c7cdd3c5a16c6";
    int kangaroo_power = 7;
    mpz_class lower_range_limit = mpz_class(1) << (puzzle - 1);
    mpz_class upper_range_limit = (mpz_class(1) << puzzle) - 1;

    mpz_class DP_rarity = mpz_class(1) << ((puzzle - 2 * kangaroo_power) / 2 - 2);
    int hop_modulo = ((puzzle - 1) / 2) + kangaroo_power;

    int Nt = 1 << kangaroo_power;
    int Nw = 1 << kangaroo_power;

    vector<mpz_class> powers_of_two = generate_powers_of_two(hop_modulo);

    mpz_class X, Y;
    if (compressed_public_key.length() == 66) {
        X = mpz_class(compressed_public_key.substr(2), 16);
        Y = X2Y(X, stoi(compressed_public_key.substr(0, 2)) - 2);
    } else {
        cout << "[error] pubkey len(66/130) invalid!" << endl;
        return 1;
    }

    Point W0 = {X, Y};
    auto starttime = chrono::high_resolution_clock::now();
    time_t currentTime = time(nullptr);
    cout << "\r\033[01;33m[+]\033[32m KANGAROO: \033[01;33m" << ctime(&currentTime) << "\033[0m" << "\r";
    cout << "[+] [Puzzle]: " << puzzle << endl;
    cout << "[+] [Lower range limit]: " << lower_range_limit << endl;
    cout << "[+] [Upper range limit]: " << upper_range_limit << endl;
    cout << "[+] [EC Point Coordinate X]: " << X << endl;
    cout << "[+] [EC Point Coordinate Y]: " << Y << endl;
    mpz_class expected_hops = 2.2 * sqrt(mpz_class(1) << (puzzle - 1));
    double log_expected_hops = log2(expected_hops.get_d());
    cout << "[+] [Expected Hops: 2^" << fixed << setprecision(2)
            << log_expected_hops << " (" << expected_hops << ")]" << endl;

    mpz_class start = lower_range_limit;  
    mpz_class end = upper_range_limit;    
    Point start_point = mul(start);
    Point end_point = mul(end);
    mpz_class inverse_find = start + end;
    Point inverse_find_point = add(start_point, end_point);
    Point W1 = point_subtraction(inverse_find_point, W0);

    vector<Point> P = {PG};
    P.reserve(256);
    for (int k = 0; k < 255; ++k) {
        P.push_back(add(P[k], P[k]));
    }

    unsigned long seed = static_cast<unsigned long>(time(nullptr));
    gmp_randclass rand(gmp_randinit_default);
    rand.seed(seed);

    kangaroo(P, W0, W1, DP_rarity, Nw, Nt, hop_modulo, upper_range_limit, lower_range_limit, powers_of_two);
 
    cout << "\r[+] Average time to solve: " <<
chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - starttime).count() << " sec" << endl;

    return 0;
}