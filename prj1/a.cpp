#include <cstdio>
#include <chrono>
#include <random>
#include <algorithm>
#include <functional>

const long G = 1000000000;
const int MMAX = 5000;
int m;
int s[MMAX], e[MMAX], x[MMAX];
std::default_random_engine gen(std::chrono::system_clock::now().time_since_epoch().count());
std::default_random_engine gens[8];
std::uniform_int_distribution<long> longdist;
std::uniform_int_distribution<int> intdist;

template <int N>
struct longbits {
    static int n;
    static int nbody;
    static long tailmask;

    long d[N];

    static void init(int _n) {
        n = _n;
        nbody = (n - 1) / (sizeof(long) * 8);
        tailmask = (1l << ((n - 1) % (sizeof(long) * 8) + 1)) - 1;
    }

    void set(int i) {
        d[i / (sizeof(long) * 8)] |= 1l << (i % (sizeof(long) * 8));
    }

    void flip(int i) {
        d[i / (sizeof(long) * 8)] ^= 1l << (i % (sizeof(long) * 8));
    }

    bool get(int i) const {
        return d[i / (sizeof(long) * 8)] & (1l << (i % (sizeof(long) * 8)));
    }

    void rand() {
        for (int i = 0; i < N; ++i) {
            d[i] = longdist(gens[i]);
        }
    }

    longbits crossover(const longbits& rhs) const {
        longbits res;
        for (int i = 0; i < N; ++i) {
            res.d[i] = ((d[i] ^ rhs.d[i]) & longdist(gens[i])) | (d[i] & rhs.d[i]);
        }
        return res;
    }

    void mutation(int k) {
        //*
        for (int i = 0; i < k; ++i) {
            flip(intdist(gen) % n);
        }
        //*/
        /*
        for (int i = 0; i < n; ++i) {
            if (intdist(gen) % n < 1) {
                flip(i);
            }
        }
        //*/
    }

    int evaluate() const {
        bool a[N * sizeof(long) * 8];
        long t = 0;
        for (int i = 0; i < n; ++i) {
            if (i % (sizeof(long) * 8) == 0) {
                t = d[i / (sizeof(long) * 8)];
            }
            a[i] = t & 1;
            t >>= 1;
        }
        int w = 0;
        for (int i = 0; i < m; ++i) {
            if (a[s[i]] != a[e[i]]) {
                w += x[i];
            }
        }
        return w;
    }

    int dist(const longbits& rhs) const {
        int s = 0;
        for (int i = 0; i < nbody; ++i) {
            s += __builtin_popcount(d[i] ^ rhs.d[i]);
        }
        s += __builtin_popcount((d[nbody] ^ rhs.d[nbody]) & tailmask);
        return s;
    }

    void printans(FILE *fout) const {
        for (int j = 0; j < n; ++j) {
            if (!get(j)) {
                fprintf(fout, "%d ", j + 1);
            }
        }
        fprintf(fout, "\n");
    }
};

template <int N> int longbits<N>::n;
template <int N> int longbits<N>::nbody;
template <int N> long longbits<N>::tailmask;

typedef longbits<8> b512;

const int NMAX = 500;
const int NCUR = 256; // 256 is best for L1$
//const int NCUR = 2048; // 2048 is best for L2$
const int NNEW = 1024;
const int METAGEN = 8;

b512 *curgen, *curgen_another, curgen_buf0[NCUR + NNEW], curgen_buf1[NCUR + NNEW];
int *curgen_score, *curgen_score_another, curgen_score_buf0[NCUR + NNEW], curgen_score_buf1[NCUR + NNEW];
int curgen_index[NCUR + NNEW];
long curgen_fit[NCUR + NNEW];
b512 bestgen[NCUR];

void reportcurgen() {
    int t = 5, q = 10;
    printf("=== curgen ===\n");
    for (int i = 0; i < t; ++i) {
        printf("[T%d] %d\n", i, curgen_score[i]);
    }
    for (int i = 1; i < q; ++i) {
        printf("[Q%d] %d\n", i, curgen_score[NCUR / q * i]);
    }
}

void sortcurgen() {
    for (int i = 0; i < NCUR; ++i) {
        curgen_index[i] = i;
    }
    std::sort(curgen_index, curgen_index + NCUR, [](int x, int y) { return curgen_score[x] > curgen_score[y]; });
    for (int i = 0; i < NCUR; ++i) {
        curgen_another[i] = curgen[curgen_index[i]];
        curgen_score_another[i] = curgen_score[curgen_index[i]];
    }
    std::swap(curgen, curgen_another);
    std::swap(curgen_score, curgen_score_another);
}

long gettime() {
    return std::chrono::system_clock::now().time_since_epoch().count();
}

int main(int argc, char *argv[]) {
    long st = gettime(), et;
    FILE *fin = fopen(argv[1], "r");
    FILE *fout = fopen(argv[2], "w");
    for (int i = 0; i < 8; ++i) {
        gens[i] = std::default_random_engine(gettime() + i);
    }

    int n;
    fscanf(fin, "%d%d", &n, &m);
    for (int i = 0; i < m; ++i) {
        fscanf(fin, "%d%d%d", &s[i], &e[i], &x[i]);
        --s[i], --e[i];
    }
    b512::init(n);
    curgen = curgen_buf0, curgen_another = curgen_buf1;
    curgen_score = curgen_score_buf0, curgen_score_another = curgen_score_buf1;

    for (int i = 0; i < NCUR; ++i) {
        long s = 0;
        s = NCUR - i;
        if (i == 0) {
            curgen_fit[i] = s;
        } else {
            curgen_fit[i] = curgen_fit[i - 1] + s;
        }
    }
    for (int gg = 0; gg <= METAGEN; ++gg) {
        // make gen0
        for (int i = 0; i < NCUR; ++i) {
            if (gg == METAGEN) {
                curgen[i] = bestgen[i];
            } else {
                curgen[i].rand();
            }
            curgen_score[i] = curgen[i].evaluate();
        }
        sortcurgen(); //printf("GEN -1\n"); reportcurgen();

        long convst = gettime(), convet;
        for (int g = 0; ; ++g) {
            // generate newgen
            for (int i = 0; i < NCUR; ++i) {
                curgen_another[i] = curgen[i];
                curgen_score_another[i] = curgen_score[i];
            }
            for (int i = 0; i < NNEW; ++i) {
                int j = std::upper_bound(curgen_fit, curgen_fit + NCUR, longdist(gen) % curgen_fit[NCUR - 1]) - curgen_fit;
                int k = std::upper_bound(curgen_fit, curgen_fit + NCUR, longdist(gen) % curgen_fit[NCUR - 1]) - curgen_fit;
                b512 chr = curgen[j].crossover(curgen[k]);
                chr.mutation(intdist(gen) % (n / 100 + 2));
                int chr_score = chr.evaluate();
                int mindist = n + 1, minl = -1;
                for (int l = 0; l < NCUR; ++l) {
                    int dist = chr.dist(curgen[l]);
                    if (mindist > dist) {
                        mindist = dist;
                        minl = l;
                    }
                }
                if (curgen_score_another[minl] < chr_score) {
                    curgen_score_another[minl] = chr_score;
                    curgen_another[minl] = chr;
                }
            }
            std::swap(curgen, curgen_another);
            std::swap(curgen_score, curgen_score_another);
            sortcurgen();
            if (gg < METAGEN) {
                convet = gettime();
                if ((convet - convst) / G >= 8) {
                    for (int i = 0; i < NCUR / METAGEN; ++i) {
                        bestgen[NCUR / METAGEN * gg + i] = curgen[i];
                    }
                    break;
                }
            }
            et = gettime();
            if ((et - st) / G >= 175) {
                break;
            }
        }
        et = gettime();
        if ((et - st) / G >= 175) {
            break;
        }
    }
    curgen[0].printans(fout);
    //fprintf(fout, "%d\n", curgen[0].evaluate());

    return 0;
}
