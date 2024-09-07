#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define n 4 // number of bits




//                          UTILITIES

#define pow(i) (1 << (i))                               // 2^i
#define llpow(i) ((long long int)1 << (i))              // 2^i
#define N pow(n)                                        // N = 2^n
#define bit(x,i) (((x) >> (i)) % 2)                     // read i-th bit of x
#define node(d, t) ((t)*pow(n-(d)) + pow(n-(d)-1) - 1)  // linear index of a node in binary tree (in-order depth-first, starting from zero)
#define jthfather(d, t, j) (node((d)-j, (t)/pow(j)))    // linear index of the j-th father of a node in binary tree

void print_function(int* f) {
    for(int i = 0; i < N; i++) printf("%d", f[i]);
    printf("\n");
}

void print_tree(int* tree) {
    for(int i = 0; i < N - 1; i++) printf("%c", 97 + tree[i]);  // x0=a, x1=b, x2=c, x3=d, ...
    printf("\n");
}

int products_in_DNF(int* g) {
    int count = 0; 
    for(long long unsigned int i = 0; i < N; i++) { count += g[i]; }
    return (n-1)*count;
}

int products_in_ANF(int* g) {
    int count = 0;
    for(long long unsigned int i = 0; i < N; i++) { 
        int terms = 0;
        for(int j = 0; j < n; j++) terms += bit(i, j);
        count += terms > 1 ? g[i] * (terms - 1) : 0; 
    }
    return count;
}




//                          MOBIUS TRANSFORM

void mobius_transform(int* g) {
    int f[N];
    for(int i = 0; i < n; i++) {    
        for(int x = 0; x < N; x++) {
            if( bit(x,i) ) {
                f[x] = g[x] ^ g[x - pow(i)];
            }else{
                f[x] = g[x];
            }
        }
        for(int x = 0; x < N; x++) g[x] = f[x];
    }
}




//                          GREEDY ALGORITHM

int find_variable_to_collect(int g[], int t, int variables) {
    const int length = pow(variables);
    const int* table = &g[t*length];

    int best_v = 0;
    int best_count = 0;
    for(int v = 0; v < variables; v++) {
        int count_v = 0;
        for(int i = 0; i < pow(variables-v-1); i++) {
            for(int j = 0; j < pow(v); j++) {
                count_v += table[pow(v)*(2*i+1) + j];
            }
        }
        if(count_v > best_count) {
            best_count = count_v;
            best_v = v;
        }
    } 

    return best_v;
}

int collect(int f[], int g[], int t, int m, int v) {
    const int length = pow(m);
    int* table_f = &f[t*length];
    int* table_g = &g[t*length];

    int no_products = 1;
    for(int i = 0; i < pow(m-v-1); i++) {
        for(int j = 0; j < pow(v); j++) {
            table_f[pow(v)*i + j] = table_g[pow(v)*2*i + j];
            table_f[length/2 + pow(v)*i + j] = table_g[pow(v)*(2*i+1) + j];
            if(!(i == 0 && j == 0) && table_g[pow(v)*(2*i+1) + j]) no_products = 0; 
        }
    }

    return !no_products;
}

int minimize_products(int* g, int* tree) {
    if(tree) memset(tree, 0, (N-1)*sizeof(int));
    
    int products = 0;
    
    int f[N];
    int all_variables[n]; 
    for(int step = 0; step < n; step++) {    
        int m = n - step;
        for(int t = 0; t < pow(step); t++) {
            // Collect variable
            int variable_to_collect = find_variable_to_collect(g, t, m);
            products += collect(f, g, t, m, variable_to_collect);

            // Save in the tree which of the n original variables we have collected
            if(tree) {
                for(int i = 0; i < n; i++) all_variables[i] = 1;
                for(int j = 1; j <= step; j++) all_variables[tree[jthfather(step, t, j)]] = 0;
                int count = 0;
                for(int i = 0; i < n; i++) {
                    if(count == variable_to_collect && all_variables[i] == 1){
                        tree[node(step, t)] = i;
                        break;
                    }
                    count += all_variables[i];
                }
            }
        }
        for(int x = 0; x < N; x++) g[x] = f[x]; // could use a memcpy
    }

    return products;
}




//                          COLLECT IN ORDER

int minimize_products_collect_in_order(int* g, int* tree) {
    if(tree) memset(tree, 0, (N-1)*sizeof(int));
    
    int products = 0;
    
    int f[N];
    int all_variables[n]; 
    for(int step = 0; step < n; step++) {    
        int m = n - step;
        for(int t = 0; t < pow(step); t++) {
            // Collect variable
            int variable_to_collect = 0;
            products += collect(f, g, t, m, variable_to_collect);

            // Save in the tree which of the n original variables we have collected
            if(tree) {
                for(int i = 0; i < n; i++) all_variables[i] = 1;
                for(int j = 1; j <= step; j++) all_variables[tree[jthfather(step, t, j)]] = 0;
                int count = 0;
                for(int i = 0; i < n; i++) {
                    if(count == variable_to_collect && all_variables[i] == 1){
                        tree[node(step, t)] = i;
                        break;
                    }
                    count += all_variables[i];
                }
            }
        }
        for(int x = 0; x < N; x++) g[x] = f[x]; // could use a memcpy
    }

    return products;
}




//                          OUTPUT A READABLE EQUATION

typedef struct {
    int constant;       // -1 not a constant, 0 constant zero, 1 constant one
    int terms;          // number of summed terms in the expression
    char* string;       // string representation if constant = -1, NULL otherwise
    int length;         // length of the string, \0 excluded but required in string
} equation;
equation human_readable_equation(int* leaves, int* tree, int depth, int t){
    // XOR(+), AND(juxtaposition)
    // x0=a, x1=b, x2=c, x3=d, ...
    char variable = 97 + tree[node(depth, t)]; 
    equation outside;
    equation inside;
    if(depth == n - 1) {
        outside = (equation) { leaves[2*t + 0], 1, NULL, 0};
        inside = (equation) { leaves[2*t + 1], 1, NULL, 0};
    }else{
        outside = human_readable_equation(leaves, tree, depth + 1, 2*t + 0);
        inside  = human_readable_equation(leaves, tree, depth + 1, 2*t + 1);
    }
    if(outside.constant == 0 && inside.constant == 0) {
        return (equation) { 0, 0, NULL, 0 };
    }else if(outside.constant == 0 && inside.constant == 1) {
        char* string = malloc(2 * sizeof(char)); 
        string[0] = variable; string[1] = '\0';
        return (equation) { -1, 1, string, 1 };
    }else if(outside.constant == 1 && inside.constant == 0) {
        return (equation) { 1, 1, NULL, 0 };
    }else if(outside.constant == 1 && inside.constant == 1) {
        char* string = malloc(4 * sizeof(char));
        string[0] = '1'; string[1] = '+'; string[2] = variable; string[3] = '\0';
        return (equation) { -1, 2, string, 3 };
    }else if(outside.constant == 0 && inside.constant == -1) {
        if(inside.terms > 1) {
            char* string = malloc((inside.length+4) * sizeof(char));
            sprintf(string, "%c(%s)", variable, inside.string);
            return (equation) { -1, 1, string, inside.length+3 };
        }else{
            char* string = malloc((inside.length+2) * sizeof(char));
            sprintf(string, "%c%s", variable, inside.string);
            return (equation) { -1, 1, string, inside.length+1 };
        }
    }else if(outside.constant == -1 && inside.constant == 0) {
        return outside;
    }else if(outside.constant == 1 && inside.constant == -1) {
        if(inside.terms > 1) {
            char* string = malloc((inside.length+6) * sizeof(char));
            sprintf(string, "1+%c(%s)", variable, inside.string);
            return (equation) { -1, 2, string, inside.length+5 };
        }else{
            char* string = malloc((inside.length+4) * sizeof(char));
            sprintf(string, "1+%c%s", variable, inside.string);
            return (equation) { -1, 2, string, inside.length+3 };
        }
    }else if(outside.constant == -1 && inside.constant == 1) {
        char* string = malloc((outside.length+3) * sizeof(char));
        sprintf(string, "%s+%c", outside.string, variable);
        return (equation) { -1, outside.terms + 1, string, outside.length+2 };
    }else if(outside.constant == -1 && inside.constant == -1) {
        if(inside.terms > 1) {
            char* string = malloc((outside.length+inside.length+6) * sizeof(char));
            sprintf(string, "%s+%c(%s)", outside.string, variable, inside.string);
            return (equation) { -1, outside.terms + inside.terms, string, outside.length+inside.length+5 };
        }else{
            char* string = malloc((outside.length+inside.length+3) * sizeof(char));
            sprintf(string, "%s+%c%s", outside.string, variable, inside.string);
            return (equation) { -1, outside.terms + inside.terms, string, outside.length+inside.length+2 };
        }
    }else{
        exit(1);
    }
}




//                          BENCHMARKS

void sample_all_functions(int histograms[][4]) {
    if(n > 4) exit(1); // llpow doesn't work anymore 

    // Try them all
    int g[N];

    for(long long int x = 0; x < llpow(N); x++) {
        // DNF corresponding to binary representation of x
        for(long long int j = 0; j < N; j++) g[j] = (int) bit(x, j);
        histograms[products_in_DNF(g)][2] += 1;
        
        // Calculate ANF
        mobius_transform(g);
        histograms[products_in_ANF(g)][3] += 1;

        // Minimize products with greedy approach
        int* f = malloc(N*sizeof(int)); memcpy(f, g, N*sizeof(int));
        int products_greedy = minimize_products(f, NULL);
        histograms[products_greedy][0] += 1;

        // Minimize products blindly collecting the first variable every time
        int products_in_order = minimize_products_collect_in_order(g, NULL);
        histograms[products_in_order][1] += 1;

        if(x % 500000 == 0) printf("\rDone: %.5f%c  ", ((double)x) / llpow(N), '%'); fflush(stdout); 
    }
    printf("\rDone: 100%c     \n", '%'); fflush(stdout); 
}

void random_truth_table(int* g) {
    // sizeof(int) * 8 > 16 (always, usually = 32)
    // => 15 random bits (1 bit goes for the sign)

    for(int i = 0; i < N; i+=15) {
        unsigned int r = rand();
        for(int ii = 0; ii < 15; ii++) {
            if(i + ii < N) g[i+ii] = bit(r,ii);
        }
    }
}

void sample_functions_randomly(int histograms[][4], int samples) {
    int g[N];
    for(int x = 0; x < samples; x++) {
        // Random DNF 
        // We do not sample ANFs directly only to check the statistics of the DNF
        random_truth_table(g);
        histograms[products_in_DNF(g)][2] += 1;
        
        // Calculate ANF
        mobius_transform(g);
        histograms[products_in_ANF(g)][3] += 1;

        // Minimize products with greedy approach
        int* f = malloc(N*sizeof(int)); memcpy(f, g, N*sizeof(int));
        int products_greedy = minimize_products(f, NULL);
        histograms[products_greedy][0] += 1;

        // Minimize products blindly collecting the first variable every time
        int products_in_order = minimize_products_collect_in_order(g, NULL);
        histograms[products_in_order][1] += 1;

        // Log progress
        if(x % (samples / 1000) == 0) printf("\rDone: %.2f%c  ", ((double)x) / samples * 100.0, '%'); fflush(stdout); 
    }
    printf("\rDone: 100%c (%d samples) \n", '%', samples); fflush(stdout); 
}

void save_histogram(int histograms[][4], int BINS, FILE* file) {
    for(int i = 0; i < BINS; i++) {
        //                                Products  Greedy            CollectInOrder    DNF               ANF
        fprintf(file, "%d %d %d %d %d\n", i,        histograms[i][0], histograms[i][1], histograms[i][2], histograms[i][3]);
    }
}

void histogram_statistics(int histograms[][4], int samples) {
    long long unsigned int temp_greedy = 0;
    long long unsigned int temp_in_order = 0;
    long long unsigned int temp_DNF = 0;
    long long unsigned int temp_ANF = 0;
    for(long long int i = 0; i < n*N; i++) {
        temp_greedy += i * histograms[i][0];
        temp_in_order += i * histograms[i][1];
        temp_DNF += i * histograms[i][2];
        temp_ANF += i * histograms[i][3];
    }
    long double mean_greedy = (long double) temp_greedy / samples;
    long double mean_in_order = (long double) temp_in_order / samples;
    long double mean_DNF = (long double) temp_DNF / samples;
    long double mean_ANF = (long double) temp_ANF / samples;

    long double var_greedy = 0.0;
    long double var_in_order = 0.0;
    long double var_DNF = 0.0;
    long double var_ANF = 0.0;
    for(long long int i = 0; i < n*N; i++) {
        var_greedy += histograms[i][0] * (i - mean_greedy) * (i - mean_greedy) / ((long double) samples);
        var_in_order += histograms[i][1] * (i - mean_in_order) * (i - mean_in_order) / ((long double) samples);
        var_DNF += histograms[i][2] * (i - mean_DNF) * (i - mean_DNF) / ((long double) samples);
        var_ANF += histograms[i][3] * (i - mean_ANF) * (i - mean_ANF) / ((long double) samples);
    }
    printf("mean_greedy: %LF ± %LF \n", mean_greedy, sqrtl(var_greedy / ((long double) samples)) );
    printf("mean_in_order: %LF ± %LF \n", mean_in_order, sqrtl(var_in_order / ((long double) samples)) );
    printf("mean_ANF: %LF ± %LF \n", mean_ANF, sqrtl(var_ANF / ((long double) samples)) );
    printf("mean_DNF: %LF ± %LF \n", mean_DNF, sqrtl(var_DNF / ((long double) samples)) );
}


//                          TEST FUNCTIONS

// g(x, y) = (x > y)
void greater_than_truth_table(int* g) {
    if(n % 2 == 1) { printf("Error: greater_than_truth_table."); exit(1); }
    if(n > 8*sizeof(long long unsigned int)) { printf("Error: greater_than_truth_table."); exit(2); }
    for(long long unsigned int i = 0; i < N; i++) {
        int first = 0;
        int second = 0;
        for(int j = 0; j < n/2; j++) {
            first += pow(j) * bit(i, j);
            second += pow(j) * bit(i, n/2 + j);
        }
        g[i] = first > second;
    }
}

// g(x, y) = (x == y)
void equals_to_truth_table(int* g) {
    if(n % 2 == 1) { printf("Error: equals_to_truth_table."); exit(1); }
    if(n > 8*sizeof(long long unsigned int))  { printf("Error: equals_to_truth_table."); exit(2); }
    for(long long unsigned int i = 0; i < N; i++) {
        int first = 0;
        int second = 0;
        for(int j = 0; j < n/2; j++) {
            first += pow(j) * bit(i, j);
            second += pow(j) * bit(i, n/2 + j);
        }
        g[i] = first == second;
    }
}


// g(x, y) = x^y
void elevation_truth_table(int* g) {
    if(n % 2 == 1) { printf("Error: elevation_truth_table."); exit(1); }
    if(n > 8*sizeof(long long unsigned int)) { printf("Error: elevation_truth_table."); exit(2); }
    for(long long unsigned int i = 0; i < N; i++) {
        g[i] = 1;
        for(int j = 0; j < n/2; j++) {
            if(bit(i, n/2 + j) == 1 && bit(i, j) == 0) {
                g[i] = 0;
                break;
            }
        }
    }
}

int main() {
    printf("n = %d\n", n);



    // TEST A SPECIFIC FUNCTION

    // Input the function truth table manually
    //               000  001  010  011  100  101  110  111
    //               ZYX  ZYx  ZyX  Zyx  zYX  zYx  zyX  zyx        // DNF
    // int g[N] = {    0,   1,   1,   0,   0,   0,   1,   1    }; 
    
    // Use one of the implemented structured functions
    int g[N]; 
    equals_to_truth_table(g);
    // greater_than_truth_table(g);
    // elevation_truth_table(g);

    int products_DNF = products_in_DNF(g);
    printf("DNF: ");print_function(g);

    // Compute the Mobius Transform
    //               000  001  010  011  100  101  110  111
    //                 1    x    y   yx    z   zx   yz  zyx        // ANF
    mobius_transform(g);
    int products_ANF = products_in_ANF(g);
    printf("ANF: ");print_function(g);

    // Minimize products
    int tree[N - 1];
    int products = minimize_products(g, tree);
    printf("Tree: ");print_tree(tree);
    printf("Leaves: ");print_function(g);

    // Output human readable equation
    equation g_equation = human_readable_equation(g, tree, 0, 0);
    printf("Equation: %s\n", g_equation.string);
    printf("Products: (DNF %d) (ANF %d) (GREEDY %d)\n", products_DNF, products_ANF, products );



    // TEST A SAMPLE OF FUNCTIONS

    int BINS = n*N;
    int histograms[BINS][4];
    memset(histograms, 0, 4*BINS*sizeof(int));

    // int samples = pow(N);
    // sample_all_functions(histograms);
    int samples = 5e3;
    sample_functions_randomly(histograms, samples);

    char file_histogram_name[50]; sprintf(file_histogram_name, "data/histogram%d.dat", n);
    FILE* file_histogram = fopen(file_histogram_name, "w+");
    save_histogram(histograms, BINS, file_histogram);

    histogram_statistics(histograms, samples);



    return 0;
}
