#include <stdio.h>
#include <math.h>
#include <time.h>
#define E 0.001
#define N 512
#define w 1.99

double present[N+1][N+1];
double grid[N+1][N+1];
double X[N+1];
double past[N+1][N+1];

double f(double a, double b) {
	double fr = 0;
	fr = 4 + 2 * pow(a,2) - 2 * a + 2 * pow(b,2) - 2 * b;
	return fr;
}

double max(double a, double b) {
	double m = 0;
	if (a > b)
		m = a;
	else
		m = b;
	return m;
}

double u(double a, double b) {
	double c = 0;
	c = (pow(a,2) - a + 1) * (pow(b,2) - b + 1);
	return c;
}

int main(void) {
	time_t start, end;
	double h = 1.0 / N;
	int iter = 0;

	start = time(NULL);
	for (int i = 0; i < N+1; i++) {
		X[i] = h * i;
	}

	for (int i = 0; i < N + 1; i++) {
		for (int j = 0; j < N + 1; j++) {
				grid[i][j] = u(X[j], X[i]);
			if (i == 0 || j == 0 || i == N || j == N) {
					present[i][j] = grid[i][j];
					past[i][j] = present[i][j];
			}
		}
	}

	double maxe;


	//3 метод
/*	do {
                int i;
                int j;
		maxe = 0;
		for (i = 1; i < N; i++) {
			if (i % 2 == 1) {
				for (j = 1; j < N; j += 2) {
					present[i][j] = (w / 4.0) * (present[i - 1][j] + past[i + 1][j] + present[i][j - 1] + past[i][j + 1] - pow(h, 2) * f(X[j], X[i])) + (1 - w) * present[i][j];
					maxe = max(maxe, fabs(past[i][j] - grid[i][j]));
					past[i][j] = present[i][j];
			}
               }
			if (i % 2 == 0) {
				for (j = 2; j < N; j += 2) {
					present[i][j] = (w / 4.0) * (present[i - 1][j] + past[i + 1][j] + present[i][j - 1] + past[i][j + 1] - pow(h, 2) * f(X[j], X[i])) + (1 - w) * present[i][j];
					maxe = max(maxe, fabs(past[i][j] - grid[i][j]));
					past[i][j] = present[i][j];
			}
		}
        }
		for (i = 1; i < N; i++) {
			if (i % 2 == 1) {
				for (j = 2; j < N; j += 2) {
					present[i][j] = (w / 4.0) * (present[i - 1][j] + past[i + 1][j] + present[i][j - 1] + past[i][j + 1] - pow(h, 2) * f(X[j], X[i])) + (1 - w) * present[i][j];
					maxe = max(maxe, fabs(past[i][j] - grid[i][j]));
					past[i][j] = present[i][j];
			}
                }
			if (i % 2 == 0) {
				for (j = 1; j < N; j += 2) {
					present[i][j] = (w / 4.0) * (present[i - 1][j] + past[i + 1][j] + present[i][j - 1] + past[i][j + 1] - pow(h, 2) * f(X[j], X[i])) + (1 - w) * present[i][j];
					maxe = max(maxe, fabs(past[i][j] - grid[i][j]));
					past[i][j] = present[i][j];
			}
		}
        }
		iter++;
		printf("%d %8.7f\n", iter, maxe);
	} while (maxe >= E);
	end = time(NULL);*/

	// 2 метод
/*		do {
		maxe = 0;
		for (int i = 1; i < N; i++) {
			if (i % 2 == 1) {
				for (int j = 1; j < N; j += 2) {
					present[i][j] = (present[i - 1][j] + past[i + 1][j] + present[i][j - 1] + past[i][j + 1] - pow(h, 2) * f(X[j], X[i])) / 4.0;
					maxe = max(maxe, fabs(past[i][j] - grid[i][j]));
					past[i][j] = present[i][j];
				}
			}
			if (i % 2 == 0) {
				for (int j = 2; j < N; j += 2) {
					present[i][j] = (present[i - 1][j] + past[i + 1][j] + present[i][j - 1] + past[i][j + 1] - pow(h, 2) * f(X[j], X[i])) / 4.0;
					maxe = max(maxe, fabs(past[i][j] - grid[i][j]));
					past[i][j] = present[i][j];
				}
			}
		}
		for (int i = 1; i < N; i++) {
			if (i % 2 == 1) {
				for (int j = 2; j < N; j += 2) {
					present[i][j] = (present[i - 1][j] + past[i + 1][j] + present[i][j - 1] + past[i][j + 1] - pow(h, 2) * f(X[j], X[i])) / 4.0;
					maxe = max(maxe, fabs(past[i][j] - grid[i][j]));
					past[i][j] = present[i][j];
				}
			}
			if (i % 2 == 0) {
				for (int j = 1; j < N; j += 2) {
					present[i][j] = (present[i - 1][j] + past[i + 1][j] + present[i][j - 1] + past[i][j + 1] - pow(h, 2) * f(X[j], X[i])) / 4.0;
					maxe = max(maxe, fabs(past[i][j] - grid[i][j]));
					past[i][j] = present[i][j];
				}
			}
		}
		iter++;
		printf("%d %8.7f\n", iter, maxe);
	} while (maxe >= E);
	end = time(NULL);*/

	// 1 метод
	do {
		maxe = 0;
		for (int i = 1; i < N; i++) {
			for (int j = 1; j < N; j++) {
				present[i][j] = (past[i - 1][j] + past[i + 1][j] + past[i][j - 1] + past[i][j + 1] - pow(h, 2) * f(X[j], X[i])) / 4.0;
				maxe = max(maxe, fabs(past[i][j] - grid[i][j]));
			}
		}
		for (int i = 1; i < N; i++) {
			for (int j = 1; j < N; j++) {
				past[i][j] = present[i][j];
			}
		}
		iter++;
		printf("%d %8.7f\n", iter, maxe);
	} while (maxe >= E);
	end = time(NULL);

/*	FILE* fout;
	fout = fopen("Out.txt", "w");
	for (int i = 1; i < N+1; i++)
		for (int j = 1; j < N+1; j++) {
			fprintf(fout, "%5.4f\t %5.4f\t %5.4f\t\n", X[j], X[i], present[i][j]);
		}
	fclose(fout);*/
	printf("Number of iterations:\n");
	printf("%d\n", iter);
	printf("Time, s:\n");
	printf("%4.2f\n", difftime(end, start));
	return 0;
}