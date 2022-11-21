#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 1599          // xの要素数
#define p 11            // xのパラメータ数
#define C 6             // クラスタ数
#define max_iter 10000  //最大試行回数
#define theta 2.0       //ファジィ度
#define epsilon 0.001

double calc_dist(double a[][p], double b[][p], int j, int c) {
  double dist = 0;
  for (int i = 0; i < p; i++) {
    dist += (a[j][i] - b[c][i]) * (a[j][i] - b[c][i]);
  }
  return dist;
}

double calc_bc(double x[][p], double u[][N], int P, int c) {
  double bunsi = 0;
  double bunbo = 0;

  for (int i = 0; i < N; i++) {
    bunbo += pow(u[c][i], theta);
    bunsi += pow(u[c][i], theta) * x[i][P];
  }
  return bunsi / bunbo;
}

double calc_centre(double x[][p], double u[][N], int P, int c) {
  double bunsi = 0;
  double bunbo = 0;

  for (int i = 0; i < N; i++) {
    bunbo += u[c][i];
    bunsi += u[c][i] * x[i][P];
  }
  return bunsi / bunbo;
}

double calc_uci(double x[][p], double b[][p], int i, int c) {
  double temp = 0;
  double temp_bunsi = calc_dist(x, b, i, c);

  for (int l = 0; l < C; l++) {
    temp += pow((temp_bunsi / calc_dist(x, b, i, l)), (1 / (theta - 1)));
  }
  return 1 / temp;
}

double calc_J_fcm(double u[][N], double x[][p], double b[][p]) {
  double temp = 0;

  for (int c = 0; c < C; c++) {
    for (int i = 0; i < N; i++) {
      temp += pow(u[c][i], theta) * calc_dist(x, b, i, c);
    }
  }
  return temp;
}

double calc_J_km(double u[][N], double x[][p], double b[][p]) {
  double temp = 0;

  for (int c = 0; c < C; c++) {
    for (int i = 0; i < N; i++) {
      temp += calc_dist(x, b, i, c);
    }
  }
  return temp;
}

int main(void) {
  /*read data*/
  double x[N][p];
  FILE *fp;
  char fname[] = "winequality-red.csv";
  char str[16];
  double f[p + 1];

  fp = fopen(fname, "r");  // ファイルを開く。失敗するとNULLを返す。
  if (fp == NULL) {
    printf("%s file not open!\n", fname);
    return -1;
  }

  int count = 0;
  while (fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &f[0],
                &f[1], &f[2], &f[3], &f[4], &f[5], &f[6], &f[7], &f[8], &f[9],
                &f[10], &f[11]) != EOF) {
    for (int i = 0; i < p; i++) {
      x[count][i] = f[i];
    }
    count++;
  }

  fclose(fp);  // ファイルを閉じる

  for (int iter = 0; iter < 1; iter++) {
    double _j = 1e+8;
    double u[C][N] = {0};
    double k[C][N] = {0};
    double b[C][p] = {0};
    double b_km[C][p] = {0};
    int k_res[N] = {0};

    /*initialize u,k*/
    srand((unsigned int)time(NULL));
    for (int i = 0; i < C; i++) {
      for (int j = 0; j < N; j++) {
        u[i][j] = rand() / (double)RAND_MAX;
        k[i][j] = u[i][j];
      }
    }

    /*Fuzzy c-means*/

    for (int t = 0; t < max_iter; t++) {
      for (int c = 0; c < C; c++) {
        for (int P = 0; P < p; P++) {
          b[c][P] = calc_bc(x, u, P, c);
        }
      }

      for (int c = 0; c < C; c++) {
        for (int i = 0; i < N; i++) {
          u[c][i] = calc_uci(x, b, i, c);
        }
      }
      double conv = calc_J_fcm(u, x, b);
      if (fabs(_j - conv) < epsilon) {
        break;
      } else {
        _j = conv;
      }
    }

    /*k-means*/
    _j = 1e+8;
    for (int t = 0; t < max_iter; t++) {
      for (int c = 0; c < C; c++) {
        for (int P = 0; P < p; P++) {
          b_km[c][P] = calc_centre(x, k, P, c);
        }
      }

      /*ｋの更新式*/
      for (int i = 0; i < N; i++) {  //すべてのｘに対して
        int assigned = -1;
        double temp_dist = 1e+8;
        for (int c = 0; c < C; c++) {
          if (temp_dist > calc_dist(x, b_km, i, c)) {
            assigned = c;
            temp_dist = calc_dist(x, b_km, i, c);
          }
        }
        for (int c = 0; c < C; c++) {
          if (c == assigned) {
            k[c][i] = 1;
          } else {
            k[c][i] = 0;
          }
        }
      }

      /*収束判定*/
      double convergence = calc_J_km(k, x, b_km);
      if (fabs(_j - convergence) < epsilon) {
        // printf("%d\n", t);
        break;
      } else {
        _j = convergence;
      }
    }

    /*結果 u*/
    printf("/*結果 u*/\n");
    for (int i = 0; i < C; i++) {
      for (int j = 0; j < N; j++) {
        printf("%.3lf,", u[i][j]);
      }
      printf("\n");
    }
    printf("\n\n\n\n");

    printf("/*結果 b クラスタ中心*/\n");
    for (int i = 0; i < C; i++) {
      for (int j = 0; j < p; j++) {
        printf("%.2lf, ", b[i][j]);
      }
      printf("\n");
    }

    /*結果 k*/
    // printf("/*結果 k*/\n");
    /*for (int i = 0; i < N; i++) {
      for (int j = 0; j < C; j++) {
        if (k[j][i] == 1) {
          k_res[i] = j;
        }
      }
      printf("%d,", k_res[i]);
    }*/

    printf("\n");
  }
  return 0;
}