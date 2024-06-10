

#include <stdio.h>
#include <math.h>
#include <string.h>




void save_to_csv(const char *filename, double *t, char* t_name, double *y, char* y_name, int size) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Failed to open file");
        return;
    }
    fprintf(file, t_name);
    fprintf(file, ",");
    fprintf(file, y_name);
    fprintf(file, "\n");
    for (int i = 0; i < size; ++i)
    {
        fprintf(file, "%.6e,%.6e\n", t[i], y[i]);
    }
    fclose(file);
}
void matrixVectorMultiply(double matrix[2][2], double vector[2], double result[2]) {
    result[0] = matrix[0][0] * vector[0] + matrix[0][1] * vector[1];
    result[1] = matrix[1][0] * vector[0] + matrix[1][1] * vector[1];
}
void vectorMatrixMultiply(double vector[2], double matrix[2][2], double result[2]) {
    result[0] = vector[0] * matrix[0][0] + vector[1] * matrix[1][0];
    result[1] = vector[0] * matrix[0][1] + vector[1] * matrix[1][1];
}
void vectorVectorMultiply(double row_vector[2], double column_vector[2], double *result) {
    *result = row_vector[0] * column_vector[0] + row_vector[1] * column_vector[1];
}
void invertMatrix(double A[2][2], double A_inv[2][2]) {
    double determinant = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    
    if (determinant == 0) {
        printf("Error: Matrix is singular, cannot be inverted.\n");
        return;
    }
    
    double inv_det = 1.0 / determinant;
    
    A_inv[0][0] = A[1][1] * inv_det;
    A_inv[0][1] = -A[0][1] * inv_det;
    A_inv[1][0] = -A[1][0] * inv_det;
    A_inv[1][1] = A[0][0] * inv_det;
}
void changeDelimiter(const char *original_delim, const char *new_delim, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) 
    {
        printf("Error: Unable to open file %s\n", filename);
        return;
    }

    FILE *tempFile = tmpfile();
    if (tempFile == NULL) 
    {
        printf("Error: Unable to create temporary file\n");
        fclose(file);
        return;
    }

    char line[1000];
    while (fgets(line, sizeof(line), file) != NULL) 
    {
        char *delim = strstr(line, original_delim);
        while (delim != NULL) 
        {
            memcpy(delim, new_delim, strlen(new_delim));
            delim = strstr(delim + strlen(new_delim), original_delim);
        }
        fputs(line, tempFile);
    }

    fclose(file);
    rewind(tempFile);

    file = fopen(filename, "w");
    if (file == NULL) 
    {
        printf("Error: Unable to open file %s\n", filename);
        fclose(tempFile);
        return;
    }

    int c;
    while ((c = fgetc(tempFile)) != EOF)
    {
        fputc(c, file);
    }

    fclose(tempFile);
    fclose(file);
}




int main() 
{
    printf("\n");
    printf(" _____________________________________________ \n");
    printf("|                                             |\n");
    printf("|          MECHATRONIKA HAZI FELADAT          |\n");
    printf("|             DC MOTOR MODELLEZES             |\n");
    printf("|_____________________________________________|\n");
    printf("\n\n\n\n\n");





    // ##################### ADATOK ###################

    // Aramkor parameterek
    double u_n = 42;        // [V]
    double i_n = 0.659;     // [A]
    double R = 16.6;        // [ohm]
    double L = 2.22e-3;     // [H]
    // Motor parameterek
    double B = 2e-6;        // [Nms/rad]
    double k_m = 70.4e-3;   // [Nm/A]
    double J_r = 43.8e-7;   // [kgm2]
    // Bojgomu parameterek
    double N = 35;             // [-]
    double J_t = 200e-7;    // [kgm2]
    // Egyeb
    double dt = 0.050;      // [s]

    // Adatok kiiratasa
    printf("========== ADATOK ==========\n");
    printf("\n\tAramkor\n");
    printf("\t\tR  = %E [ohm]\n", R);
    printf("\t\tL  = %E [H]\n", L);
    printf("\tMotor\n");
    printf("\t\tB  = %E [Nms/rad]\n", B);
    printf("\t\tJr = %E [kgm2]\n", J_r);
    printf("\t\tkm = %E [Nm/A]\n", k_m);
    printf("\tBolygomu\n");
    printf("\t\tN  = %E [-]\n", N);
    printf("\tTarcsa\n");
    printf("\t\tJt = %E [kgm2]\n", J_t);
    printf("\n\n\n\n\n");





    // ############# UO ATVITELI FUGGVENY #############

    /*
        ss * b0 + s * b1 + 1 * b2
        -------------------------
        ss * a0 + s * a1 + 1 * a2
    */

    // Atviteli fuggveny parameterek
    double b0 = k_m * N;
    double a2 = J_t * L + N * N * J_r * L;
    double a1 = J_t * R + N * N * J_r * R + N * N * B * L;
    double a0 = k_m * k_m * N * N + N * N * B * R;
    printf("========== UO ATVITELI FUGGVENY ==========\n");
    printf("\n\tW_UO(s) = ____ss(%E)_+_s(%E)_+_(%E)____\n", 0.0, 0.0, b0);
    printf(  "\t              ss(%E) + s(%E) + (%E)\n", a2, a1, a0);

    // Diszkriminans
    double D = a1 * a1 - 4 * a0 * a2;

    // Polusok
    double p1 = (-a1 + sqrt(D)) / (2 * a2);
    double p2 = (-a1 - sqrt(D)) / (2 * a2);
    printf("\n\tA rendszer polusai:\n");
    printf("\t\tp1 = %E\n", p1);
    printf("\t\tp2 = %E\n", p2);

    // Atviteli fuggveny egyutthatoi
    double v1 = b0 / ((p1 - p2) * p1);
    double v2 = -b0 / ((p1 - p2) * p2);
    double v3 = b0 / (p1 * p2);
    printf("\n\tAz atviteli fuggveny alakja:\n");
    printf("\n\t\tw_uo(t) = %E * exp(%E * t) + %E * exp(%E * t) + %E\n", v1, p1, v2, p2, v3);

    // Idoallandok
    double T1 = 1 / fabs(p1);
    double T2 = 1 / fabs(p2);
    printf("\n\tAz atviteli fuggveny idoallandoi:\n");
    printf("\t\tT1 = %E [s]\n\t\tT2 = %E [s]\n", T1, T2);
    printf("\n\n\n\n\n");
    




    // ########### ATVITELI FUGGVENYEK ALLANDOSULT ALLAPOTAI ###########

    double W_UO_Stat = k_m * N / a0;
    double W_MO_Stat = N * N * B / a0;
    double W_UI_Stat = -R / a0;
    double W_MI_Stat = k_m * N / a0;
    printf("========== ATVITELI FUGGVENYEK ALLANDOSULT ALLAPOTAI ==========\n");
    printf("\n\tlim(s->0)  W_UO(s) = %E [rad/s]\n\tlim(s->0)  W_UI(s) = %E [A]\n\tlim(s->0)  W_MO(s) = %E [rad/s]\n\tlim(s->0)  W_MI(s) = %E [A]\n",  W_UO_Stat, W_MO_Stat, W_UI_Stat, W_MI_Stat);
    printf("\n\n\n\n\n");





    // ########### MAXIMALIS TERHELO NYOMATEK ###########

    double M_t_max = (i_n - W_MO_Stat * u_n) / W_MI_Stat;
    printf("========== MAXIMALIS TERHELO NYOMATEK ==========\n");
    printf("\n\tM_t_max = %E [Nm]\n", M_t_max);
    printf("\n\n\n\n\n");





    // ########### ALLAPOTTER MODELL ###########
    
    // Elore tarto Euler
    double A_fe[2][2] = {
        {1 - (dt * R / L), -dt * k_m / L},
        {k_m * dt / J_r, 1 - (B * dt) / J_r}
    };
    double B_fe[2] = {dt / L, 0};
    double C_fe[2] = {0, 1};
    double D_fe = 0;
    printf("========== ALLAPOTTER MODELL ==========\n");
    printf("\n\tElore tarto (forward) Euler-formulaval\n");
    printf("\t\t     _                                   _ \n");
    printf("\t\tA = |   %14.6E   %14.6E   |\n", A_fe[0][0], A_fe[0][1]);
    printf("\t\t    |_  %14.6E   %14.6E  _|\n", A_fe[1][0], A_fe[1][1]);
    printf("\t\t     _                  _ \n");
    printf("\t\tB = |   %14.6E   |\n", B_fe[0]);
    printf("\t\t    |_  %14.6E  _|\n", B_fe[1]);
    printf("\n");
    printf("\t\tC = [   %14.6E   %14.6E   ]\n", C_fe[0], C_fe[1]);
    printf("\n");
    printf("\t\tD = [   %14.6E   ]\n", D_fe);

    // Hatra tarto Euler (be ~ backward Euler)
    double E_be[2][2] = {{1,0}, {0,1}};
    double A_be_inv[2][2];
    double A_be[2][2];
    double B_be[2];
    double C_be[2];
    double D_be;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            A_be_inv[i][j] = 2 * E_be[i][j] - A_fe[i][j];
        }
    }
    invertMatrix(A_be_inv, A_be);               
    matrixVectorMultiply(A_be, B_fe, B_be);
    vectorMatrixMultiply(C_fe, A_be, C_be);
    vectorVectorMultiply(C_fe, B_be, &D_be);
    printf("\n\tHatra tarto (backward) Euler-formulaval\n");
    printf("\t\t     _                                   _ \n");
    printf("\t\tA = |   %14.6E   %14.6E   |\n", A_be[0][0], A_be[0][1]);
    printf("\t\t    |_  %14.6E   %14.6E  _|\n", A_be[1][0], A_be[1][1]);
    printf("\t\t     _                  _ \n");
    printf("\t\tB = |   %14.6E   |\n", B_be[0]);
    printf("\t\t    |_  %14.6E  _|\n", B_be[1]);
    printf("\n");
    printf("\t\tC = [   %14.6E   %14.6E   ]\n", C_be[0], C_be[1]);
    printf("\n");
    printf("\t\tD = [   %14.6E   ]\n", D_be);
    printf("\n\n\n\n\n");





    // ########### GRAFIKON ADATOK KIMENTESE ###########

    // w_uo atmeneti fuggveny
    int data_points = 1e3;
    double t[data_points];
    double w_uo[data_points];
    for (int i = 0; i < data_points; ++i) {
        t[i] = i * 1e-4;
        w_uo[i] = v1 * exp(p1 * t[i]) + v2 * exp(p2 * t[i]) + v3;
    }
    printf("========== ADATOK MENTESE ==========\n");
    save_to_csv("w_uo_atmeneti_fgv.csv", t, "t",w_uo, "w_uo", data_points);
    printf("\n\tA w_uo atmeneti fuggveny el van mentve!\n\t\tw_uo_atmeneti_fgv.csv\n");

    // Eloretarto (forward) Euler
    double y[data_points];
    double u[data_points];
    for (int i = 0; i < data_points; ++i) {
        u[i] = 42;
    }
    double X_1[2] = {0, 0};
    double X_0[2] = {0, 0};
    for (int k = 0; k < data_points; ++k) {
        X_1[0] = A_fe[0][0] * X_0[0] + A_fe[0][1] * X_0[1] + B_fe[0] * u[k];
        X_1[1] = A_fe[1][0] * X_0[0] + A_fe[1][1] * X_0[1] + B_fe[1] * u[k];
        y[k] = C_fe[0] * X_0[0] + C_fe[1] * X_0[1] + D_fe * u[k];
        t[k] = k*dt;
        X_0[0] = X_1[0];
        X_0[1] = X_1[1];
    }
    save_to_csv("forward_euler.csv", t, "t [s]", y, "O [rad/s]", data_points);
    printf("\n\tA w_uo atmeneti fuggveny el van mentve!\n\t\tforward_euler.csv\n");

    // Hatratarto (backward) Euler
    for (int i = 0; i < 2; i++)
    {
        X_0[i] = 0;
        X_1[i] = 0;
    }
    
    for (int k = 0; k < data_points; ++k)
    {
        X_1[0] = A_be[0][0] * X_0[0] + A_be[0][1] * X_0[1] + B_be[0] * u[k];
        X_1[1] = A_be[1][0] * X_0[0] + A_be[1][1] * X_0[1] + B_be[1] * u[k];
        y[k] = C_be[0] * X_0[0] + C_be[1] * X_0[1] + D_be * u[k];
        t[k] = k*dt;
        X_0[0] = X_1[0];
        X_0[1] = X_1[1];
    }
    save_to_csv("backward_euler.csv", t, "t [s]", y, "O [rad/s]", data_points);
    printf("\n\tA w_uo atmeneti fuggveny el van mentve!\n\t\tbackward_euler.csv\n");
    printf("\n\n\n\n\n");
    
    // CSV MAGYARRA ALLITASA, MERT NEM KEPES RA AZ EXCEL
    changeDelimiter(",",";", "w_uo_atmeneti_fgv.csv");
    changeDelimiter(".",",", "w_uo_atmeneti_fgv.csv");
    changeDelimiter(",",";", "forward_euler.csv");
    changeDelimiter(".",",", "forward_euler.csv");
    changeDelimiter(",",";", "backward_euler.csv");
    changeDelimiter(".",",", "backward_euler.csv");





    // ########### PARAMETERVALTOZAS HATASA ###########
    double tmp = B;
    B = B*10;
    a0 = k_m * k_m * N * N + N * N * B * R;
    W_UO_Stat = k_m * N / a0;
    W_MO_Stat = N * N * B / a0;
    W_UI_Stat = -R / a0;
    W_MI_Stat = k_m * N / a0;
    printf("========== B PARAMETER VALTOZASANAK HATASA ==========\n");
    printf("\nB * 0.1 = 2E-7 [Nms/rad] eseten:\n");
    M_t_max = (i_n - W_MO_Stat * u_n) / W_MI_Stat;
    printf("\n\tM_t_max = %E [Nm]\n", M_t_max);
    B = tmp;
    a0 = k_m * k_m * N * N + N * N * B * R;
    W_UO_Stat = k_m * N / a0;
    W_MO_Stat = N * N * B / a0;
    W_UI_Stat = -R / a0;
    W_MI_Stat = k_m * N / a0;
    printf("\n");

    printf("\nB * 1 = 2E-6 [Nms/rad] eseten:\n");
    M_t_max = (i_n - W_MO_Stat * u_n) / W_MI_Stat;
    printf("\n\tM_t_max = %E [Nm]\n", M_t_max);
    B = tmp/10;
    a0 = k_m * k_m * N * N + N * N * B * R;
    W_UO_Stat = k_m * N / a0;
    W_MO_Stat = N * N * B / a0;
    W_UI_Stat = -R / a0;
    W_MI_Stat = k_m * N / a0;
    printf("\n");
    
    printf("\nB * 10 = 2E-5 [Nms/rad] eseten:\n");
    M_t_max = (i_n - W_MO_Stat * u_n) / W_MI_Stat;
    printf("\n\tM_t_max = %E [Nm]\n", M_t_max);
    printf("\n\n\n\n\n");





    // ########### VEGE ###########
    printf("\n");
    printf(" ______________________________________________ \n");
    printf("|                                              |\n");
    printf("|             DOBONDI-REISZ MENDEL             |\n");
    printf("|                    C8GJJW                    |\n");
    printf("|                 2024.05.16.                  |\n");
    printf("|______________________________________________|\n");
    printf("\n\n\n\n\n");

    return 0;
}

