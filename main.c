main.c
O
Тип
С
Размер
5 КБ (5 073 байта)
Использовано
5 КБ (5 073 байта)
Расположение
chmi_2sem
Владелец
я
Изменено
мной 15:06
Открыто
мной 18:13
Создано
15:06 в приложении ZIP Extractor
Добавить описание
Читателям разрешено скачивать файл

#include "stdio.h"
#include <stdlib.h>
#include "math.h"
#include "main.h"

double V_nplus1_m(double v_n_mplus1, double v_n_m){
    return (tao/(2*h))*(v_n_mplus1 - v_n_m) + v_n_m;
}


int stepx_by_x(double x){
    return (int)(fabs(x - x_begin)/h);
}

int stept_by_t(double t){
    return (int)(fabs(t - t_begin)/tao);
}

double t_by_stept(int step){
    return (double)(t_begin + step*tao);
}

double x_by_stepx(int step){
    return (double)(x_begin + step*h);
}

void print_block(double * block, int length){
    for (int i = 0; i < length; ++i)
    {
        printf("%-4.5g ", block[i]);
    }
    printf("\n");
}


double delta(double v_n_m, double x, double t){
    double delta;
    if (x<=-0.5*t){
        delta  = fabs(v_n_m - 1);
        // printf("delta %g\n", delta);
    }
    else{
        delta = fabs(v_n_m);
        // printf("delta %g\n", delta);
    }
    return delta;
}


void norma_c(double delta, double *max){
    if (fabs(delta) >  *max) *max = fabs(delta);
}


void norma_l(double v_n_m, double * norma){
    *norma += fabs(v_n_m);
}


void print_matrix(double * matrix, int length){
    printf("\nMATRIX PRINT\n");
    for (int i = 0; i < length; i++)
        for (int j = 0; j < length; j++){
            printf("%-4.5g ", *(matrix + i*length + j));
    }
    printf("\n");
}


void print_block_to_file(double * block, int length, FILE * fout){
    for (int i = 0; i < length; ++i)
    {
        fprintf(fout, "%-8g ", block[i]);
    }
    fprintf(fout, "\n");
}

void begin_value(double * block, int length){
    for (int i = 0; i < length; ++i)
    {
        if (x_by_stepx(i) <= 0)
            block[i] = 1.0;
        else
            block[i] = 0.0;
    }
}

void characteristic(double * block){
    block[0] = 1;
    block[m-1] = 0;
}

void next_block(double * block, int length, double step_by_t,
    double * norma_DELTA_c_line, double * norma_delta_c_line, //[c]
    double * norma_DELTA_l_line, double * norma_delta_l_line  //[l]
    ){

    double x = x_begin;
    double del;
    *norma_DELTA_c_line = 0; *norma_delta_c_line = 0; //начальное максимальное значение [c]
    *norma_DELTA_l_line = 0; *norma_delta_l_line = 0; //начальное значение [l]

    characteristic(block);
    for (int i = 0; i < m-1; ++i)
    {
        block[i] = V_nplus1_m(block[i+1], block[i]);

        del = delta(block[i], x, step_by_t*tao); // разность между численным решением и точным решением
        norma_c(del, norma_DELTA_c_line);        //[c]
        norma_c(block[i], norma_delta_c_line);   //[c]
        norma_l(del, norma_DELTA_l_line);        //[l]
        norma_l(block[i], norma_delta_l_line);

        x+=h;
    }
        // printf("NORMA 1 %g\n", *norma_DELTA_c_line);
        // printf("NORMA 2 %g\n", *norma_delta_c_line);

}

void main_loop_lin_simple(){
    // double previous_block[m];
    double block[m];
    double norma_DELTA_c = 0; double norma_delta_c = 0;   //по всей сетке [c]
    double norma_DELTA_c_line; double norma_delta_c_line; //по линиям
    double norma_DELTA_l = 0; double norma_delta_l = 0;   //по всей сетке [l]
    double norma_DELTA_l_line; double norma_delta_l_line; //по линиям
    FILE * fout = fopen("output_1.txt", "w+");

    int length = (int)(sizeof(block)/sizeof(double));
    begin_value(block, length);

    for (int i = 0; i < n; ++i){
        next_block(block, length, i, &norma_DELTA_c_line, &norma_delta_c_line, &norma_DELTA_l_line, &norma_delta_l_line);
        norma_c(norma_DELTA_c_line, &norma_DELTA_c);
        norma_c(norma_delta_c_line, &norma_delta_c);
        norma_l(norma_DELTA_l_line, &norma_DELTA_l);
        norma_l(norma_delta_l_line, &norma_delta_l);
        // printf("NORMA_delta_l %g\n", norma_delta_l_line);

        // print_block(block, length);
        print_block_to_file(block, length, fout);

    }

    printf("NORMA_DELTA_c %g\n", norma_DELTA_c);
    printf("NORMA_delta_c %g\n", norma_DELTA_c/norma_delta_c);
    printf("NORMA_DELTA_l %g\n", h*norma_DELTA_l);
    printf("NORMA_delta_l %g\n", norma_DELTA_l/norma_delta_l);

}


/////////////////////
void main_loop_lin_dif(){

    double * matrix; // указатель на массив
    matrix = (double *)malloc(m*m * sizeof(double));

    FILE * fout = fopen("output_2.txt", "w+");

    int length = (int)(sizeof(matrix)/sizeof(double));
    // print("length %g", )
    for (int i = 0; i < length; i++)
        for (int k = 0; k <= 2; k++)
            if ((i+k) < m) {
                *(matrix + i*m + k) = (double)(rand() % 10);
                // printf();
            }

    print_matrix(matrix, length);


    // begin_value(block, length);

    // for (int i = 0; i < n; ++i){
    //     next_block(block, length);
    //     // print_block(block, length);
    //     print_block_to_file(block, length, fout);
    // }
}



///////////////////////

int main(){
    printf("%d x %d\n", m,n);
    main_loop_lin_simple();
    // main_loop_lin_dif();
    return 0;
}
