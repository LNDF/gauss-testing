/*
Features:
    - Determinant
    - Gauss
    - Determinant (Gauss)
    - Range (Gauss)
    - Inverse matrix (Gauss)
    - Gauss-Jordan
*/


#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#define GAUSS_INVERSE        1
#define GAUSS_DIAG           2
#define GAUSS_PIVOT          4
#define GAUSS_IGNORE_LAST    8

typedef uint8_t u8; //integer 8-bit
typedef uint16_t u16; //integer 16-bit
typedef uint32_t u32; //integer 32-bit

u8 steps = 0; //step-by-step

//Estructura de matriz
typedef struct matrix {
    u32 width;
    u32 height;
    double* mem;
    double** rows;
} matrix;

/*

  Allocate matrix in memory
  - mtx: pointer to matrix struct
  - width: matrix width
  - height: matrix height

*/
void allocMatrix(matrix* mtx, u32 width, u32 height) {
    double *start, *end;
    u32 i = 0;
    mtx->mem = (double*) malloc(sizeof(double) * width * height);
    mtx->rows = (double**) malloc(sizeof(double*) * width);
    mtx->width = width;
    mtx->height = height;
    memset(mtx->mem, 0, sizeof(double) * width * height);
    start = mtx->mem;
    end = start + (width * height);
    for (; start < end; start += width, ++i) {
        mtx->rows[i] = start;
    }
}

/*

  Allocate identity matrix in memory
  - mtx: pointer to matrix struct
  - size: width and height

*/
void allocIdentityMatrix(matrix* mtx, u32 size) {
    u32 i = 0;
    allocMatrix(mtx, size, size);
    for (; i < size; ++i) {
        mtx->rows[i][i] = 1;
    }
}

/*

  Copy a matrix from mtx to dest.
  - mtx: pointer to a matrix struct
  - dest: pointer to another matrix struct

*/
void copyMatrix(matrix* source, matrix* dest) {
    allocMatrix(dest, source->width, source->height);
    memcpy(dest->mem, source->mem, sizeof(double) * source->width * source->height);
}

/*

  Copy a matrix from mtx to submtx, removing a row and a column.
  - mtx: pointer to a matrix struct
  - dest: pointer to another matrix struct
  - col: column to remove
  - row: row to remove

*/
void allocSubmatrix(matrix* mtx, matrix* submtx, u32 col, u32 row) {
    u32 i, c = 0, posPart2 = col + 1;
    size_t sizePart1 = sizeof(double) * col, sizePart2 = sizeof(double) * (mtx->width - col - 1);
    allocMatrix(submtx, mtx->width - 1, mtx->height - 1);
    for (i = 0; i < mtx->height; ++i) {
        if (i == row) continue;
        memcpy(submtx->rows[c], mtx->rows[i], sizePart1);
        memcpy(submtx->rows[c] + col, mtx->rows[i] + posPart2, sizePart2);
        ++c;
    }
}


/*

  Swap two rows of a matrix.
  - mtx: pointer to a matrix struct
  - r1: the row to swap
  - r2: the other row to swap with

*/
void swapRowsMatrix(matrix* mtx, u32 r1, u32 r2) {
    double* swap = mtx->rows[r1];
    mtx->rows[r1] = mtx->rows[r2];
    mtx->rows[r2] = swap;
}


/*

  Free memory used by a matrix struct.
  - mtx: pointer to a matrix struct

*/
void freeMatrix(matrix* mtx) {
    free(mtx->mem);
    free(mtx->rows);
}

/*

  Print the matrix
  - mtx: pointer to a matrix struct

*/
void printMatrix(matrix* mtx) {
    u32 i, j;
    for (i = 0; i < mtx->height; ++i) {
        for (j = 0; j < mtx->width; ++j) {
            printf("%lf  ", mtx->rows[i][j]);
        }
        printf("\n");
    }
}

void printGaussStep(matrix* gauss, matrix* inverse, u8 isInverse) {
    if (steps) {
        printf("Gauss step:\n");
        printMatrix(gauss);
        if (isInverse) {
            printf("\nInverse step:\n");
            printMatrix(inverse);
        }
        printf("\n");
    }
}

/*

  Gauss alorithm
  - mtx: pointer to a matrix struct
  - out: pointer to a matrix struct. This struct will contain the output.
  - flags:
      - GAUSS_INVERSE: calc inverse matrix (used with GAUSS_DIAG | GAUSS_PIVOTES_AL_1)
      - GAUSS_DIAG: also make zeros the numbers avobe the matrix disg
      - GAUSS_PIVOT: make all pivots 1
      - GAUSS_IGNORE_LAST: ignore last column

    For Gauss-Jordan use GAUSS_DIAG | GAUSS_PIVOT | GAUSS_IGNORE_LAST

    Returnse: 1 if the sign for the determinant should be changed, 0 otherwise

*/
u8 gaussMatrix(matrix* mtx, matrix* out, u8 flags) {
    u32 i, j, k, columns, inverted = 0;
    double multiplier, pivotNumber;
    matrix tmp;
    u8 doInverse = flags & GAUSS_INVERSE;
    u8 pivotsOne = flags & GAUSS_PIVOT;
    u8 doDiag = flags & GAUSS_DIAG;
    u8 ignoreLast = flags & GAUSS_IGNORE_LAST;
    columns = ignoreLast ? mtx->width - 1 : mtx->width;
    if (columns > mtx->height) columns = mtx->height;
    if (doInverse) {
        if (mtx->width != mtx->height) return 2;
        allocIdentityMatrix(out, mtx->width);
        copyMatrix(mtx, &tmp);
    } else {
        copyMatrix(mtx, out);
        memcpy(&tmp, out, sizeof(matrix));
    }
    for (i = 0; i < columns; ++i) {
        if (tmp.rows[i][i] == 0) {
            for (j = i + 1; j < tmp.height; ++j) {
                if (tmp.rows[j][i] != 0) {
                    swapRowsMatrix(&tmp, i, j);
                    if (doInverse == 1) swapRowsMatrix(out, i, j);
                    inverted = !inverted;
                    printGaussStep(&tmp, out, doInverse);
                    break;
                }
            }
        }
        if (tmp.rows[i][i] == 0) continue;
        pivotNumber = tmp.rows[i][i];
        if (pivotsOne) {
            for (j = 0; j < tmp.width; ++j) {
                tmp.rows[i][j] /= pivotNumber;
                if (doInverse) out->rows[i][j] /= pivotNumber;
            }
            printGaussStep(&tmp, out, doInverse);
            pivotNumber = 1;
        }
        for (j = doDiag ? 0 : i + 1; j < tmp.height; ++j) {
            if ((tmp.rows[j][i] == 0 && !doInverse) || j == i) continue;
            multiplier = tmp.rows[j][i] / pivotNumber;
            for (k = doInverse ? 0 : i; k < tmp.width; ++k) {
                tmp.rows[j][k] -= tmp.rows[i][k] * multiplier;
                if (doInverse) out->rows[j][k] -= out->rows[i][k] * multiplier;
            }
        }
        printGaussStep(&tmp, out, doInverse);
    }
    if (doInverse) {
        freeMatrix(&tmp);
    }
    return inverted;
}

/*

  Calc matrix determinant
  - mtx: pointer to a matrix struct

  Returns: Determinant of the matrix

*/
double determinant(matrix* mtx) {
    if (mtx->width != mtx->height) return;
    if (mtx->width == 1) return mtx->mem[0];
    u32 i;
    matrix submtx;
    double detFinal = 0, tmp;
    for (i = 0; i < mtx->width; ++i) {
        allocSubmatrix(mtx, &submtx, i, 0);
        tmp = mtx->rows[0][i] * determinant(&submtx);;
        if (i % 2 != 0) tmp *= -1;
        detFinal += tmp;
        freeMatrix(&submtx);
    }
    return detFinal;
}

/*

  Calc matrix determinant using Gauss
  - mtx: pointer to a matrix struct

  Returns: Determinant of the matrix

*/
double determinantGauss(matrix* mtx) {
    matrix gaussTmp;
    u32 i;
    double finalDet = 1;
    if (mtx->height != mtx->width) return;
    u32 inverted = gaussMatrix(mtx, &gaussTmp, 0);
    for (i = 0; i < mtx->width; ++i) {
        finalDet *= gaussTmp.rows[i][i];
    }
    freeMatrix(&gaussTmp);
    return inverted == 1 ? -finalDet : finalDet;
}

/*

  Calc matrix range using Gauss
  - mtx: pointer to a matrix struct

  Devuelve: Range of the matrix

*/
u32 range(matrix* mtx) {
    matrix gaussTmp;
    u32 r = 0, i, j;
    gaussMatrix(mtx, &gaussTmp, 0);
    for (i = 0; i < gaussTmp.height; ++i) {
        for (j = 0; j < gaussTmp.width; ++j) {
            if (gaussTmp.rows[i][j] != 0) {
                ++r;
                break;
            }
        }
    }
    freeMatrix(&gaussTmp);
    return r;
}

/*

  Prints solutions for the previusly Gauss-Jordaned matrix
  - mtx: pointer to a matrix struct

*/
void printGaussJordan(matrix* mtx) {
    u32 i, j, ew = mtx->width - 1, params = 0, param = 0, p = 0;
    u8 zeros;
    for (i = 0; i < ew; ++i, p = 0) {
        for (j = i + 1; j < ew; ++j) {
            if (mtx->rows[i][j] != 0.0) ++p;
        }
        if (p > params) params = p;
    }
    for (i = 0; i < mtx->height; ++i) {
        if (mtx->rows[i][ew] != 0.0) {
            zeros = 1;
            for (j = 0; j < ew; ++j) {
                if (mtx->rows[i][j] != 0.0) {
                    zeros = 0;
                    break;
                }
            }
            if (zeros) {
                printf("There is no solution.\n");
                return;
            }
        }
    }
    for (i = 0; i < ew - params; ++i) {
        printf("x%d=%lf", i + 1, mtx->rows[i][ew]);
        for (j = i + 1; j < ew; ++j) {
            if (mtx->rows[i][j] != 0) {
                if (mtx->rows[i][j] < 0) {
                    printf("+");
                }
                printf("%lfx%d", -mtx->rows[i][j], j + 1);
            }
        }
        printf("\n");
    }
    for (; i < ew; ++i) {
        printf("x%d=k%d\n", i + 1, ++param);
    }
}

/*

  The entrypoint

  Returns: 0

*/
int main() {
    matrix mtx, gauss;
    u32 width, height, rg, i, j, op = 0;
    double res;
    printf("Enter the width of the matrix: ");
    scanf("%u", &width);
    printf("Enter the height of the matrix: ");
    scanf("%u", &height);
    allocMatrix(&mtx, width, height);
    for (i = 0; i < mtx.height; ++i) {
        for (j = 0; j < mtx.width; ++j) {
            printf("Enter the number for the position [%ux%u]: ", i + 1, j + 1);
            scanf("%lf", &mtx.rows[i][j]);
        }
    }
    printf("\n");
    printMatrix(&mtx);
    printf("\n");
    while (op != 8) {
        printf("\n");
        printf("Select an option:\n 1. Determinant\n 2. Gauss\n 3. Determinant (Gauss)\n 4. Range\n 5. Inverse\n 6. Solve as equation\n 7. Switch step-by-step flag (currently %u)\n 8. Quit\nEnter a number: ", steps);
        scanf("%u", &op);
        printf("\n");
        switch (op) {
            case 1:
            case 3:
                if (mtx.width != mtx.height) {
                    printf("This is not an sqare matrix.\n");
                    break;
                }
                if (op == 1) {
                    res = determinant(&mtx);
                } else {
                    res = determinantGauss(&mtx);
                }
                printf("Determinant of the matrix: %lf\n", res);
                break;
            case 2:
                gaussMatrix(&mtx, &gauss, 0);
                printf("Gauss:\n");
                printMatrix(&gauss);
                freeMatrix(&gauss);
                break;
            case 4:
                rg = range(&mtx);
                printf("The range is %u\n", rg);
                break;
            case 5:
                if (mtx.width != mtx.height) {
                    printf("The matrix should be an sqare matrix.\n");
                    break;
                }
                gaussMatrix(&mtx, &gauss, GAUSS_INVERSE | GAUSS_PIVOT | GAUSS_DIAG);
                printf("Inverse matrix:\n");
                printMatrix(&gauss);
                freeMatrix(&gauss);
                break;
            case 6:
                if (mtx.width <= 1) {
                    printf("The matrix should have at least two columns.\n");
                    break;
                }
                gaussMatrix(&mtx, &gauss, GAUSS_IGNORE_LAST | GAUSS_DIAG | GAUSS_PIVOT);
                printf("Gauss-Jordan:\n");
                printMatrix(&gauss);
                printf("\nSolution:\n");
                printGaussJordan(&gauss);
                freeMatrix(&gauss);
                break;
            case 7:
                steps = !steps;
            default:
                break;
        }
    }
    freeMatrix(&mtx);
    return 0;
}
